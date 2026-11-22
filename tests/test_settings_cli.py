"""Strict CLI/file parsing, exact input paths, and cleanup ownership regressions."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


def main():
    binary = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix='zuck-settings-') as scratch:
        root = Path(scratch)
        cycles = root / 'empty cycles.txt'
        cycles.write_text('')
        output = root / 'result'
        output.mkdir()
        sentinel = output / 'entries.tsv'
        sentinel.write_text('previous result\n')
        config = root / 'settings.txt'

        def run(arguments, success=True):
            result = subprocess.run([str(binary), *arguments], capture_output=True,
                                    text=True, timeout=60)
            assert (result.returncode == 0) == success, result.stdout + result.stderr
            return result

        common = ['--cycles', str(cycles), '--output-folder', str(output), '--threads', '1', '--ram', '1G']
        bad_values = {
            'threads': ['1junk', '-1', '18446744073709551616'],
            'cycle-max-length': ['90junk', '-1', '2147483648'],
            'cycle-min-length': ['2.5', '-1'],
            'threshold-multiplicity': ['20junk', '-1', '18446744073709551616'],
            'ram': ['1Gjunk', '1.2.3G', 'NaN', 'inf', '0G', '1e9'],
        }
        for name, values in bad_values.items():
            for value in values:
                run([*common, '--' + name, value], success=False)
                assert sentinel.read_text() == 'previous result\n'
                config.write_text(f'{name}={value}\n')
                # Do not mask invalid settings with a CLI override for that key.
                run(['--cycles', str(cycles), '--output-folder', str(output),
                     '--settings', str(config)], success=False)
                assert sentinel.read_text() == 'previous result\n'
        for name in ('retain-graph', 'retain-cycles', 'experimental-alternative-views'):
            config.write_text(f'{name}=maybe\n')
            run([*common, '--settings', str(config)], success=False)
        config.write_text('threads=999999\nram=2G\n')
        run([*common, '--settings', str(config)])  # CLI wins irrespective of order.
        run([*common, '--settings', str(config), '--settings', str(config)], success=False)
        run([*common, '--input-files'], success=False)
        run([*common, '--benchmark', 'expected.txt'], success=False)
        run([*common, '--output-folder', '--threads', '1'], success=False)

        # A configuration file may itself be a result-name or hard-link alias.
        settings_result = output / 'arrays.tsv'
        settings_result.write_text('threads=1\nram=1G\n')
        before_settings = settings_result.read_bytes()
        run([*common, '--settings', str(settings_result)], success=False)
        assert settings_result.read_bytes() == before_settings
        alias = root / 'settings-alias'
        alias.hardlink_to(settings_result)
        run([*common, '--settings', str(alias)], success=False)
        assert settings_result.read_bytes() == before_settings

        # Full graph builds must preserve whitespace and '#' in literal paths,
        # including when both the input and graph output directories use spaces.
        reads = root / 'read inputs'
        reads.mkdir()
        source = Path(__file__).resolve().parents[1] / 'libs/megahit/test_data'
        first = reads / 'read #1.fa'
        second = reads / 'read #2.fa'
        shutil.copyfile(source / 'r3_1.fa', first)
        shutil.copyfile(source / 'r3_2.fa', second)
        built = root / 'built result'
        config.write_text(f'input-files="{first}", "{second}" # paired reads\n')
        run(['--settings', str(config), '--output-folder', str(built), '--retain-graph',
             '--threads', '1', '--ram', '1G'])
        assert (built / 'graph/graph.sdbg_info').exists()
        from_cli = root / 'cli result'
        run(['--input-files', str(first), str(second), '--output-folder', str(from_cli),
             '--threads', '1', '--ram', '1G'])
        assert (from_cli / 'entries.tsv').read_bytes() == (built / 'entries.tsv').read_bytes()

        # Supplied graphs survive replay even in a directory formerly removed by cleanup.
        original_graph = built / 'graph'
        graph_files = {p.name: p.read_bytes() for p in original_graph.iterdir()}
        run(['--graph', str(original_graph), '--output-folder', str(built),
             '--threads', '1', '--ram', '1G'])
        assert {p.name: p.read_bytes() for p in original_graph.iterdir()} == graph_files
        alias_output = root / 'graph-alias-output'
        alias_output.mkdir()
        (alias_output / 'cycles.txt').hardlink_to(original_graph / 'graph.sdbg_info')
        run(['--graph', str(original_graph), '--output-folder', str(alias_output),
             '--threads', '1', '--ram', '1G'], success=False)
        assert (original_graph / 'graph.sdbg_info').read_bytes() == graph_files['graph.sdbg_info']
        nested = root / 'nested'
        nested.mkdir()
        shutil.copytree(original_graph, nested / 'cycles')
        run(['--graph', str(nested / 'cycles'), '--output-folder', str(nested),
             '--cleanup-cycles', '--threads', '1', '--ram', '1G'])
        assert (nested / 'cycles/graph.sdbg_info').read_bytes() == graph_files['graph.sdbg_info']

        # Custom graph cleanup removes only generated files and preserves unrelated data.
        custom = root / 'custom graph'
        custom.mkdir()
        (custom / 'notes.txt').write_text('keep me')
        config.write_text(f'graph-folder={custom}\n')
        run(['--input-files', str(first), str(second), '--output-folder', str(from_cli),
             '--settings', str(config), '--threads', '1', '--ram', '1G'])
        assert list(custom.iterdir()) == [custom / 'notes.txt']
        assert (custom / 'notes.txt').read_text() == 'keep me'

        # Reject input aliases before the generated cycle file can truncate them.
        conflicting = root / 'conflict'
        conflicting.mkdir()
        raw = conflicting / 'cycles.txt'
        shutil.copyfile(first, raw)
        before = raw.read_bytes()
        run(['--input-files', str(raw), '--output-folder', str(conflicting),
             '--threads', '1', '--ram', '1G'], success=False)
        assert raw.read_bytes() == before
        graph = conflicting / 'graph'
        graph.mkdir()
        for name in ('data.lib', 'outfile_prefix.bin', 'graph.sdbg.0'):
            for hardlink in (False, True):
                destination = graph / name
                if hardlink:
                    destination.hardlink_to(first)
                else:
                    destination.symlink_to(first)
                original = first.read_bytes()
                run(['--input-files', str(first), '--output-folder', str(conflicting),
                     '--threads', '1', '--ram', '1G'], success=False)
                assert first.read_bytes() == original
                destination.unlink()
    print('Strict settings, exact paths, precedence and cleanup ownership passed')


if __name__ == '__main__':
    main()
