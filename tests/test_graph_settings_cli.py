"""Real SDBG: singleton policy and thread-invariant cycle records and calls."""
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def main():
    binary = Path(sys.argv[1]).resolve()
    rng = random.Random(231113)
    dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    repeat = dna(40)
    spacers = ['ACGT'[i % 4] + dna(31) + 'TGCA'[i % 4] for i in range(9)]
    available = os.cpu_count() or 1
    with tempfile.TemporaryDirectory(prefix='zuck-graph-settings-') as temporary:
        root = Path(temporary)
        reads = root / 'reads.fa'
        sequences = [repeat + spacer + repeat for i, spacer in enumerate(spacers)
                     for _ in range(1 if i == 8 else 12)]
        reads.write_text(''.join(f'>r{i}\n{s}\n' for i, s in enumerate(sequences)))

        def run(name, workers, extra=(), env=None, build=False):
            output = root / name
            source = ['--input-files', str(reads), '--retain-graph'] if build else ['--graph', str(root / 'built/graph')]
            result = subprocess.run([str(binary), *source, '--output-folder', str(output),
                '--threads', str(workers), '--ram', '1G', '--read-linked-groups', 'false',
                '--confidence-scores', 'false', *extra], env={**os.environ, **(env or {})},
                capture_output=True, text=True, timeout=60)
            assert result.returncode == 0, result.stdout + result.stderr
            return output

        built = run('built', min(2, available), build=True)
        # A cycle destination must not follow links into an unrelated file.
        protected = root / 'protected-input.txt'
        protected.write_text('preserve this file\n')
        obstructed = root / 'linked-output'
        obstructed.mkdir()
        (obstructed / 'cycles.txt').symlink_to(protected)
        failed = subprocess.run([str(binary), '--graph', str(built / 'graph'),
            '--output-folder', str(obstructed), '--threads', '1', '--ram', '1G'],
            capture_output=True, text=True, timeout=60)
        assert failed.returncode != 0, failed.stdout + failed.stderr
        assert protected.read_text() == 'preserve this file\n'
        assert (obstructed / 'cycles.txt').is_symlink()
        assert not list(obstructed.glob('.Zuck-cycles-*'))

        rare = spacers[-1]
        reverse = rare.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        has_rare = lambda path: any(rare in c+c or reverse in c+c for c in path.read_text().splitlines())
        assert not has_rare(built / 'cycles.txt')
        preserved = run('preserved', min(2, available), ['--low-abundance', 'true'])
        assert has_rare(preserved / 'cycles.txt'), 'Retaining multiplicity-one edges must recover the planted rare path'
        configurations = [(w, {}) for w in (1, 2, 4, 8) if w <= available]
        configurations += [(min(8, available), {'OMP_THREAD_LIMIT': '1'}),
                           (min(8, available), {'OMP_DYNAMIC': 'TRUE'})]
        for workers, env in configurations:
            output = run('workers-' + str(workers) + '-' + str(len(env)) + '-' + str(bool(env.get('OMP_DYNAMIC'))), workers, env=env)
            for name in ['cycles.txt', 'arrays.tsv', 'entries.tsv']:
                assert (output / name).read_bytes() == (built / name).read_bytes(), (workers, env, name)
        config = root / 'settings.txt'
        config.write_text('low-abundance=true\ntraversal-partitions=8\n')
        from_file = run('file', 1, ['--settings', str(config)])
        assert (from_file / 'cycles.txt').read_bytes() == (preserved / 'cycles.txt').read_bytes()
        overridden = run('override', 1, ['--settings', str(config), '--low-abundance', 'false'])
        assert (overridden / 'cycles.txt').read_bytes() == (built / 'cycles.txt').read_bytes()
        manifest = json.loads((built / 'run_manifest.json').read_text())['cycle_search_settings']
        assert manifest['traversal_partitions'] == 8 and manifest['low_abundance'] is False
        assert manifest['applied_in_this_run'] is True and manifest['singleton_edges_removed'] is True
        config.write_text('low-abundance=maybe\n')
        invalid = subprocess.run([str(binary), '--cycles', str(built / 'cycles.txt'), '--settings', str(config),
                                  '--output-folder', str(root / 'invalid')], capture_output=True, text=True, timeout=30)
        assert invalid.returncode != 0
    print('SDBG singleton retention, settings precedence, worker and OpenMP-team invariance passed')


if __name__ == '__main__':
    main()
