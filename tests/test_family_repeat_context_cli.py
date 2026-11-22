"""Isolate mixed repeat-copy support using replay; discovery is held fixed."""
import csv
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    binary = Path(sys.argv[1]).resolve()
    rng = random.Random(781425)
    dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    repeat = dna(36)
    alternate = repeat[:5] + ('A' if repeat[5] != 'A' else 'C') + repeat[6:]
    spacers = ['ACGT'[i % 4] + dna(30) + 'ACGT'[(i // 4 + i) % 4] for i in range(12)]
    with tempfile.TemporaryDirectory(prefix='zuck-family-context-') as temporary:
        root = Path(temporary)
        cycles = root / 'cycles.txt'
        cycles.write_text(''.join(r + s + r[:22] + '\n'
            for i, s in enumerate(spacers) for r in [alternate if i == 1 else repeat]))
        reads = root / 'reads.fastq'

        def run(name, options=()):
            out = root / name
            proc = subprocess.run([str(binary), '--cycles', str(cycles), '--output-folder', str(out),
                '--threads', '1', '--ram', '1G', '--read-linked-groups', 'false', *options],
                capture_output=True, text=True, timeout=45)
            assert proc.returncode == 0, proc.stdout + proc.stderr
            return out

        absent = run('absent')
        entries = table(absent / 'entries.tsv')
        assert len(entries) == 12, entries
        assert {row['repeat_variant'] for row in entries} == {repeat, alternate}
        assert {row['spacer'] for row in entries} == set(spacers)
        spans = []
        for row in entries:
            following = alternate if row['repeat_variant'] == repeat else repeat
            context = row['repeat_variant'] + row['spacer'] + following[:16]
            spans.extend([context + 'AA', 'CCC' + context])

        def write(quality):
            reads.write_text(''.join(f'@{i}\n{s}\n+\n{quality * len(s)}\n' for i, s in enumerate(spans)))

        write('I')
        exact = run('exact', ['--grouping-reads', str(reads)])
        enabled = run('enabled', ['--grouping-reads', str(reads), '--family-repeat-context', 'true'])
        assert not table(exact / 'supported_entries.tsv')
        selected = table(enabled / 'supported_entries.tsv')
        assert len(selected) == 12
        assert {row['spacer'] for row in selected} == set(spacers)
        for name in ('arrays.tsv', 'entries.tsv'):
            assert (absent / name).read_bytes() == (exact / name).read_bytes() == (enabled / name).read_bytes()
        manifest = json.loads((enabled / 'run_manifest.json').read_text())
        assert manifest['settings']['family_repeat_context'] is True
        assert manifest['counts']['confidence_family_context_alternative_rows'] == 12
        assert 'same_family' in manifest['supported_call_interpretation']['read_span']
        assert json.loads((exact / 'run_manifest.json').read_text())['settings']['family_repeat_context'] is False

        config = root / 'settings.conf'
        config.write_text('family_repeat_context=true\n')
        configured = run('configured', ['--settings', str(config), '--grouping-reads', str(reads)])
        assert (configured / 'supported_entries.tsv').read_bytes() == (enabled / 'supported_entries.tsv').read_bytes()
        overridden = run('overridden', ['--family-repeat-context', 'false', '--settings', str(config),
                                      '--grouping-reads', str(reads)])
        assert (overridden / 'supported_entries.tsv').read_bytes() == (exact / 'supported_entries.tsv').read_bytes()
        write('!')
        low_quality = run('low_quality', ['--grouping-reads', str(reads), '--family-repeat-context', 'true'])
        assert not table(low_quality / 'supported_entries.tsv')
        print(json.dumps({'fixture': 'mixed_observed_repeat_alleles', 'generating_spacers': 12,
            'default_supported': 0, 'family_context_supported': len(selected),
            'low_quality_supported': 0, 'inclusive_rows_preserved': True}))


if __name__ == '__main__':
    main()
