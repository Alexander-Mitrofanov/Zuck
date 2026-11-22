"""Confidence is a complete, transaction-managed annotation over unchanged calls."""
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
    rng = random.Random(173924)
    dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    repeat = dna(30)
    cores = [repeat + 'ACGT'[i % 4] + dna(31) + 'ACGT'[(i // 4 + i) % 4]
             for i in range(12)]
    with tempfile.TemporaryDirectory(prefix='zuck-confidence-cli-') as temporary:
        root = Path(temporary)
        cycles = root / 'cycles.txt'
        cycles.write_text(''.join(c + c[:22] + '\n' for c in cores))
        reads = root / 'reads.fastq'
        spans = [prefix + c + repeat[:16] + suffix for c in cores
                 for prefix, suffix in [('', 'AA'), ('CCC', '')]]
        reads.write_text(''.join(f'@{i}\n{s}\n+\n' + 'I' * len(s) + '\n'
                                 for i, s in enumerate(spans)))
        def run(output, extra=(), success=True):
            result = subprocess.run([str(binary), '--cycles', str(cycles),
                '--output-folder', str(output), '--threads', '1', '--ram', '1G',
                '--read-linked-groups', 'false', *extra], capture_output=True, text=True, timeout=30)
            assert (result.returncode == 0) == success, result.stdout + result.stderr
            return result
        off, on = root / 'off', root / 'on'
        run(off, ['--confidence-scores', 'false'])
        run(on, ['--grouping-reads', str(reads)])
        for name in ['arrays.tsv', 'entries.tsv']:
            assert (off / name).read_bytes() == (on / name).read_bytes()
        entries = table(on / 'entries.tsv')
        scores = table(on / 'element_confidence.tsv')
        assert entries and len(entries) == len(scores)
        fields = ['array_id', 'entry_index', 'repeat_variant', 'spacer']
        assert [[r[k] for k in fields] for r in entries] == [[r[k] for k in fields] for r in scores]
        score_bytes = (on / 'element_confidence.tsv').read_bytes()
        run(on, ['--grouping-reads', str(reads)])
        assert score_bytes == (on / 'element_confidence.tsv').read_bytes()
        manifest = json.loads((on / 'run_manifest.json').read_text())
        assert manifest['scoring']['primary_rows_preserved'] is True
        assert manifest['scoring']['calibrated_probability'] is False
        assert 'element_confidence.tsv' in {f['path'] for f in manifest['files']}
        before = {p.name: p.read_bytes() for p in on.iterdir() if p.is_file()}
        reads.write_text('@truncated\nACGT\n+\nI\n')
        run(on, ['--grouping-reads', str(reads)], False)
        assert before == {p.name: p.read_bytes() for p in on.iterdir() if p.is_file()}
        # Input aliases are protected even when read-linked grouping is off.
        alias = root / 'input-alias'
        alias.mkdir()
        alias_reads = alias / 'element_confidence.tsv'
        alias_reads.write_text('>read\nACGT\n')
        run(alias, ['--grouping-reads', str(alias_reads)], False)
        assert alias_reads.read_text() == '>read\nACGT\n'
        run(on, ['--confidence-scores', 'false'])
        assert not (on / 'element_confidence.tsv').exists()
        assert json.loads((on / 'run_manifest.json').read_text())['scoring'] is None
        run(root / 'invalid', ['--confidence-scores', 'maybe'], False)
        empty = root / 'empty-output'
        cycles.write_text('')
        run(empty)
        assert table(empty / 'element_confidence.tsv') == []
        print('Confidence row coverage, isolation, repeatability, rollback and stale cleanup passed')


if __name__ == '__main__':
    main()
