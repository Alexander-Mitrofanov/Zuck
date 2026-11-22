"""Optional output limits apply to actual rows, consensus and aligned sidecars."""
import csv
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def main():
    binary = Path(sys.argv[1]).resolve()
    rng = random.Random(231113)
    dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    repeat = dna(40)
    spacers = ['ACGT'[i % 4] + dna(29 + i % 3) + 'TGCA'[i % 4] for i in range(12)]
    with tempfile.TemporaryDirectory(prefix='zuck-lengths-') as temporary:
        root = Path(temporary)
        cycles = root / 'cycles.txt'
        cycles.write_text(''.join((repeat + s + repeat[:22] + '\n') * 3 for s in spacers))
        database = root / 'database.tsv'
        database.write_text('array_id\taccession\tcategory\tconfidence\trepeat_count\tspacer_count\trepeats\tspacers\n'
                           + f'known\tGCA_1.1\tBona-fide\t0.9\t13\t12\t{repeat}\t' + ','.join(spacers) + '\n')
        config = root / 'settings.txt'
        config.write_text('spacer-min-length=40\nrepeat-min-length=45\n')
        reads = root / 'reads.fastq'
        spans = [prefix + repeat + s + repeat[:16] + suffix for s in spacers
                 for prefix, suffix in [('', 'AA'), ('CCC', '')]]
        reads.write_text(''.join(f'@{i}\n{s}\n+\n' + 'I'*len(s) + '\n' for i, s in enumerate(spans)))

        def run(name, extra=(), success=True):
            output = root / name
            result = subprocess.run([str(binary), '--cycles', str(cycles), '--output-folder',
                str(output), '--threads', '1', '--ram', '1G', '--grouping-reads', str(reads), *extra],
                capture_output=True, text=True, timeout=30)
            assert (result.returncode == 0) == success, result.stdout + result.stderr
            return output

        baseline = run('baseline')
        entries = table(baseline / 'entries.tsv')
        assert len(entries) == 12 and {len(r['spacer']) for r in entries} == {31, 32, 33}
        for name, flag, value in [('smin', 'spacer-min-length', '34'), ('smax', 'spacer-max-length', '30'),
                                  ('rmin', 'repeat-min-length', '41'), ('rmax', 'repeat-max-length', '39')]:
            output = run(name, ['--' + flag, value, '--candidate-evidence', 'true'])
            assert table(output / 'entries.tsv') == []
            assert table(output / 'arrays.tsv') == []
            assert table(output / 'element_confidence.tsv') == []
            assert all(r['passes_length_constraints'] == '0' for r in table(output / 'candidate_alleles.tsv'))
        assert table(run('settings', ['--settings', str(config)]) / 'entries.tsv') == []
        override = run('override', ['--settings', str(config), '--spacer-min-length', '0', '--repeat-min-length', '0'])
        assert (override / 'entries.tsv').read_bytes() == (baseline / 'entries.tsv').read_bytes()
        # Conserved spacer ends cause a different de novo cut. The catalogue
        # supplies exact observed pairs absent from the de novo output, so its
        # evidence vector really is exercised by filtering mixed spacer sizes.
        db_repeat = dna(40)
        db_spacers = ['A' + dna(29 + i % 3) + 'T' for i in range(12)]
        with cycles.open('a') as stream:
            stream.write(''.join((db_repeat + s + db_repeat[:22] + '\n') * 3 for s in db_spacers))
        database.write_text('array_id\taccession\tcategory\tconfidence\trepeat_count\tspacer_count\trepeats\tspacers\n'
                           + f'known\tGCA_1.1\tBona-fide\t0.9\t13\t12\t{db_repeat}\t' + ','.join(db_spacers) + '\n')
        with reads.open('a') as stream:
            for i, s in enumerate(db_spacers):
                for j, (prefix, suffix) in enumerate([('', 'AA'), ('CCC', '')]):
                    seq = prefix + db_repeat + s + db_repeat[:16] + suffix
                    stream.write(f'@db-{i}-{j}\n{seq}\n+\n' + 'I'*len(seq) + '\n')
            linked = [s for s in spacers if len(s) == 32][:3]
            for i in range(2):
                for flank in (20, 21):
                    seq = linked[i][-flank:] + repeat + linked[i + 1][:flank]
                    stream.write(f'@link-{i}-{flank}\n{seq}\n+\n' + 'I'*len(seq) + '\n')
        bounded = run('bounded', ['--spacer-min-length', '32', '--spacer-max-length', '32',
                                 '--repeat-min-length', '40', '--repeat-max-length', '40',
                                 '--array-database', str(database)])
        rows = table(bounded / 'entries.tsv')
        assert rows and all(len(r['spacer']) == 32 and len(r['repeat_variant']) == 40 for r in rows)
        allowed = {(r['array_id'], r['entry_index']): (r['repeat_variant'], r['spacer']) for r in rows}
        for name in ['element_confidence.tsv', 'database_evidence.tsv']:
            evidence = table(bounded / name)
            assert evidence, name
            for row in evidence:
                assert allowed[(row['array_id'], row['entry_index'])] == (row['repeat_variant'], row['spacer'])
        membership = table(bounded / 'read_linked_groups.tsv')
        assert len(membership) == 3
        assert {allowed[(r['parent_array_id'], r['entry_index'])][1] for r in membership} == set(linked)
        manifest = json.loads((bounded / 'run_manifest.json').read_text())
        assert manifest['output_length_constraints']['applied'] is True
        assert manifest['output_length_constraints']['spacer_minimum'] == 32
        for index, args in enumerate([
            ['--spacer-min-length', '-1'], ['--repeat-max-length', '40junk'],
            ['--spacer-max-length', '2147483648'], ['--repeat-min-length', '42', '--repeat-max-length', '41'],
            ['--spacer-min-length', '33', '--spacer-max-length', '32'], ['--traversal-partitions', '0'],
            ['--low-abundance', 'maybe']]):
            run('invalid-' + str(index), args, False)
    print('All length bounds, settings precedence, CLI validation and sidecar alignment passed')


if __name__ == '__main__':
    main()
