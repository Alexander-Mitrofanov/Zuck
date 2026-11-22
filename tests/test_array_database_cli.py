"""Independent synthetic observations, opt-in isolation, and output integrity."""
from collections import Counter
import csv
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def main():
    binary = Path(sys.argv[1]).resolve()
    rng = random.Random(61293)
    seq = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    repeat = seq(30)
    spacers = ['A' + seq(30) + 'T' for _ in range(4)]
    header = 'array_id\taccession\tcategory\tconfidence\trepeat_count\tspacer_count\trepeats\tspacers\n'
    with tempfile.TemporaryDirectory(prefix='zuck-database-cli-') as tmp:
        root = Path(tmp)
        cycles = root / 'cycles.txt'
        cores = [repeat + s for s in spacers]
        cycles.write_text(''.join(c + c[:22] + '\n' for c in cores))
        database = root / 'database.tsv'
        database.write_text(header + f'known\tGCA_1.1\tBona-fide\t0.9\t5\t4\t{repeat}\t' + ','.join(spacers) + '\n')
        reads = root / 'reads.fastq'
        generated_reads = [prefix + repeat + s + repeat[:16] + suffix
                           for s in spacers for prefix, suffix in [('', 'AA'), ('CCC', '')]]
        reads.write_text(''.join(f'@{i}\n{s}\n+\n' + 'I'*len(s) + '\n'
                                 for i,s in enumerate(generated_reads)))
        def run(folder, extra=(), success=True, with_reads=True):
            read_args = ['--grouping-reads', str(reads)] if '--array-database' in extra and with_reads else []
            result = subprocess.run([str(binary), '--cycles', str(cycles), '--output-folder',
                str(folder), '--threads', '1', '--ram', '1G', *read_args, *extra], capture_output=True, text=True, timeout=30)
            assert (result.returncode == 0) == success, result.stdout + result.stderr
            return result
        def rows(folder):
            with (folder / 'entries.tsv').open() as handle:
                return Counter((r['array_id'], r['entry_index'], r['repeat_variant'], r['spacer'])
                               for r in csv.DictReader(handle, delimiter='\t'))
        off, on = root / 'off', root / 'on'
        run(off)
        run(on, ['--array-database', str(database)])
        assert not rows(off) - rows(on), 'DB guidance must preserve all original rows and IDs'
        evidence = list(csv.DictReader((on / 'database_evidence.tsv').open(), delimiter='\t'))
        assert evidence, 'Fixture must recover exact boundaries over conserved spacer ends'
        assert {r['spacer'] for r in evidence} == set(spacers)
        assert all(r['repeat_variant'] == repeat for r in evidence)
        assert all(r['read_validation_applied'] == '1' and int(r['read_distinct_sequence_units']) >= 2 and int(r['read_distinct_start_offsets']) >= 2 for r in evidence)
        assert all(r['evidence_kind'] == 'known_same_array' and r['read_span_rule'] == 'whole_R_S_nextR16_q10_ee1' for r in evidence)
        run(root / 'no-reads', ['--array-database', str(database)], False, with_reads=False)
        run(root / 'diagnostic', ['--array-database', str(database), '--database-read-validation', 'false'], with_reads=False)
        diagnostic = ['--array-database', str(database), '--database-read-validation', 'false']
        run(root / 'diagnostic-reads', diagnostic)
        run(root / 'diagnostic-candidates', diagnostic + ['--candidate-evidence', 'true'])
        assert rows(root / 'diagnostic-reads') == rows(root / 'diagnostic-candidates')
        with (root / 'diagnostic-candidates/boundary_evidence.tsv').open() as f:
            assert list(csv.DictReader(f, delimiter='\t')) == [], 'Diagnostic collection cannot enable promotion'
        run(root / 'invalid-novel', ['--array-database', str(database), '--database-read-validation', 'false', '--database-novel-spacers', 'true'], False, with_reads=False)
        manifest = json.loads((on / 'run_manifest.json').read_text())
        assert manifest['array_database_input']['path'] == str(database)
        assert manifest['settings']['database_read_validation'] is True
        assert manifest['settings']['database_include_possible'] is False
        assert manifest['settings']['database_novel_spacers'] is False
        assert manifest['counts']['database_groups_added'] == 1
        assert manifest['counts']['database_new_canonical_spacers'] == 4
        before = {p.name: p.read_bytes() for p in on.iterdir() if p.is_file()}
        database.write_text(header + 'malformed\n')
        run(on, ['--array-database', str(database)], False)
        assert before == {p.name: p.read_bytes() for p in on.iterdir() if p.is_file()}, 'Failure must retain previous outputs'
        database.write_text(header)
        run(root / 'empty-db', ['--array-database', str(database)])
        assert rows(root / 'empty-db') == rows(off)
        run(root / 'legacy', ['--array-database', str(database), '--boundary-mode', 'legacy'], False)
        run(root / 'missing', ['--array-database', str(root / 'absent')], False)
        alias = root / 'alias'
        alias.mkdir()
        alias_db = alias / 'entries.tsv'
        alias_db.write_text(header)
        run(alias, ['--array-database', str(alias_db)], False)
        assert alias_db.read_text() == header, 'Database must not alias managed output'
        hard = root / 'hard.tsv'
        hard.hardlink_to(alias_db)
        run(alias, ['--array-database', str(hard)], False)
        # Fail before graph/cycle generation can overwrite an input database.
        generated = root / 'generated'
        generated.mkdir()
        generated_db = generated / 'cycles.txt'
        generated_db.write_text(header)
        for source in (generated_db, root / 'generated-hard.tsv'):
            if source != generated_db:
                source.hardlink_to(generated_db)
            result = subprocess.run([str(binary), '--input-files', str(cycles),
                '--output-folder', str(generated), '--array-database', str(source),
                '--threads', '1', '--ram', '1G'], capture_output=True, text=True, timeout=30)
            assert result.returncode != 0 and 'conflicts with generated cycles.txt' in result.stderr
            assert generated_db.read_text() == header
        run(on)
        assert not (on / 'database_evidence.tsv').exists(), 'DB-off rerun removes stale evidence'
        assert not (on / 'database_groups.tsv').exists() and not (on / 'database_group_members.tsv').exists()
        print('Database CLI preserves de novo evidence, recovers observed boundaries, and protects results')


if __name__ == '__main__':
    main()
