"""Exercise native call profiles, evidence missingness, and result isolation."""
import csv
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


FILES = ('call_decisions.tsv', 'supported_entries.tsv',
         'supported_arrays.tsv', 'supported_spacers.fasta')
IDENTITY = ('array_id', 'entry_index', 'repeat_variant', 'spacer')


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def canonical(sequence):
    return min(sequence, sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1])


def main():
    binary = Path(sys.argv[1]).resolve()
    rng = random.Random(987431)
    dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    repeat = dna(30)
    spacers = ['ACGT'[i % 4] + dna(n - 2) + 'ACGT'[(i // 4 + i) % 4]
               for i, n in enumerate([25, 26, 32, 33, 49, 50, 51, 32, 33, 34, 35, 36])]
    with tempfile.TemporaryDirectory(prefix='zuck-supported-cli-') as temporary:
        root = Path(temporary)
        cycles = root / 'cycles.txt'
        cycles.write_text(''.join(repeat + spacer + repeat[:22] + '\n' for spacer in spacers))

        def run(output, extra=(), success=True):
            result = subprocess.run([str(binary), '--cycles', str(cycles),
                '--output-folder', str(output), '--threads', '1', '--ram', '1G',
                '--read-linked-groups', 'false', *extra], capture_output=True, text=True, timeout=30)
            assert (result.returncode == 0) == success, result.stdout + result.stderr
            return result

        def ids(rows):
            return [tuple(row[key] for key in IDENTITY) for row in rows]

        def audit(output):
            entries = table(output / 'entries.tsv')
            decisions = table(output / 'call_decisions.tsv')
            selected = table(output / 'supported_entries.tsv')
            assert ids(entries) == ids(decisions)
            assert ids(selected) == ids([row for row in decisions if row['selected'] == '1'])
            wanted = sorted({canonical(row['spacer']) for row in selected})
            lines = (output / 'supported_spacers.fasta').read_text().splitlines()
            assert lines[1::2] == wanted
            assert lines[::2] == [f'>Zuck_supported_{i + 1}' for i in range(len(wanted))]
            by_array = {}
            for row in selected:
                by_array.setdefault(row['array_id'], []).append(row)
            groups = table(output / 'supported_arrays.tsv')
            assert {row['array_id'] for row in groups} == set(by_array)
            for group in groups:
                rows = by_array[group['array_id']]
                distinct = len({canonical(row['spacer']) for row in rows})
                assert int(group['entry_count']) == len(rows)
                assert int(group['distinct_canonical_spacer_count']) == distinct
                assert group['group_status'] == ('singleton' if distinct == 1 else 'multiple_spacers')
            return decisions, selected

        absent = root / 'absent'
        run(absent)
        decisions, selected = audit(absent)
        assert decisions and not selected
        assert all(row['context_status'] == 'unknown' and row['context_observability'] == 'no_reads'
                   and 'no_reads' in row['reason'] for row in decisions)
        assert all(row['profile'] == 'standard' for row in decisions)
        manifest = json.loads((absent / 'run_manifest.json').read_text())
        assert manifest['settings']['call_profile'] == 'standard'
        assert manifest['supported_call_interpretation']['primary_rows_preserved'] is True
        assert manifest['supported_call_interpretation']['singletons_retained_and_labeled'] is True
        assert set(FILES) <= {row['path'] for row in manifest['files']}

        dedup_empty = root / 'dedup-empty'
        run(dedup_empty, ['--circular-unit-dedup', 'true'])
        dedup_manifest = json.loads((dedup_empty / 'run_manifest.json').read_text())
        assert (dedup_empty / 'circular_unit_dedup.tsv').exists()
        assert dedup_manifest['settings']['circular_unit_dedup'] is True
        assert dedup_manifest['supported_call_interpretation']['circular_unit_dedup']['enabled'] is True
        assert 'circular_unit_dedup.tsv' in {row['path'] for row in dedup_manifest['files']}

        # Generate exact supporting contexts for the emitted rows, at distinct
        # offsets and with distinct full sequences. This tests CLI plumbing and
        # serialization independently of repeat-boundary choices on this fixture.
        reads = root / 'reads.fastq'
        contexts = [row['repeat_variant'] + row['spacer'] + row['repeat_variant'][:16]
                    for row in decisions]
        spans = [prefix + context + suffix for context in contexts
                 for prefix, suffix in [('', 'TT'), ('CCC', 'A')]]
        def fastq(quality='I'):
            reads.write_text(''.join(f'@{i}\n{s}\n+\n{quality * len(s)}\n'
                                     for i, s in enumerate(spans)))
        fastq()
        standard, context = root / 'standard', root / 'context'
        run(standard, ['--grouping-reads', str(reads)])
        standard_decisions, standard_rows = audit(standard)
        assert standard_rows
        assert all(26 <= len(row['spacer']) <= 50 for row in standard_rows)
        assert all(row['context_status'] == 'supported' for row in standard_decisions)
        assert all(row['selected'] == str(int(26 <= len(row['spacer']) <= 50))
                   for row in standard_decisions)
        assert any(row['reason'] == 'spacer_length_outside_profile' for row in standard_decisions)
        run(context, ['--call-profile', 'context', '--grouping-reads', str(reads)])
        context_decisions, context_rows = audit(context)
        assert len(context_rows) == len(context_decisions)
        assert all(row['reason'] == 'passes_profile' for row in context_decisions)
        for name in ('arrays.tsv', 'entries.tsv'):
            assert (absent / name).read_bytes() == (standard / name).read_bytes() == (context / name).read_bytes()

        # The call profile computes evidence even with the annotation disabled.
        unannotated = root / 'unannotated'
        run(unannotated, ['--confidence-scores', 'false', '--grouping-reads', str(reads)])
        assert not (unannotated / 'element_confidence.tsv').exists()
        for name in FILES:
            assert (unannotated / name).read_bytes() == (standard / name).read_bytes()
        assert json.loads((unannotated / 'run_manifest.json').read_text())['scoring'] is None

        # Settings are loaded first even when --settings follows a CLI override.
        config = root / 'settings.conf'
        config.write_text('call_profile=off\nconfidence-scores=false\n')
        configured = root / 'configured'
        run(configured, ['--settings', str(config)])
        assert not any((configured / name).exists() for name in FILES)
        run(configured, ['--call-profile', 'context', '--settings', str(config), '--grouping-reads', str(reads)])
        assert len(audit(configured)[1]) == len(context_rows)
        assert json.loads((configured / 'run_manifest.json').read_text())['settings']['call_profile'] == 'context'

        # Successful disable removes stale files, but retains the inclusive rows.
        original_entries = (configured / 'entries.tsv').read_bytes()
        run(configured, ['--call-profile', 'off', '--confidence-scores', 'false'])
        assert not any((configured / name).exists() for name in FILES)
        assert (configured / 'entries.tsv').read_bytes() == original_entries
        assert json.loads((configured / 'run_manifest.json').read_text())['supported_call_interpretation']['enabled'] is False

        for value in ('unknown', 'true', 'Standard'):
            run(root / ('invalid-' + value), ['--call-profile', value], False)
        run(root / 'missing-profile', ['--call-profile'], False)
        run(root / 'dedup-disabled-profile',
            ['--call-profile', 'off', '--circular-unit-dedup', 'true'], False)
        config.write_text('call-profile=invalid\n')
        run(root / 'invalid-settings', ['--settings', str(config)], False)
        run(root / 'missing-reads', ['--confidence-scores', 'false', '--grouping-reads', str(root / 'absent.fastq')], False)

        # Unknown evidence differs from observable evidence that fails quality.
        fasta = root / 'reads.fasta'
        fasta.write_text(''.join(f'>{i}\n{s}\n' for i, s in enumerate(spans)))
        no_quality = root / 'no-quality'
        run(no_quality, ['--grouping-reads', str(fasta)])
        unknown, selected = audit(no_quality)
        assert not selected and all(row['context_status'] == 'unknown' and
                                    row['context_observability'] == 'no_quality' for row in unknown)
        reads.write_text('@short\nACGT\n+\nIIII\n')
        short = root / 'short'
        run(short, ['--grouping-reads', str(reads)])
        unknown, selected = audit(short)
        assert not selected and all(row['context_observability'] == 'span_exceeds_max_quality_read'
                                    and row['context_status'] == 'unknown' for row in unknown)
        fastq('!')
        low_quality = root / 'low-quality'
        run(low_quality, ['--grouping-reads', str(reads)])
        failed, selected = audit(low_quality)
        assert not selected and all(row['context_status'] == 'insufficient' and
                                    row['context_observability'] == 'observable' for row in failed)

        # A malformed input must not partially replace any previous result.
        before = {p.name: p.read_bytes() for p in standard.iterdir() if p.is_file()}
        reads.write_text('@truncated\nACGT\n+\nI\n')
        run(standard, ['--grouping-reads', str(reads)], False)
        assert before == {p.name: p.read_bytes() for p in standard.iterdir() if p.is_file()}
        for name in FILES:
            output = root / ('alias-' + name)
            output.mkdir()
            source = output / name
            source.write_text('>read\nACGT\n')
            run(output, ['--grouping-reads', str(source)], False)
            assert source.read_text() == '>read\nACGT\n'

        cycles.write_text('')
        empty = root / 'empty'
        run(empty)
        assert audit(empty) == ([], [])
        print('Supported call profiles, evidence states, membership, isolation and transactions passed')


if __name__ == '__main__':
    main()
