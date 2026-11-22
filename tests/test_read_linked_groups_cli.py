"""Exercise grouping evidence, input parsing and transactional CLI semantics."""
import csv
import gzip
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def main():
    binary = Path(sys.argv[1]).resolve()
    rng = random.Random(872381)
    sequence = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    canonical = lambda s: min(s, rc(s))
    with tempfile.TemporaryDirectory(prefix='zuck-read-links-cli-') as tmp:
        root = Path(tmp)
        repeat = sequence(30)
        spacers = ['ACGT'[i % 4] + sequence(30) + 'ACGT'[(i // 4 + i) % 4] for i in range(12)]
        cycles = root / 'cycles.txt'
        cycles.write_text('\n'.join(repeat + s + (repeat + s)[:22] for s in spacers) + '\n')
        output = root / 'result'

        def run(extra=(), dest=output, raw=False, good=True):
            command = [str(binary), '--output-folder', str(dest), '--threads', '1', '--ram', '1G']
            if not raw:
                command += ['--cycles', str(cycles)]
            result = subprocess.run(command + list(extra), capture_output=True, text=True, timeout=60)
            if good:
                assert result.returncode == 0, result.stdout + result.stderr
            else:
                assert result.returncode != 0, 'Invalid grouping input succeeded'
            return result

        run()
        baseline = {name: (output / name).read_bytes() for name in ('arrays.tsv', 'entries.tsv')}
        entries = rows(output / 'entries.tsv')
        arrays = rows(output / 'arrays.tsv')
        chosen = None
        for array in arrays:
            members = [e for e in entries if e['array_id'] == array['array_id']]
            distinct = {}
            for entry in members:
                if len(entry['spacer']) >= 21:
                    distinct.setdefault(canonical(entry['spacer']), entry['spacer'])
            if len(distinct) >= 4:
                chosen = array, members, list(distinct.values())[:3]
                break
        assert chosen, 'Fixture needs a parent with more than three spacers'
        array, members, selected = chosen
        observed_repeat = array['repeat_consensus']
        reads = [selected[i][-n:] + observed_repeat + selected[i + 1][:n]
                 for i in range(2) for n in (20, 21)]
        text = ''.join(f'@read{i}\n{s}\n+\n' + 'I' * len(s) + '\n' for i, s in enumerate(reads))
        fastq = root / 'reads.fq'
        fastq.write_text(text)
        run(['--grouping-reads', str(fastq)])
        for name, data in baseline.items():
            assert (output / name).read_bytes() == data, 'Linking changed primary evidence'
        membership = rows(output / 'read_linked_groups.tsv')
        edge_rows = rows(output / 'read_linked_edges.tsv')
        assert len({r['group_id'] for r in membership}) == 1
        expected = {e['entry_index'] for e in members if canonical(e['spacer']) in {canonical(s) for s in selected}}
        assert {r['entry_index'] for r in membership} == expected
        assert all(r['parent_array_id'] == array['array_id'] for r in membership)
        assert len(edge_rows) == 2
        assert all(int(e['distinct_sequence_units']) == 2 and int(e['distinct_start_offsets']) == 2
                   and e['quality_mode'] == 'ee1_q10_exact_flanks' for e in edge_rows)
        manifest = json.loads((output / 'run_manifest.json').read_text())
        assert manifest['settings']['read_linked_groups_applied']
        assert manifest['counts']['read_linked_groups'] == 1
        assert manifest['read_linkage_interpretation']['independent_molecules'] is False
        assert {'read_linked_groups.tsv', 'read_linked_edges.tsv'} <= {f['path'] for f in manifest['files']}
        assert manifest['postprocessing_thresholds']['read_linked_minimum_start_offsets'] == 2
        before = {p.name: p.read_bytes() for p in output.iterdir() if p.is_file()}
        bad = root / 'bad.fq'
        for content, extra in [('@x\nACGT\n+\nII\n', []), ('@x\nACGT\n', []),
                               ('@different\nACGT\n+\nIIII\n', [str(fastq)])]:
            bad.write_text(content)
            run(['--grouping-reads', *extra, str(bad)], good=False)
            assert {p.name: p.read_bytes() for p in output.iterdir() if p.is_file()} == before
        run(['--grouping-reads', str(root / 'absent.fq')], good=False)
        assert {p.name: p.read_bytes() for p in output.iterdir() if p.is_file()} == before
        zipped = root / 'reads.fq.gz'
        with gzip.open(zipped, 'wt') as stream:
            stream.write(text)
        run(['--grouping-reads', str(zipped)])
        assert rows(output / 'read_linked_groups.tsv') == membership
        assert rows(output / 'read_linked_edges.tsv') == edge_rows
        # A missing gzip footer can leave every FASTQ record intact, so record
        # parsing alone cannot establish successful input decoding.
        footerless = root / 'no-footer.fq.gz'
        with gzip.open(footerless, 'wt') as stream:
            stream.write(text * 2000)
        footerless.write_bytes(footerless.read_bytes()[:-8])
        current = {p.name: p.read_bytes() for p in output.iterdir() if p.is_file()}
        run(['--grouping-reads', str(footerless)], good=False)
        assert {p.name: p.read_bytes() for p in output.iterdir() if p.is_file()} == current
        config = root / 'settings.txt' 
        config.write_text(f'read-linked-groups=false\ngrouping-reads={fastq}\n')
        run(['--settings', str(config), '--read-linked-groups', 'true'])
        assert rows(output / 'read_linked_groups.tsv') == membership
        run(['--settings', str(config)])
        assert not (output / 'read_linked_groups.tsv').exists()
        assert not (output / 'read_linked_edges.tsv').exists()
        assert not json.loads((output / 'run_manifest.json').read_text())['settings']['read_linked_groups_applied']
        run(['--grouping-reads', str(fastq)])
        run(['--grouping-reads', str(fastq), '--read-linked-groups', 'false'])
        assert not (output / 'read_linked_groups.tsv').exists(), 'Disabled run retained stale views'
        fasta = root / 'reads.fa'
        fasta.write_text(''.join(f'>read{i}\n{s}\n' for i, s in enumerate(reads)))
        run(['--grouping-reads', str(fasta)])
        manifest = json.loads((output / 'run_manifest.json').read_text())
        assert manifest['counts']['read_linked_records_without_quality'] == len(reads)
        assert not rows(output / 'read_linked_groups.tsv'), 'FASTA invented quality-supported links'
        # Read inputs may not alias files that the result transaction replaces.
        # Check both spelling and filesystem identity while preserving the bundle.
        for indirect in (False, True):
            collision = root / ('collision-link' if indirect else 'collision-direct')
            collision.mkdir()
            source = collision / 'read_linked_groups.tsv'
            source.write_text(text)
            (collision / 'keep.txt').write_text('prior user content')
            protected = root / 'input-hardlink.fq' if indirect else source
            if indirect:
                protected.hardlink_to(source)
            prior = {p.name: p.read_bytes() for p in collision.iterdir()}
            run(['--grouping-reads', str(protected)], dest=collision, good=False)
            assert {p.name: p.read_bytes() for p in collision.iterdir()} == prior
            assert protected.read_text() == text
        raw_output = root / 'raw' 
        run(['--input-files', str(fastq)], dest=raw_output, raw=True)
        raw_manifest = json.loads((raw_output / 'run_manifest.json').read_text())
        assert raw_manifest['grouping_read_inputs'] == [str(fastq)]
        assert raw_manifest['settings']['read_linked_groups_applied']
        assert raw_manifest['counts']['read_linked_records'] == len(reads)
        print('Read-linked grouping preserves rows, validates reads and replaces sidecars transactionally')


if __name__ == '__main__':
    main()
