"""Retain a planted internal SNP before collapse without changing primary calls."""
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
    spacers = ['ACGT'[i % 4] + dna(31) + 'TGCA'[i % 4] for i in range(8)]
    original = spacers[0]
    spacers.append(original[:25] + 'ACGT'[('ACGT'.index(original[25]) + 1) % 4] + original[26:])
    with tempfile.TemporaryDirectory(prefix='zuck-candidates-') as temporary:
        root = Path(temporary)
        cycles = root / 'cycles.txt'
        records = [repeat + s + repeat[:22] + '\n' for s in spacers]
        cycles.write_text(''.join(records * 3))
        reads = root / 'reads.fastq'
        spans = [prefix + repeat + s + repeat[:16] + suffix for s in spacers
                 for prefix, suffix in [('', 'AA'), ('CCC', '')]]
        reads.write_text(''.join(f'@{i}\n{s}\n+\n' + 'I'*len(s) + '\n' for i, s in enumerate(spans)))

        def run(output, extra=(), success=True):
            result = subprocess.run([str(binary), '--cycles', str(cycles),
                '--output-folder', str(output), '--threads', '1', '--ram', '1G',
                '--read-linked-groups', 'false', '--confidence-scores', 'false', *extra],
                capture_output=True, text=True, timeout=30)
            assert (result.returncode == 0) == success, result.stdout + result.stderr
            return result

        off, on = root / 'off', root / 'on'
        run(off)
        args = ['--candidate-evidence', 'true', '--grouping-reads', str(reads)]
        run(on, args)
        for name in ['arrays.tsv', 'entries.tsv']:
            assert (off / name).read_bytes() == (on / name).read_bytes()
        entries = table(on / 'entries.tsv')
        candidates = table(on / 'candidate_alleles.tsv')
        assert len(entries) == 8 and len(candidates) == 9
        boundary = table(on / 'boundary_candidates.tsv')
        assert len(boundary) == 9 * 49
        assert {(int(r['repeat_start_shift']), int(r['repeat_end_shift'])) for r in boundary} == {
            (a,b) for a in range(-3,4) for b in range(-3,4)}
        assert all(float(r['same_array_score']) == 0 for r in boundary)
        run(root / 'radius-one', args + ['--boundary-radius', '1'])
        assert len(table(root / 'radius-one/boundary_candidates.tsv')) == 9 * 9
        run(root / 'bad-radius', ['--boundary-radius', '4'], False)
        run(root / 'bad-policy', ['--boundary-policy', 'truth'], False)
        assert {(r['repeat_variant'], r['spacer']) for r in candidates} == {(repeat, s) for s in spacers}
        omitted = [r for r in candidates if r['emitted_primary_pair'] == '0']
        assert len(omitted) == 1 and omitted[0]['spacer'] in (original, spacers[-1])
        assert all(r['context_score'] == '0.5' and r['context_sequence_units'] == '2'
                   and r['context_start_offsets'] == '2' for r in candidates)
        assert all(r['traversal_observations'] == '3' and r['distinct_circular_cores'] == '1'
                   for r in candidates)
        manifest = json.loads((on / 'run_manifest.json').read_text())
        assert manifest['counts']['candidate_context_supported_unemitted_pairs'] == 1
        assert manifest['candidate_evidence_interpretation']['selection_uses_read_scores'] is False
        assert 'candidate_alleles.tsv' in {r['path'] for r in manifest['files']}
        cycles.write_text(''.join(records * 6))
        run(root / 'duplicated', args)
        duplicated = table(root / 'duplicated/candidate_alleles.tsv')
        for first, second in zip(candidates, duplicated):
            assert second['traversal_observations'] == '6'
            assert {k: v for k, v in first.items() if k != 'traversal_observations'} == {
                k: v for k, v in second.items() if k != 'traversal_observations'}
        run(root / 'no-reads', ['--candidate-evidence', 'true'])
        assert all(r['context_score'] == 'NA' for r in table(root / 'no-reads/candidate_alleles.tsv'))
        fasta = root / 'reads.fa'
        fasta.write_text(''.join(f'>r{i}\n{s}\n' for i, s in enumerate(spans)))
        run(root / 'no-quality', ['--candidate-evidence', 'true', '--grouping-reads', str(fasta)])
        assert all(r['context_score'] == 'NA' for r in table(root / 'no-quality/candidate_alleles.tsv'))
        before = {p.name: p.read_bytes() for p in on.iterdir() if p.is_file()}
        reads.write_text('@truncated\nACGT\n+\nI\n')
        run(on, args, False)
        assert before == {p.name: p.read_bytes() for p in on.iterdir() if p.is_file()}
        alias = root / 'alias'
        alias.mkdir()
        alias_reads = alias / 'candidate_alleles.tsv'
        alias_reads.write_text('>read\nACGT\n')
        run(alias, ['--candidate-evidence', 'true', '--grouping-reads', str(alias_reads)], False)
        assert alias_reads.read_text() == '>read\nACGT\n'
        run(on)
        assert not (on / 'candidate_alleles.tsv').exists()
        assert not (on / 'boundary_candidates.tsv').exists()
        assert not (on / 'boundary_evidence.tsv').exists()
        run(root / 'legacy', ['--candidate-evidence', 'true', '--boundary-mode', 'legacy'], False)
        run(root / 'invalid', ['--candidate-evidence', 'maybe'], False)
        cycles.write_text('')
        run(root / 'empty', ['--candidate-evidence', 'true'])
        assert table(root / 'empty/candidate_alleles.tsv') == []
        # Shared spacer termini are absorbed into the inferred repeat. The
        # real catalogue cut needs both ends of the spacer lengthened by 3 nt.
        shared_spacers = ['AAA' + base + dna(24) + end + 'CCC'
                          for base, end in zip('ACGT', 'TGCA')]
        cycles.write_text(''.join(repeat+s+repeat[:22]+'\n' for s in shared_spacers))
        spans = [prefix+repeat+s+repeat[:16]+suffix for s in shared_spacers
                 for prefix,suffix in [('', 'AA'),('CCC','')]]
        reads.write_text(''.join(f'@{i}\n{s}\n+\n'+ 'I'*len(s)+'\n' for i,s in enumerate(spans)))
        db = root/'catalogue.tsv'
        db.write_text('array_id\taccession\tcategory\tconfidence\trepeat_count\tspacer_count\trepeats\tspacers\n'
                      +f'known\tGCA_synthetic\tBona-fide\t0.9\t5\t4\t{repeat}\t'+','.join(shared_spacers)+'\n')
        run(root/'longer-spacers',args+['--array-database',str(db)])
        boundary = table(root/'longer-spacers/boundary_candidates.tsv')
        recovered = [r for r in boundary if r['repeat_variant']==repeat and r['spacer'] in shared_spacers]
        assert {r['spacer'] for r in recovered} == set(shared_spacers)
        assert all((int(r['repeat_start_shift']),int(r['repeat_end_shift']))==(3,-3) for r in recovered)
        assert all(float(r['same_array_score'])>=.5 for r in recovered)
        # Varying repeat termini produce the opposite error: each initial
        # spacer contains three repeat bases at each end and must be shortened.
        repeat_variants = [a*3+repeat[3:-3]+b*3 for a,b in zip('ACGT','TGCA')]
        true_spacers = [a+dna(30)+b for a,b in zip('ACGT','TGCA')]
        true_pairs = set(zip(repeat_variants,true_spacers))
        cycles.write_text(''.join(r+s+r[:22]+'\n' for r,s in true_pairs))
        spans = [prefix+r+s+r[:16]+suffix for r,s in true_pairs
                 for prefix,suffix in [('', 'AA'),('CCC','')]]
        reads.write_text(''.join(f'@{i}\n{s}\n+\n'+'I'*len(s)+'\n' for i,s in enumerate(spans)))
        db.write_text('array_id\taccession\tcategory\tconfidence\trepeat_count\tspacer_count\trepeats\tspacers\n'
                      +'known\tGCA_synthetic\tBona-fide\t0.9\t5\t4\t'+','.join(repeat_variants)+'\t'+','.join(true_spacers)+'\n')
        run(root/'shorter-spacers',args+['--array-database',str(db)])
        boundary = table(root/'shorter-spacers/boundary_candidates.tsv')
        recovered = [r for r in boundary if (r['repeat_variant'],r['spacer']) in true_pairs]
        assert {(r['repeat_variant'],r['spacer']) for r in recovered} == true_pairs
        assert all((int(r['repeat_start_shift']),int(r['repeat_end_shift']))==(-3,3) for r in recovered)
        assert all(float(r['same_array_score'])>=.5 for r in recovered)
    print('Pre-collapse SNP evidence, score independence, isolation and bundle integrity passed')


if __name__ == '__main__':
    main()
