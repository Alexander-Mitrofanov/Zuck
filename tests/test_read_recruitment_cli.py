"""Exercise discovery, clipped junctions, and bundle transactions through the CLI."""
import csv
import json
from pathlib import Path
import random
import struct
import subprocess
import sys
import tempfile


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def canonical(sequence):
    return min(sequence, sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1])


def xxh64(data):
    """Independent seed-zero byte checksum for the native bundle inventory."""
    mask = (1 << 64) - 1
    p1, p2, p3, p4, p5 = (11400714785074694791, 14029467366897019727,
                          1609587929392839161, 9650029242287828579, 2870177450012600261)
    def rot(x, n):
        x &= mask
        return ((x << n) | (x >> (64-n))) & mask
    def step(x, word):
        return (rot(x + word*p2, 31)*p1) & mask
    position = 0
    if len(data) >= 32:
        accumulators = [(p1+p2) & mask, p2, 0, (-p1) & mask]
        while position + 32 <= len(data):
            accumulators = [step(a, b) for a, b in zip(accumulators, struct.unpack_from('<4Q', data, position))]
            position += 32
        value = sum(rot(x, n) for x, n in zip(accumulators, (1, 7, 12, 18))) & mask
        for x in accumulators:
            value = ((value ^ step(0, x))*p1 + p4) & mask
    else:
        value = p5
    value = (value + len(data)) & mask
    while position + 8 <= len(data):
        value = (rot(value ^ step(0, struct.unpack_from('<Q', data, position)[0]), 27)*p1+p4) & mask
        position += 8
    if position + 4 <= len(data):
        value = (rot(value ^ ((struct.unpack_from('<I', data, position)[0]*p1) & mask), 23)*p2+p3) & mask
        position += 4
    for byte in data[position:]:
        value = (rot(value ^ ((byte*p5) & mask), 11)*p1) & mask
    value ^= value >> 33
    value = (value*p2) & mask
    value ^= value >> 29
    value = (value*p3) & mask
    value ^= value >> 32
    return f'{value:016x}'


def main():
    assert xxh64(b'') == 'ef46db3751d8e999'
    binary = Path(sys.argv[1]).resolve()
    rng = random.Random(873109)
    dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    repeat = dna(36)
    spacers = ['ACGT'[i % 4] + dna(30) + 'ACGT'[(i // 4 + i) % 4] for i in range(12)]
    truth = {canonical(s) for s in spacers}
    with tempfile.TemporaryDirectory(prefix='zuck-read-recruitment-') as temporary:
        root = Path(temporary)
        cycles, reads = root/'cycles.txt', root/'reads.fastq'
        cycles.write_text(''.join(repeat+s+repeat[:22]+'\n' for s in spacers[:10]))

        def write_reads(sequences):
            reads.write_text(''.join(f'@{i}\n{s}\n+\n'+ 'I'*len(s)+'\n' for i,s in enumerate(sequences)))

        def run(name, options=(), success=True, linked=False):
            out = root/name
            proc = subprocess.run([str(binary), '--cycles', str(cycles), '--output-folder', str(out),
                '--threads', '1', '--ram', '1G', '--read-linked-groups', str(linked).lower(), *options],
                capture_output=True, text=True, timeout=60)
            assert (proc.returncode == 0) == success, proc.stdout+proc.stderr
            return out

        def inventory(out):
            return {canonical(row['spacer']) for row in table(out/'supported_entries.tsv')}

        write_reads([prefix+repeat+s+repeat+'TT' for s in spacers for prefix in ('A','CCC')])
        options = ['--grouping-reads', str(reads)]
        baseline = run('baseline', options)
        enabled = run('enabled', [*options, '--read-recruitment', 'true'])
        assert inventory(baseline) == {canonical(s) for s in spacers[:10]}
        assert inventory(enabled) == {canonical(s) for s in spacers[:10]}
        assert {canonical(r['spacer']) for r in table(enabled/'entries.tsv')} == truth
        old = {(r['array_id'],r['entry_index']):(r['repeat_variant'],r['spacer']) for r in table(baseline/'entries.tsv')}
        new = {(r['array_id'],r['entry_index']):(r['repeat_variant'],r['spacer']) for r in table(enabled/'entries.tsv')}
        assert all(new[k] == v for k,v in old.items())
        proofs = table(enabled/'read_recruitment.tsv')
        candidates = table(enabled/'read_recruitment_candidates.tsv')
        assert any(r['added_row']=='1' for r in proofs)
        assert all(r['selection_eligible']=='0' and r['selection_status']=='novel_boundary_unresolved'
                   for r in proofs if r['added_row']=='1')
        snapshot = table(enabled/'read_recruitment_original_rows.tsv')
        assert {(r['array_id'],r['entry_index']):(r['repeat_variant'],r['spacer']) for r in snapshot} == old
        assert table(enabled/'read_recruitment_context.tsv')
        for proof in proofs:
            candidate = candidates[int(proof['candidate_id'])-1]
            assert candidate['pair_only']=='1' and candidate['accepted']=='1'
            assert candidate['source_candidate_ids']
        manifest = json.loads((enabled/'run_manifest.json').read_text())
        assert manifest['settings']['read_recruitment'] is True
        assert manifest['read_recruitment_interpretation']['legacy_context_measurements_preserved'] is True
        interpretation = manifest['read_recruitment_interpretation']
        assert interpretation['schema_version'] == 4
        policy = interpretation['row_promotion_model']
        assert policy['scope'] == 'pre_recruitment_exact_array_oriented_repeat_spacer_rows'
        assert policy['owner_scope'] == 'every_declared_candidate_owner'
        assert policy['candidate_evidence'] == 'schema3_tandem_and_error_eligibility_unchanged'
        assert policy['recursive_context_pruning'] is False and policy['minimum_spacer_distance'] == .18
        assert policy['read_linkage'] == 'inferred_before_recruitment'
        assert policy['supplemental_views'] == 'recruitment_added_rows_excluded'
        assert policy['family_context'] == 'recruitment_added_rows_cannot_supply_repeat_alleles'
        file_inventory = {row['path']: row for row in manifest['files']}
        for name in ('read_recruitment_original_rows.tsv', 'read_recruitment_context.tsv'):
            data = (enabled/name).read_bytes()
            assert file_inventory[name]['bytes'] == len(data)
            assert file_inventory[name]['checksum'] == xxh64(data)
        config = root/'settings.conf'; config.write_text('read_recruitment=true\n')
        configured = run('configured', [*options,'--settings',str(config)])
        assert (configured/'supported_entries.tsv').read_bytes() == (enabled/'supported_entries.tsv').read_bytes()
        profile_off = run('profile_off', [*options, '--read-recruitment', 'true', '--call-profile', 'off'])
        off_manifest = json.loads((profile_off/'run_manifest.json').read_text())
        assert off_manifest['counts']['read_recruitment_policy_eligible_rows'] > 0
        assert not (profile_off/'supported_entries.tsv').exists()
        before = {p.name:p.read_bytes() for p in enabled.iterdir() if p.is_file()}
        reads.write_text('@broken\nACGT\n+\nI\n')
        run('enabled', [*options,'--read-recruitment','true'], success=False)
        assert before == {p.name:p.read_bytes() for p in enabled.iterdir() if p.is_file()}
        write_reads([prefix+repeat+s+repeat+'TT' for s in spacers for prefix in ('A','CCC')])
        restored = run('enabled', options)
        assert not (restored/'read_recruitment.tsv').exists()
        assert not (restored/'read_recruitment_candidates.tsv').exists()
        assert not (restored/'read_recruitment_original_rows.tsv').exists()
        assert not (restored/'read_recruitment_context.tsv').exists()
        for name in ('entries.tsv','arrays.tsv','supported_entries.tsv','call_decisions.tsv','supported_spacers.fasta'):
            assert (restored/name).read_bytes() == (baseline/name).read_bytes(), name
        run('missing', ['--read-recruitment','true'], success=False)
        run('invalid', [*options,'--read-recruitment','maybe'], success=False)

        # A novel bridge would extend the original proper-subset linked group
        # if recruitment-appended rows participated in linkage discovery.
        singles = [prefix+repeat+s+repeat+'TT' for s in spacers for prefix in ('A','CCC')]
        chains = [(0, 1, 2), (2, 10, 3)]
        write_reads(singles + [prefix+repeat+repeat.join(spacers[i] for i in chain)+repeat+'TT'
                              for chain in chains for prefix in ('A','CCC')])
        combined = [*options, '--family-repeat-context', 'true', '--experimental-alternative-views']
        combo_baseline = run('combo_baseline', combined, linked=True)
        combo_enabled = run('combo_enabled', [*combined, '--read-recruitment', 'true'], linked=True)
        assert table(combo_baseline/'read_linked_groups.tsv')
        for name in ('read_linked_groups.tsv', 'read_linked_edges.tsv', 'supported_entries.tsv'):
            assert (combo_baseline/name).read_bytes() == (combo_enabled/name).read_bytes(), name
        assert any(r['added_row']=='1' for r in table(combo_enabled/'read_recruitment.tsv'))
        for path in (combo_baseline, combo_enabled):
            assert json.loads((path/'run_manifest.json').read_text())['counts']['supplemental_arrays'] == 0

        # A novel repeat allele must not lend its prefix to old literal calls,
        # even though independently measured recruitment can support those rows.
        alternate = repeat[:5] + ('A' if repeat[5] != 'A' else 'C') + repeat[6:]
        write_reads([prefix+repeat+s+alternate[:16] for s in spacers[:10] for prefix in ('A','CCC')]
                    + [prefix+alternate+s+alternate+'TT' for s in spacers[10:] for prefix in ('A','CCC')])
        family_options = [*options, '--family-repeat-context', 'true']
        family_baseline = run('family_baseline', family_options)
        family_enabled = run('family_enabled', [*family_options, '--read-recruitment', 'true'])
        before_context = {(r['array_id'],r['entry_index']):r for r in table(family_baseline/'call_decisions.tsv')}
        after_context = {(r['array_id'],r['entry_index']):r for r in table(family_enabled/'call_decisions.tsv')}
        for key, before in before_context.items():
            for column in ('context_status','context_observability','context_sequence_units','context_start_offsets'):
                assert after_context[key][column] == before[column], (key, column)
        assert any(r['added_row']=='1' and r['repeat_variant']==alternate
                   for r in table(family_enabled/'read_recruitment.tsv'))

        cycles.write_text(''.join(repeat+s+repeat[:22]+'\n' for s in spacers))
        write_reads([repeat[-n:]+s+repeat[:n] for s in spacers for n in (16,18)])
        partial_baseline = run('partial_baseline', options)
        partial = run('partial', [*options,'--read-recruitment','true'])
        assert not inventory(partial_baseline)
        assert inventory(partial) == truth
        decisions = table(partial/'call_decisions.tsv')
        selected = [r for r in decisions if r['selected']=='1']
        assert selected and all(r['selected_evidence']=='read_recruitment' for r in selected)
        assert all(r['context_status']!='supported' for r in selected)
        print(json.dumps(dict(discovery_baseline=10, discovery_selected=10, discovery_inclusive=12,
                              clipped_baseline=0, clipped_recruited=12, transaction_preserved=True)))


if __name__ == '__main__':
    main()
