"""Opt-in repeat-mixture reconstruction, bounded metadata and default isolation."""
import csv
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def canonical(s):
    return min(s, s.translate(str.maketrans('ACGT', 'TGCA'))[::-1])


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    binary = Path(sys.argv[1]).resolve()
    rng = random.Random(961827)
    dna = lambda n: ''.join(rng.choice('ACGT') for _ in range(n))
    common = dna(32)
    # Differences span the full five-base block (no shared terminal suffix).
    repeats = ['AACCG' + common, 'GATCA' + common]
    pairs = [(r, 'ACGT'[i % 4] + dna(30) + 'TGCA'[(i // 4 + i) % 4])
             for r in repeats for i in range(12)]
    truth = {canonical(s) for r, s in pairs}
    with tempfile.TemporaryDirectory(prefix='zuck-repeat-mixtures-') as temporary:
        root = Path(temporary)
        source = root / 'cycles.txt'
        source.write_text(''.join(r + s + r[:22] + '\n' for r, s in pairs))

        def run(name, flags=(), success=True):
            out = root / name
            proc = subprocess.run([str(binary), '--cycles', str(source), '--output-folder', str(out),
                '--threads', '1', '--ram', '1G', '--confidence-scores', 'false', '--call-profile', 'off',
                '--read-linked-groups', 'false', *flags], capture_output=True, text=True, timeout=45)
            assert (proc.returncode == 0) == success, proc.stdout + proc.stderr
            return out, proc

        base, _ = run('default', ['--boundary-mode', 'branch-aware'])
        off, _ = run('explicit_off', ['--boundary-mode', 'branch-aware', '--repeat-mixture-refinement', 'false'])
        refined, proc = run('refined', ['--boundary-mode', 'branch-aware', '--repeat-mixture-refinement', 'true'])
        assert (base / 'entries.tsv').read_bytes() == (off / 'entries.tsv').read_bytes()
        baseline = {canonical(row['spacer']) for row in table(base / 'entries.tsv')}
        observed = {canonical(row['spacer']) for row in table(refined / 'entries.tsv')}
        assert not baseline & truth
        assert observed == truth, (truth - observed, observed - truth)
        manifest = json.loads((refined / 'run_manifest.json').read_text())
        assert manifest['settings']['repeat_mixture_refinement'] is True
        assert manifest['repeat_mixture_interpretation']['maximum_partition_levels'] == 2
        assert manifest['repeat_mixture_interpretation']['reference_sequences_used'] is False
        assert 'Experimental repeat-mixture refinement' in proc.stdout
        assert '# Mixtures  : nested repeat refinement' in (refined / 'CRISPR_Arrays_1.txt').read_text()
        _, failed = run('incompatible', ['--repeat-mixture-refinement', 'true'], success=False)
        assert 'requires --boundary-mode branch-aware' in failed.stderr
        config = root / 'settings.conf'
        config.write_text('boundary_mode=branch-aware\nrepeat_mixture_refinement=true\n')
        configured, _ = run('configured', ['--settings', str(config)])
        assert (configured / 'entries.tsv').read_bytes() == (refined / 'entries.tsv').read_bytes()

        # Identical observations are one circular core, regardless of depth;
        # they cannot create the required distinct spacer contexts.
        r, s = pairs[0]
        source.write_text((r + s + r[:22] + '\n') * 100)
        negative, _ = run('duplicates', ['--boundary-mode', 'branch-aware', '--repeat-mixture-refinement', 'true'])
        assert not table(negative / 'entries.tsv')
        print(json.dumps({'fixture': 'two_repeat_alleles_with_five_base_block', 'truth': len(truth),
                          'baseline_exact': len(baseline & truth), 'refined_exact': len(observed & truth),
                          'refined_extra': len(observed - truth), 'duplicate_core_calls': 0}))


if __name__ == '__main__':
    main()
