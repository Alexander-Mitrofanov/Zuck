"""CLI boundary regression checks with synthetic generating sequences.

These test exact repeat/spacer rows, not biological specificity or universal
invariance to phase, duplication, or the graph traversal schedule.
"""
from pathlib import Path
import random
import subprocess
import sys
import tempfile

DNA = "ACGT"
COMPLEMENT = str.maketrans(DNA, "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def pair_key(repeat, spacer):
    return min((repeat, spacer), (reverse_complement(repeat), reverse_complement(spacer)))


def output_pairs(directory):
    pairs = set()
    for filename in directory.glob("CRISPR_Arrays_*.txt"):
        consensus = None
        for line in filename.read_text().splitlines():
            if line.startswith(">"):
                consensus = None
            elif line and not line.startswith("#") and consensus is None:
                consensus = line.strip()
            elif line.startswith("        "):
                repeat, spacer = line.strip().split("\t")
                if set(repeat) == {"-"}:
                    repeat = consensus
                pairs.add(pair_key(repeat, spacer))
    return pairs


def main():
    binary = Path(sys.argv[1]).resolve()
    rng = random.Random(173924)
    sequence = lambda length: "".join(rng.choice(DNA) for _ in range(length))
    repeat = sequence(45)
    observations = []
    for index in range(12):
        spacer = DNA[index % 4] + sequence(31) + DNA[(index // 4 + index) % 4]
        allele = repeat
        if index < 3:
            allele = DNA[(DNA.index(allele[0]) + 1) % 4] + allele[1:]
        if index == 7:
            allele = allele[:-1] + DNA[(DNA.index(allele[-1]) + 2) % 4]
        observations.append((allele, spacer))
    truth = {pair_key(repeat, spacer) for repeat, spacer in observations}

    with tempfile.TemporaryDirectory(prefix="my-mcaat-boundaries-") as temporary:
        root = Path(temporary)
        runs = 0

        def run(lines, options=(), expected_success=True):
            nonlocal runs
            runs += 1
            source = root / f"cycles-{runs}.txt"
            source.write_text("\n".join(lines) + "\n")
            output = root / f"output-{runs}"
            process = subprocess.run(
                [str(binary), "--cycles", str(source), "--output-folder", str(output),
                 "--threads", "1", "--ram", "1G", *options],
                capture_output=True, text=True, timeout=30,
            )
            if (process.returncode == 0) != expected_success:
                raise AssertionError(f"unexpected exit {process.returncode}:\n{process.stdout}\n{process.stderr}")
            return output_pairs(output), process

        anchored = [r + s + (r + s)[:22] for r, s in observations]
        transformed = []
        for index, (r, s) in enumerate(observations):
            core = r + s
            if index % 2:
                core = reverse_complement(core)
            phase = (index * 7 + 13) % len(core)
            core = core[phase:] + core[:phase]
            transformed.extend([core + core[:22]] * (1 + index * 3 % 7))
        # Expected truth comes from the generating pairs, independently of the
        # caller's inferred consensus, orientation, grouping, and output order.
        for lines in (anchored, transformed, list(reversed(transformed))):
            actual, _ = run(lines, ["--boundary-mode", "two-sided"])
            if actual != truth:
                raise AssertionError(f"exact pairs differ: missing={truth - actual}, extra={actual - truth}")

        # The retained legacy k-mer format must preserve the same closed cycle.
        legacy = [" ".join(line[i:i + 23] for i in range(len(line) - 22)) for line in anchored]
        actual, _ = run(legacy, ["--boundary-mode", "two-sided"])
        if actual != truth:
            raise AssertionError("legacy and compact cycles disagree on exact generating pairs")

        for cores in ([sequence(63) for _ in range(12)], [sequence(63)] * 40):
            actual, _ = run([core + core[:22] for core in cores], ["--boundary-mode", "two-sided"])
            if actual:
                raise AssertionError("unrelated or duplicate-only cycles emitted an array")
        template = sequence(63)
        noisy = [template[:i] + DNA[(DNA.index(template[i]) + 1) % 4] + template[i + 1:]
                 for i in range(43, 55)]
        actual, _ = run([core + core[:22] for core in noisy], ["--boundary-mode", "two-sided"])
        if actual:
            raise AssertionError("single-circle substitution-error control emitted an array")

        malformed = anchored[0][:-1] + DNA[(DNA.index(anchored[0][-1]) + 1) % 4]
        _, process = run([malformed], ["--boundary-mode", "two-sided"], expected_success=False)
        if "clos" not in process.stderr.lower():
            raise AssertionError(f"malformed closure lacked diagnostic: {process.stderr}")
        _, process = run(anchored, ["--boundary-mode", "unknown"], expected_success=False)
        if "boundary" not in process.stderr.lower():
            raise AssertionError(f"invalid boundary mode lacked diagnostic: {process.stderr}")

        defaults, _ = run(anchored)
        explicit_legacy, _ = run(anchored, ["--boundary-mode", "legacy"])
        if defaults != explicit_legacy:
            raise AssertionError("default boundary behavior changed")
        settings = root / "settings.txt"
        settings.write_text("boundary-mode=legacy\n")
        actual, _ = run(anchored, ["--settings", str(settings), "--boundary-mode", "two-sided"])
        if actual != truth:
            raise AssertionError("CLI boundary mode failed to override settings")
        settings.write_text("boundary-mode=two-sided\n")
        actual, _ = run(anchored, ["--settings", str(settings), "--boundary-mode", "legacy"])
        if actual != explicit_legacy:
            raise AssertionError("CLI legacy mode failed to override two-sided settings")
    print("Two-sided exact boundary and CLI regression checks passed")


if __name__ == "__main__":
    main()
