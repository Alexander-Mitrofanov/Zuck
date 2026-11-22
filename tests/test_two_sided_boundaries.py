"""CLI boundary regression checks with synthetic generating sequences.

These test exact repeat/spacer rows, not biological specificity or universal
invariance to phase, duplication, or the graph traversal schedule.
"""
from pathlib import Path
import json
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
    help_result = subprocess.run([str(binary), "--help"], capture_output=True,
                                 text=True, timeout=30, check=True)
    assert "Repeat boundary inference [default: two-sided]" in help_result.stdout
    assert "--short-spacer-filter <true|false>" in help_result.stdout
    assert "Reject arrays explained by one short spacer [default: false]" in help_result.stdout
    assert "two-sided and branch-aware are experimental" not in help_result.stdout
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

    with tempfile.TemporaryDirectory(prefix="zuck-boundaries-") as temporary:
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
            for mode in ("two-sided", "branch-aware"):
                actual, _ = run(lines, ["--boundary-mode", mode])
                if actual != truth:
                    raise AssertionError(f"{mode} exact pairs differ: missing={truth - actual}, extra={actual - truth}")

        # Two substantially divergent repeats can share the exact seed. A
        # single group consensus truncates them; linked branch evidence must
        # recover the generating tuples, never invent cross-family pairings.
        branch_rng = random.Random(173924)
        branch_dna = lambda length: "".join(branch_rng.choice(DNA) for _ in range(length))
        common = branch_dna(23)
        suffix = branch_dna(17)
        repeats = [common + suffix,
                   common + "".join(DNA[(DNA.index(base) + 1) % 4] for base in suffix)]
        branch_rows = [(r, DNA[i % 4] + branch_dna(31) + DNA[(i // 4 + i) % 4])
                       for r in repeats for i in range(12)]
        branch_truth = {pair_key(r, s) for r, s in branch_rows}
        branch_lines = [r + s + (r + s)[:22] for r, s in branch_rows]
        old_pairs, _ = run(branch_lines, ["--boundary-mode", "two-sided"])
        if old_pairs == branch_truth:
            raise AssertionError("branch fixture no longer distinguishes single-group extraction")
        phased_lines = []
        for index, (r, s) in enumerate(branch_rows):
            core = r + s
            for strand in (core, reverse_complement(core)):
                for phase in (0, 7, 29, 51):
                    rotated = strand[phase:] + strand[:phase]
                    phased_lines.extend([rotated + rotated[:22]] * (1 + index % 4))
        for lines in (branch_lines, phased_lines, list(reversed(phased_lines))):
            actual, _ = run(lines, ["--boundary-mode", "branch-aware"])
            if actual != branch_truth:
                raise AssertionError(f"branch recovery differs: missing={branch_truth - actual}, extra={actual - branch_truth}")
            for path in (root / f"output-{runs}").glob("CRISPR_Arrays_*.txt"):
                if "# Boundary  : branch-aware (experimental)" not in path.read_text():
                    raise AssertionError("branch-aware output lacks experimental mode marker")

        # Regression for data-dependent overextension: a six-base signature
        # may be five common repeat bases plus ONE random spacer endpoint.
        # Naive branch partitioning recovered only 1/40 of these exact pairs.
        fixture_rng = random.Random(417209)
        fixture_dna = lambda length: "".join(fixture_rng.choice(DNA) for _ in range(length))
        fixture_repeat = fixture_dna(28)
        fixture_rows = []
        for i in range(40):
            spacer = DNA[i % 4] + fixture_dna(31) + DNA[(i // 4 + i) % 4]
            allele = fixture_repeat
            if i == 0:
                allele = allele[:4] + DNA[(DNA.index(allele[4]) + 1) % 4] + allele[5:]
            fixture_rows.append((allele, spacer))
        fixture_truth = {pair_key(r, s) for r, s in fixture_rows}
        fixture_lines = [r + s + (r + s)[:22] for r, s in fixture_rows]
        for mode in ("two-sided", "branch-aware"):
            actual, _ = run(fixture_lines, ["--boundary-mode", mode])
            if actual != fixture_truth:
                raise AssertionError(f"{mode} selected spacer endpoints as repeat branches")

        # Fewer than three distinct cores in either branch cannot establish
        # the new branch model. Duplicated observations are not new support.
        insufficient = branch_lines[:3] + branch_lines[12:14]
        before, _ = run(insufficient, ["--boundary-mode", "two-sided"])
        after, _ = run(insufficient * 10, ["--boundary-mode", "branch-aware"])
        if after != before:
            raise AssertionError("duplicate records created unsupported branch evidence")

        # The retained legacy k-mer format must preserve the same closed cycle.
        legacy = [" ".join(line[i:i + 23] for i in range(len(line) - 22)) for line in anchored]
        actual, _ = run(legacy, ["--boundary-mode", "two-sided"])
        if actual != truth:
            raise AssertionError("legacy and compact cycles disagree on exact generating pairs")

        for cores in ([sequence(63) for _ in range(12)], [sequence(63)] * 40):
            for mode in ("two-sided", "branch-aware"):
                actual, _ = run([core + core[:22] for core in cores], ["--boundary-mode", mode])
                if actual:
                    raise AssertionError("unrelated or duplicate-only cycles emitted an array")
        template = sequence(63)
        noisy = [template[:i] + DNA[(DNA.index(template[i]) + 1) % 4] + template[i + 1:]
                 for i in range(43, 55)]
        for mode in ("two-sided", "branch-aware"):
            actual, _ = run([core + core[:22] for core in noisy], ["--boundary-mode", mode])
            if actual:
                raise AssertionError("single-circle substitution-error control emitted an array")

        malformed = anchored[0][:-1] + DNA[(DNA.index(anchored[0][-1]) + 1) % 4]
        _, process = run([malformed], ["--boundary-mode", "two-sided"], expected_success=False)
        if "clos" not in process.stderr.lower():
            raise AssertionError(f"malformed closure lacked diagnostic: {process.stderr}")
        _, process = run(anchored, ["--boundary-mode", "unknown"], expected_success=False)
        if "boundary" not in process.stderr.lower():
            raise AssertionError(f"invalid boundary mode lacked diagnostic: {process.stderr}")

        defaults, default_process = run(anchored)
        default_output = root / f"output-{runs}"
        default_manifest = json.loads((default_output / "run_manifest.json").read_text())
        assert default_manifest["tool"] == "Zuck"
        assert default_manifest["settings"]["boundary_mode"] == "two-sided"
        assert default_manifest["settings"]["short_spacer_filter"] is False
        assert default_manifest["settings"]["short_spacer_filter_applied"] is False
        assert "two-sided (experimental)" not in default_process.stdout
        for path in default_output.glob("CRISPR_Arrays_*.txt"):
            assert "# Boundary  : two-sided\n" in path.read_text()
            assert "# Boundary  : two-sided (experimental)" not in path.read_text()
        explicit_two_sided, _ = run(anchored, ["--boundary-mode", "two-sided"])
        explicit_output = root / f"output-{runs}"
        if defaults != explicit_two_sided or defaults != truth:
            raise AssertionError("default is not the selected two-sided behavior")
        for name in ("arrays.tsv", "entries.tsv"):
            assert (default_output / name).read_bytes() == (explicit_output / name).read_bytes()
        default_branches, _ = run(branch_lines)
        if default_branches != old_pairs:
            raise AssertionError("default unexpectedly enabled branch-aware inference")
        explicit_legacy, _ = run(anchored, ["--boundary-mode", "legacy"])
        if defaults == explicit_legacy:
            raise AssertionError("default fixture must distinguish two-sided and legacy")
        settings = root / "settings.txt"
        settings.write_text("boundary-mode=legacy\n")
        actual, _ = run(anchored, ["--settings", str(settings)])
        if actual != explicit_legacy:
            raise AssertionError("legacy settings failed to override the new default")
        actual, _ = run(anchored, ["--settings", str(settings), "--boundary-mode", "two-sided"])
        if actual != truth:
            raise AssertionError("CLI boundary mode failed to override settings")
        settings.write_text("boundary-mode=two-sided\n")
        actual, _ = run(anchored, ["--settings", str(settings)])
        if actual != defaults:
            raise AssertionError("two-sided settings differ from the default")
        actual, _ = run(anchored, ["--settings", str(settings), "--boundary-mode", "legacy"])
        if actual != explicit_legacy:
            raise AssertionError("CLI legacy mode failed to override two-sided settings")
        settings.write_text("boundary-mode=branch-aware\n")
        actual, _ = run(branch_lines, ["--settings", str(settings)])
        if actual != branch_truth:
            raise AssertionError("settings file failed to enable branch-aware inference")
        actual, _ = run(branch_lines, ["--settings", str(settings), "--boundary-mode", "two-sided"])
        if actual != old_pairs:
            raise AssertionError("CLI two-sided mode failed to disable branch inference")
        actual, _ = run(anchored, ["--settings", str(settings), "--boundary-mode", "legacy"])
        if actual != explicit_legacy:
            raise AssertionError("CLI legacy mode failed to disable both modern modes")

        # Settings are applied before CLI options, independent of their order.
        settings.write_text("boundary-mode=legacy\nshort-spacer-filter=false\n")
        for options in (["--settings", str(settings)],
                        ["--settings", str(settings), "--short-spacer-filter", "true"],
                        ["--short-spacer-filter", "true", "--settings", str(settings)]):
            run(anchored, options)
            manifest = json.loads((root / f"output-{runs}" / "run_manifest.json").read_text())
            assert manifest["settings"]["boundary_mode"] == "legacy"
            assert manifest["settings"]["short_spacer_filter"] == ("true" in options)
        actual, _ = run(anchored, ["--boundary-mode", "two-sided", "--settings", str(settings)])
        if actual != truth:
            raise AssertionError("CLI mode preceding settings did not override the file")
        settings.write_text("short-spacer-filter=true\n")
        run(anchored, ["--short-spacer-filter", "false", "--settings", str(settings)])
        manifest = json.loads((root / f"output-{runs}" / "run_manifest.json").read_text())
        assert manifest["settings"]["boundary_mode"] == "two-sided"
        assert manifest["settings"]["short_spacer_filter"] is False
        _, process = run(anchored, ["--short-spacer-filter", "unknown"], expected_success=False)
        assert "short-spacer-filter" in process.stderr
        _, process = run(anchored, ["--short-spacer-filter"], expected_success=False)
        assert "Missing value" in process.stderr
        settings.write_text("short-spacer-filter=unknown\n")
        _, process = run(anchored, ["--settings", str(settings)], expected_success=False)
        assert "short-spacer-filter" in process.stderr
    print("Two-sided exact boundary and CLI regression checks passed")


if __name__ == "__main__":
    main()
