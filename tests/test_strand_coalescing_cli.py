"""Exact synthetic observations survive strand-equivalent group coalescing."""
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
    rng = random.Random(173924)
    dna = "ACGT"
    reverse = lambda s: s.translate(str.maketrans(dna, "TGCA"))[::-1]
    sequence = lambda n: "".join(rng.choice(dna) for _ in range(n))
    repeat = sequence(45)
    spacers = [dna[i % 4] + sequence(31) + dna[(i // 4 + i) % 4]
               for i in range(12)]
    expected = Counter([(repeat, s) for s in spacers] +
                       [(reverse(repeat), reverse(s)) for s in spacers])
    cores = [c for s in spacers for c in (repeat + s, reverse(repeat + s))]
    with tempfile.TemporaryDirectory(prefix="zuck-strand-coalescing-") as temporary:
        root = Path(temporary)
        source = root / "cycles.txt"
        source.write_text("\n".join(c + c[:22] for c in cores) + "\n")
        output = root / "output"
        process = subprocess.run(
            [str(binary), "--cycles", str(source), "--output-folder", str(output),
             "--threads", "1", "--ram", "1G", "--boundary-mode", "two-sided"],
            capture_output=True, text=True, timeout=30)
        if process.returncode:
            raise AssertionError(process.stdout + process.stderr)
        with (output / "entries.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len({r["array_id"] for r in rows}) == 1
        assert Counter((r["repeat_variant"], r["spacer"]) for r in rows) == expected
        manifest = json.loads((output / "run_manifest.json").read_text())
        assert manifest["counts"]["coalesced_strand_equivalent_arrays"] == 1
        assert manifest["counts"]["primary_arrays"] == 1
        assert manifest["counts"]["primary_entries"] == 24

        text_pairs = []
        consensus = None
        for path in output.glob("CRISPR_Arrays_*.txt"):
            text = path.read_text()
            assert "# Rows      : observed orientations; unordered" in text
            for line in text.splitlines():
                if line.startswith(">"):
                    consensus = None
                elif line and not line.startswith("#") and consensus is None:
                    consensus = line.strip()
                elif line.startswith("        "):
                    observed_repeat, spacer = line.strip().split("\t")
                    if set(observed_repeat) == {"-"}:
                        observed_repeat = consensus
                    text_pairs.append((observed_repeat, spacer))
        assert Counter(text_pairs) == expected
        print("Strand coalescing preserves every original oriented row")


if __name__ == "__main__":
    main()
