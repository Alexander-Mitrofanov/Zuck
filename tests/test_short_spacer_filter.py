"""Closed-cycle CLI checks for bounded short-spacer singleton evidence."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    binary = Path(sys.argv[1]).resolve()
    repeat = "AAAAAGAATAGGATTGTTTATTAAGGAGAAA"
    short = "GACTAGAATTATTTGTTGAG"
    with tempfile.TemporaryDirectory(prefix="zuck-short-spacers-") as temporary:
        root = Path(temporary)
        counter = 0

        def run(spacers, *options):
            nonlocal counter
            counter += 1
            source = root / f"cycles-{counter}.txt"
            source.write_text("\n".join(repeat + s + (repeat + s)[:22]
                                         for s in spacers) + "\n")
            output = root / f"output-{counter}"
            process = subprocess.run(
                [str(binary), "--cycles", str(source), "--output-folder", str(output),
                 "--threads", "1", "--ram", "1G", *options],
                capture_output=True, text=True, timeout=30)
            if process.returncode:
                raise AssertionError(process.stdout + process.stderr)
            with (output / "entries.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            manifest = json.loads((output / "run_manifest.json").read_text())
            return rows, manifest

        terminal_variants = [short, short + "A"]
        baseline, disabled = run(terminal_variants, "--short-spacer-filter", "false")
        assert len(baseline) == 2, baseline
        defaults, default_manifest = run(terminal_variants)
        assert defaults == baseline
        assert default_manifest["settings"]["short_spacer_filter"] is False
        filtered, enabled = run(terminal_variants, "--short-spacer-filter", "true")
        assert filtered == []
        assert enabled["counts"]["filtered_short_spacer_singletons"] == 1
        assert enabled["settings"]["short_spacer_filter_applied"] is True
        assert disabled["counts"]["filtered_short_spacer_singletons"] == 0
        assert disabled["settings"]["short_spacer_filter_applied"] is False
        assert enabled["postprocessing_thresholds"]["short_spacer_max_terminal_extension"] == 2

        # Short length alone is never a rejection criterion. Two independent
        # short spacers retain exactly the same observed rows.
        independent = [short, "TCAATTCGAGCTAGTACGAT"]
        before, _ = run(independent, "--short-spacer-filter", "false")
        after, _ = run(independent, "--short-spacer-filter", "true")
        assert before and after == before

        before, _ = run(terminal_variants, "--boundary-mode", "legacy",
                        "--short-spacer-filter", "false")
        after, manifest = run(terminal_variants, "--boundary-mode", "legacy",
                              "--short-spacer-filter", "true")
        assert after == before
        assert manifest["settings"]["short_spacer_filter_applied"] is False

        settings = root / "settings.txt"
        settings.write_text("short-spacer-filter=false\n")
        actual, _ = run(terminal_variants, "--settings", str(settings))
        assert actual == baseline
        actual, _ = run(terminal_variants, "--settings", str(settings),
                        "--short-spacer-filter", "true")
        assert actual == []
        print("Short-spacer filter CLI checks passed")


if __name__ == "__main__":
    main()
