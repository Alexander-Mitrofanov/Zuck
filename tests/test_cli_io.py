"""Check that replay input failures cannot masquerade as successful runs."""
import subprocess
import sys
import tempfile
import random
from pathlib import Path


def main() -> None:
    binary = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="my-mcaat-io-") as scratch:
        root = Path(scratch)
        output = root / "output"
        output.mkdir()
        existing = output / "CRISPR_Arrays_1.txt"
        existing.write_text("previous result\n")
        for bad_input in (root / "missing-cycles.txt", root):
            settings = root / "settings.txt"
            settings.write_text(
                f"cycles={bad_input}\noutput-folder={output}\nthreads=1\nram=1G\n"
            )
            run = subprocess.run(
                [str(binary), "--settings", str(settings)],
                capture_output=True, text=True, timeout=30,
            )
            if run.returncode == 0:
                raise AssertionError(f"invalid cycle input succeeded: {bad_input}")
            if "cycles file" not in run.stderr:
                raise AssertionError(f"missing read error: {run.stderr}")
            if existing.read_text() != "previous result\n":
                raise AssertionError("existing output changed after input error")

        # Produce a tiny cycle corpus with a shared repeat and diverse spacers.
        # The output path is a directory, so opening the result must fail.
        rng = random.Random(8472)
        sequence = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        repeat = sequence(30)
        cycles = root / "cycles.txt"
        cycles.write_text("".join(repeat + sequence(33) + repeat[:22] + "\n" for _ in range(5)))
        for stdin in ("", "y\n"):
            run = subprocess.run(
                [str(binary), "--cycles", str(cycles), "--output-folder", str(output),
                 "--threads", "999999", "--ram", "1G"],
                input=stdin, capture_output=True, text=True, timeout=30,
            )
            if run.returncode == 0 or "Invalid settings" not in run.stderr:
                raise AssertionError(f"invalid settings were not rejected: {run.stderr}")
            if existing.read_text() != "previous result\n":
                raise AssertionError("validation failure removed existing output")
            if "(y/n)" in run.stdout:
                raise AssertionError("validation failure must not offer to delete output")
        existing.unlink()
        existing.mkdir()
        run = subprocess.run(
            [str(binary), "--cycles", str(cycles), "--output-folder", str(output),
             "--threads", "1", "--ram", "1G"],
            capture_output=True, text=True, timeout=30,
        )
        if run.returncode == 0 or "arrays file" not in run.stderr:
            raise AssertionError(f"output obstruction was not reported: {run.returncode}\n{run.stdout}\n{run.stderr}")
        if Path("/dev/full").exists():
            existing.rmdir()
            existing.symlink_to("/dev/full")
            run = subprocess.run(
                [str(binary), "--cycles", str(cycles), "--output-folder", str(output),
                 "--threads", "1", "--ram", "1G"],
                capture_output=True, text=True, timeout=30,
            )
            if run.returncode == 0 or "arrays file" not in run.stderr:
                raise AssertionError(f"full-device write was not reported: {run.returncode}\n{run.stdout}\n{run.stderr}")

        # The cycle writer runs after graph construction, including when zero
        # cycles were found; it must surface an unopenable output path.
        source = Path(__file__).resolve().parents[1] / "libs/megahit/test_data"
        graph_output = root / "graph-output"
        graph_output.mkdir()
        (graph_output / "cycles.txt").mkdir()
        run = subprocess.run(
            [str(binary), "--input-files", str(source / "r3_1.fa"), str(source / "r3_2.fa"),
             "--output-folder", str(graph_output), "--threads", "1", "--ram", "1G"],
            capture_output=True, text=True, timeout=60,
        )
        if run.returncode == 0 or "Cannot write cycles" not in run.stderr:
            raise AssertionError(f"cycle writer obstruction was not reported: {run.returncode}\n{run.stdout}\n{run.stderr}")
    print("CLI I/O and validation error checks passed")


if __name__ == "__main__":
    main()
