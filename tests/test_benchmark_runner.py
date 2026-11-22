"""Check benchmark CLI wiring and failure handling without building a graph."""
import csv
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    runner = Path(__file__).resolve().parents[1] / "bench/run_benchmarks.sh"
    with tempfile.TemporaryDirectory(prefix="zuck-benchmark-test-") as temporary:
        root = Path(temporary)
        graph = root / "graph with spaces"
        graph.mkdir()
        binary = root / "fake Zuck"
        binary.write_text('''#!/usr/bin/env python3
import json
from pathlib import Path
import sys
args = sys.argv[1:]
assert args[0] == '--graph' and args[2] == '--threads' and args[4] == '--output-folder', args
assert Path(args[1]).is_dir()
output = Path(args[5])
output.mkdir()
(output / 'run_manifest.json').write_text(json.dumps({'counts': {'loaded_cycle_records': 7}}))
''')
        binary.chmod(0o755)
        output = root / "result table.csv"
        base = ["bash", str(runner), "--binary", str(binary), "--input", str(graph),
                "--threads", "1,2", "--runs", "2", "-o", str(output)]
        result = subprocess.run(base, capture_output=True, text=True, timeout=30)
        assert result.returncode == 0, result.stdout + result.stderr
        with output.open() as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == 4 and all(row["cycles"] == "7" for row in rows)
        assert [(row["threads"], row["run"]) for row in rows] == [("1", "1"), ("1", "2"), ("2", "1"), ("2", "2")]
        for args in (["--binary"], ["--help"], [*base[2:], "--runs", "0"], [*base[2:], "--runs", "9223372036854775808"], [*base[2:], "--threads", "1,,2"]):
            result = subprocess.run(["bash", str(runner), *args], capture_output=True, text=True)
            assert (result.returncode == 0) == (args == ["--help"]), result
        binary.write_text("#!/bin/sh\necho deliberate-failure >&2\nexit 7\n")
        result = subprocess.run(base, capture_output=True, text=True, timeout=30)
        assert result.returncode != 0 and "deliberate-failure" in result.stderr
    print("Benchmark runner validates options and invokes the real Zuck interface")


if __name__ == "__main__":
    main()
