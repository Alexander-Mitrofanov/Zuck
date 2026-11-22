"""Validate replay identity, machine-readable counts, and failure preservation."""
import csv
import json
import os
from pathlib import Path
import random
import signal
import shutil
import subprocess
import sys
import tempfile

from test_standalone_binary import PT_INTERP, program_header_types


def read_tsv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def normalized_text(directory):
    return {
        path.name: "\n".join(line for line in path.read_text().splitlines()
                             if not line.startswith("# Generated :"))
        for path in directory.glob("CRISPR_Arrays_*.txt")
    }


def main():
    binary = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="zuck-bundle-") as temporary:
        root = Path(temporary)
        rename_fault = None
        # LD_PRELOAD interposition is meaningful only for dynamically loaded
        # executables. Static release builds cover the same preservation logic
        # through destination obstruction, malformed input, and RLIMIT_FSIZE.
        dynamically_loaded = PT_INTERP in program_header_types(binary)
        if sys.platform == "linux" and dynamically_loaded and shutil.which("cc"):
            # Force the second staged-file installation to fail, after the old
            # bundle has been backed up and one new result installed. This
            # exercises rollback without racing the process or filling a disk.
            fault_source = root / "rename_fault.c"
            fault_source.write_text(r'''
#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <string.h>
int rename(const char *from, const char *to) {
    static int installs = 0;
    int (*real_rename)(const char *, const char *) = dlsym(RTLD_NEXT, "rename");
    if (strstr(from, ".Zuck-results-") && strstr(from, "/new/") && ++installs == 2) {
        errno = EIO;
        return -1;
    }
    return real_rename(from, to);
}
''')
            rename_fault = root / "rename_fault.so"
            subprocess.run(["cc", "-shared", "-fPIC", str(fault_source), "-ldl",
                            "-o", str(rename_fault)], check=True, capture_output=True)
        cycles = root / 'cycles "quoted".txt'
        rng = random.Random(173924)
        dna = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
        repeat = dna(30)
        lines = [repeat + "ACGT"[i % 4] + dna(31) + "ACGT"[(i // 4 + i) % 4]
                 for i in range(12)]
        cycles.write_text("# input comment\n\n" + "\n".join(s + s[:22] for s in lines) + "\n")

        def run(output, mode="legacy", source=cycles, **kwargs):
            return subprocess.run(
                [str(binary), "--cycles", str(source), "--output-folder", str(output),
                 "--threads", "1", "--ram", "1G", "--boundary-mode", mode],
                capture_output=True, text=True, timeout=30, **kwargs)

        for mode in ("legacy", "two-sided", "branch-aware"):
            output = root / mode
            completed = run(output, mode)
            assert completed.returncode == 0, completed.stderr
            manifest = json.loads((output / "run_manifest.json").read_text())
            arrays = read_tsv(output / "arrays.tsv")
            entries = read_tsv(output / "entries.tsv")
            assert arrays and entries
            assert manifest["schema_version"] == 1
            assert manifest["tool"] == "Zuck"
            assert len(manifest["source_sha256"]) == 64
            assert manifest["settings"]["boundary_mode"] == mode
            assert manifest["settings"]["short_spacer_filter"] is False
            assert manifest["settings"]["short_spacer_filter_applied"] is False
            assert manifest["settings"]["directed_primary_groups"] is False
            assert manifest["directed_primary_group_interpretation"]["enabled"] is False
            assert manifest["settings"]["cross_variant_scaffolds"] is False
            assert manifest["cross_variant_scaffold_interpretation"]["enabled"] is False
            assert manifest["settings"]["circular_unit_dedup"] is False
            assert manifest["supported_call_interpretation"]["circular_unit_dedup"]["enabled"] is False
            assert manifest["cycle_input"]["path"] == str(cycles)
            assert manifest["cycle_input"]["bytes"] == cycles.stat().st_size
            assert not manifest["cycle_search_settings"]["applied_in_this_run"]
            assert manifest["counts"]["loaded_cycle_records"] == 12
            assert manifest["counts"]["primary_arrays"] == len(arrays)
            assert manifest["counts"]["primary_entries"] == len(entries)
            assert manifest["counts"]["supplemental_arrays"] == 0
            assert manifest["interpretation"]["spacer_order"] == "unspecified"
            assert manifest["scoring"]["kind"] == "heuristic_evidence_indices"
            assert manifest["scoring"]["calibrated_probability"] is False
            assert manifest["scoring"]["primary_rows_preserved"] is True
            inventory = {item["path"]: item for item in manifest["files"]}
            expected_inventory = {"arrays.tsv", "entries.tsv", "element_confidence.tsv",
                                  "call_decisions.tsv", "supported_entries.tsv", "supported_arrays.tsv",
                                  "supported_spacers.fasta", "repeat_haplotype_groups.tsv",
                                  *normalized_text(output)}
            assert set(inventory) == expected_inventory, {
                "missing": sorted(expected_inventory - set(inventory)),
                "unexpected": sorted(set(inventory) - expected_inventory),
            }
            for name, item in inventory.items():
                assert item["bytes"] == (output / name).stat().st_size
                assert len(item["checksum"]) == 16

            # Reconstruct the text view independently and compare every row,
            # including literal repeat alleles hidden behind consensus dashes.
            text_rows = []
            for path in sorted(output.glob("CRISPR_Arrays_*.txt")):
                array_id = consensus = None
                index = 0
                for line in path.read_text().splitlines():
                    if line.startswith(">"):
                        array_id = line.split()[0][1:]
                        consensus = None
                        index = 0
                    elif line.startswith("        "):
                        allele, spacer = line.strip().split("\t")
                        index += 1
                        text_rows.append((array_id, str(index), consensus if set(allele) == {"-"} else allele, spacer))
                    elif line and not line.startswith("#") and consensus is None:
                        consensus = line
            assert text_rows == [(e["array_id"], e["entry_index"], e["repeat_variant"], e["spacer"]) for e in entries]
            for array in arrays:
                rows = [e for e in entries if e["array_id"] == array["array_id"]]
                assert int(array["entry_count"]) == len(rows)
                assert array["kind"] == "primary" and array["parent_array_id"] == ""

            before = {p.name: p.read_bytes() for p in output.iterdir()}
            another = root / (mode + "-again")
            assert run(another, mode).returncode == 0
            assert normalized_text(output) == normalized_text(another)
            for name in ("arrays.tsv", "entries.tsv"):
                assert (output / name).read_bytes() == (another / name).read_bytes()
            again = json.loads((another / "run_manifest.json").read_text())
            assert again["cycle_input"] == manifest["cycle_input"]

            # A late destination obstruction must preserve every old file.
            (output / "run_manifest.json").unlink()
            (output / "run_manifest.json").mkdir()
            failed = run(output, mode)
            assert failed.returncode != 0, failed.stdout
            for name, content in before.items():
                if name != "run_manifest.json":
                    assert (output / name).read_bytes() == content
            (output / "run_manifest.json").rmdir()
            (output / "run_manifest.json").write_bytes(before["run_manifest.json"])

            # A malformed late record must not replace any earlier results.
            invalid = root / "invalid.txt"
            invalid.write_bytes(cycles.read_bytes() + b"INVALID-DNA\n")
            assert run(output, mode, invalid).returncode != 0
            assert before == {p.name: p.read_bytes() for p in output.iterdir()}

            # Simulate disk exhaustion after opening a staged output. Ignoring
            # SIGXFSZ lets the normal stream error/rollback path run on POSIX.
            if sys.platform == "linux":
                import resource

                def limit_file_size():
                    signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
                    resource.setrlimit(resource.RLIMIT_FSIZE, (256, 256))

                failed = run(output, mode, preexec_fn=limit_file_size)
                assert failed.returncode != 0
                assert "Error writing" in failed.stderr, failed.stderr
                assert before == {p.name: p.read_bytes() for p in output.iterdir()}

            if rename_fault:
                failed = run(output, mode, env={**os.environ, "LD_PRELOAD": str(rename_fault)})
                assert failed.returncode != 0
                assert "previous results restored" in failed.stderr, failed.stderr
                assert before == {p.name: p.read_bytes() for p in output.iterdir()}

            # Empty success removes stale numbered outputs, retaining unrelated
            # files and producing a zero-count manifest and header-only tables.
            (output / "CRISPR_Arrays_999.txt").write_text("stale\n")
            unrelated = output / "CRISPR_Arrays_notes.txt"
            unrelated.write_text("keep\n")
            empty = root / "empty.txt"
            empty.write_text("")
            assert run(output, mode, empty).returncode == 0
            assert unrelated.read_text() == "keep\n"
            assert not (output / "CRISPR_Arrays_999.txt").exists()
            zero = json.loads((output / "run_manifest.json").read_text())
            assert zero["cycle_input"]["checksum"] == "ef46db3751d8e999"
            assert zero["counts"]["primary_arrays"] == 0
            assert zero["counts"]["primary_entries"] == 0
            assert not read_tsv(output / "entries.tsv")
            assert not read_tsv(output / "arrays.tsv")
            assert not list(output.glob(".Zuck-results-*"))

        # An input alias must never be deleted during successful replacement.
        output = root / "alias"
        output.mkdir()
        alias = output / "entries.tsv"
        alias.write_bytes(cycles.read_bytes())
        assert run(output, source=alias).returncode != 0
        assert alias.read_bytes() == cycles.read_bytes()
    print("Result bundle identity, parity, and preservation checks passed")


if __name__ == "__main__":
    main()
