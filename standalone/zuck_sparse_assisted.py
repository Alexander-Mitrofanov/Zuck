#!/usr/bin/env python3
"""Optionally add independently called CRASS spacers to a low-coverage adaptive Zuck inventory."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

DNA = "ACGTRYMKBDHVSWN"
COMPLEMENT = str.maketrans(DNA, "TGCAYRKMVHDBSWN")
ZERO_MARKER = "[crass_graphBuilder]: 0 CRISPRs found!"


def canonical(seq):
    seq = seq.upper()
    if not seq or any(b not in DNA for b in seq):
        raise ValueError("Invalid nucleotide sequence in spacer inventory")
    return min(seq, seq.translate(COMPLEMENT)[::-1])


def endpoint_keys(seq):
    """Strand-independent exact interiors with up to three clipped bases at each end."""
    for left in range(4):
        for right in range(4):
            if len(seq) - left - right >= 20:
                core = seq[left:len(seq)-right]
                yield min(core, core.translate(COMPLEMENT)[::-1])


def partition(native, crass):
    retained = {canonical(s) for s in native}
    native_keys = {key for seq in retained for key in endpoint_keys(seq)}
    decisions = {}
    for seq in sorted({canonical(s) for s in crass}):
        if seq in retained:
            decisions[seq] = "exact_native"
        elif any(key in native_keys for key in endpoint_keys(seq)):
            decisions[seq] = "terminal_equivalent_native"
        else:
            decisions[seq] = "crass_added"
    return retained, decisions


def merge(native, crass):
    retained, decisions = partition(native, crass)
    return retained | {seq for seq, reason in decisions.items() if reason == "crass_added"}


def fasta(path):
    seqs, current, header = set(), [], False
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    seqs.add(canonical("".join(current)))
                header, current = True, []
            else:
                if not header:
                    raise ValueError("Invalid FASTA without header")
                current.append(line)
    if header:
        seqs.add(canonical("".join(current)))
    return seqs


def identity(path):
    digest = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            digest.update(chunk)
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def crass_calls(native, log):
    xml = sorted(native.glob("*.crispr"))
    if not xml:
        with log.open(errors="replace") as stream:
            if any(ZERO_MARKER in line for line in stream):
                return set()
        raise ValueError("CRASS omitted its XML without an explicit zero-call result")
    if len(xml) != 1:
        raise ValueError("Ambiguous CRASS XML output")
    root = ET.parse(xml[0]).getroot()
    if root.tag != "crispr":
        raise ValueError("Invalid CRASS XML root")
    return {canonical(node.attrib["seq"]) for node in root.findall("group/data/spacers/spacer")}


def plain_reads(paths, work):
    result, identities = [], []
    for i, original in enumerate(paths):
        record = identity(original)
        with original.open("rb") as probe:
            compressed = probe.read(2) == b"\x1f\x8b"
        if compressed:
            path = work / f"mate{i+1}.fastq"
            with gzip.open(original, "rb") as src, path.open("xb") as dst:
                shutil.copyfileobj(src, dst, 1 << 20)
            decoded = identity(path)
        else:
            path, decoded = original, record
        result.append(path)
        identities.append({"provided": record, "decoded_sha256": decoded["sha256"], "decoded_bytes": decoded["bytes"]})
    return result, identities


def execute(args):
    native, output = args.native.resolve(strict=True), args.output_folder.resolve()
    receipt_path = native / "coverage_routing.json"
    route = json.loads(receipt_path.read_bytes())
    details = route["routing"]
    if (details["native_exit_status"] != 0 or details["policy"] != "native_sparse_junctions_if_low_v1"
            or route["regime"] not in {"low", "intermediate", "high", "mixed", "uncertain"}):
        raise ValueError("Expected a completed coverage-routed native adaptive Zuck result")
    if not (native / "run_manifest.json").is_file():
        raise ValueError("Missing native Zuck manifest")
    inventory_file = native / "read_first_spacers.fasta"
    native_calls = fasta(inventory_file)
    recorded_argv = details["argv"]
    first = recorded_argv.index("--input-files") + 1
    stop = next((i for i in range(first,len(recorded_argv)) if recorded_argv[i].startswith("--")), len(recorded_argv))
    archived = [Path(p) for p in recorded_argv[first:stop]]
    paths = [p.resolve(strict=True) for p in (args.input_files or archived)]
    if not 1 <= len(paths) == len(archived) <= 2 or len(set(paths)) != len(paths):
        raise ValueError("Native Zuck and CRASS must process the same FASTQ input list")
    relocated = any(old.resolve() != new for old, new in zip(archived, paths))
    if relocated and not args.native_read_sha256:
        raise ValueError("Relocated inputs require the original decoded FASTQ SHA-256 values in mate order")
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise ValueError("Choose a new assisted output folder")
    with tempfile.TemporaryDirectory(dir=output.parent, prefix=".zuck-assisted-") as tmp:
        scratch = Path(tmp)
        bundle = scratch / "results"
        bundle.mkdir()
        crass, command, binary, library_dir = set(), None, None, None
        low = route["regime"] == "low"
        if low:
            delivered, identities = plain_reads(paths, scratch)
        else:
            identities = [{"provided": identity(p)} for p in paths]
            if relocated or args.native_read_sha256:
                _, identities = plain_reads(paths, scratch)
        if args.native_read_sha256 and [i["decoded_sha256"] for i in identities] != args.native_read_sha256:
            raise ValueError("Input bytes differ from the original native read bindings")
        if low:
            exe = shutil.which(args.crass)
            if not exe:
                raise ValueError("Low-depth assisted mode requires the CRASS executable")
            executable = Path(exe).resolve()
            binary = identity(executable)
            dest = bundle / "crass"
            dest.mkdir()
            command = [str(executable), "-r", "-o", str(dest), *map(str, delivered)]
            env = dict(os.environ)
            if args.crass_lib_dir:
                library_dir = str(args.crass_lib_dir.resolve(strict=True))
                env["LD_LIBRARY_PATH"] = library_dir + (":" + env["LD_LIBRARY_PATH"] if env.get("LD_LIBRARY_PATH") else "")
            log = bundle / "crass.log"
            with log.open("wb") as stream:
                subprocess.run(command, cwd=scratch, env=env, stdout=stream, stderr=subprocess.STDOUT, check=True)
            crass = crass_calls(dest, log)
        retained, decisions = partition(native_calls, crass)
        additions = {s for s in crass if decisions[s] == "crass_added"}
        with (bundle / "assisted_spacers.fasta").open("w") as f:
            for i, s in enumerate(sorted(retained | additions),1):
                f.write(f">spacer_{i} {'native' if s in retained else 'crass_added'}\n{s}\n")
        with (bundle / "crass_decisions.tsv").open("w") as f:
            f.write("canonical_spacer\tdecision\n")
            for seq, reason in decisions.items():
                f.write(f"{seq}\t{reason}\n")
        record = {"policy": "add_crass_if_low_exact_core20_terminal3", "regime": route["regime"],
                  "native_routing": identity(receipt_path), "native_inventory": identity(inventory_file),
                  "original_native_read_paths": list(map(str, archived)), "input_files": identities,
                  "native_read_sha256_supplied": args.native_read_sha256, "crass_executed": low,
                  "crass_binary": binary, "crass_library_directory": library_dir, "crass_argv": command,
                  "native_count": len(retained), "crass_count": len(crass), "added": len(additions),
                  "combined_count": len(retained | additions)}
        (bundle / "assisted.json").write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
        if output.exists():
            raise ValueError("Assisted output folder was concurrently created")
        bundle.rename(output)
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native", type=Path, required=True, help="Successful zuck_depth.py adaptive output directory.")
    p.add_argument("--output-folder", type=Path, required=True)
    p.add_argument("--crass", default="crass", help="CRASS 1.0.1 executable (only invoked on low input).")
    p.add_argument("--crass-lib-dir", type=Path, help="Shared libraries for a separately installed CRASS.")
    p.add_argument("--input-files", type=Path, nargs="+", help="Optional relocated original FASTQs in native order.")
    p.add_argument("--native-read-sha256", nargs="+", help="Decoded original FASTQ checksums, in native mate order.")
    print(json.dumps(execute(p.parse_args()), sort_keys=True))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, ET.ParseError, subprocess.CalledProcessError) as error:
        print(f"CRASS-assisted Zuck failed: {error}", file=sys.stderr)
        sys.exit(1)
