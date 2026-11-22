#!/usr/bin/env python3
"""Build a bounded, confidence-preserving CRISPRidentify catalogue projection.

Only cycle-derived sequence seeds select arrays. Every source shard is verified;
reference annotations are never accepted as input. The resulting TSV retains
whole-array co-occurrence, and native Zuck verifies complete sequence matches.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import sqlite3
import time

HEADER = "array_id\taccession\tcategory\tconfidence\trepeat_count\tspacer_count\trepeats\tspacers\n"
DNA = frozenset("ACGT")
MAX_SHARD_BYTES = 128 * 1024 * 1024
MAX_COMPRESSED_BYTES = 64 * 1024 * 1024
_SEEDS = frozenset()
_ROOT = None


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def identity(path):
    path = Path(path).resolve()
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def glue_cycle(line):
    kmers = line.upper().split()
    if not kmers:
        return ""
    sequence = kmers[0]
    for kmer in kmers[1:]:
        overlap = next((n for n in range(min(len(sequence), len(kmer)), 0, -1)
                        if sequence[-n:] == kmer[:n]), 0)
        sequence += kmer[overlap:]
    return sequence


def cycle_seeds(path):
    seeds = set()
    records = 0
    with Path(path).open() as handle:
        for number, line in enumerate(handle, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            sequence = glue_cycle(line)
            if set(sequence) - set("ACGTN"):
                raise ValueError(f"Invalid cycle DNA at line {number}")
            if len(sequence) < 44 or sequence[:22] != sequence[-22:]:
                raise ValueError(f"Missing exact 22-base cycle closure at line {number}")
            core = sequence[:-22]
            for oriented in (core, reverse_complement(core)):
                circular = oriented + oriented[:19]
                for pos in range(len(oriented)):
                    seed = circular[pos:pos + 20]
                    if set(seed) <= DNA:
                        seeds.add(seed)
            records += 1
    return frozenset(seeds), records


def eligible(sequence, minimum, maximum):
    return isinstance(sequence, str) and minimum <= len(sequence) <= maximum and set(sequence) <= DNA


def canonical_signature(category, repeats, spacers):
    forward = (tuple(sorted(set(repeats))), tuple(sorted(set(spacers))))
    reverse = (tuple(sorted({reverse_complement(s) for s in repeats})),
               tuple(sorted({reverse_complement(s) for s in spacers})))
    rs, ss = min(forward, reverse)
    return category, rs, ss


def initialize_worker(root, seeds):
    global _ROOT, _SEEDS
    _ROOT, _SEEDS = Path(root), seeds


def scan_shard(row):
    shard_id, relative, expected_hash, expanded_bytes, expected_counts_text = row
    if relative != f"shards/{shard_id:08d}.sqlite3.gz":
        raise ValueError(f"Unexpected shard path: {relative}")
    path = _ROOT / relative
    if path.is_symlink() or path.stat().st_size > MAX_COMPRESSED_BYTES:
        raise ValueError(f"Invalid or oversized shard: {relative}")
    compressed = path.read_bytes()
    if hashlib.sha256(compressed).hexdigest() != expected_hash:
        raise ValueError(f"Shard checksum mismatch: {relative}")
    if not 0 < expanded_bytes <= MAX_SHARD_BYTES:
        raise ValueError(f"Invalid declared shard expansion: {relative}")
    with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as handle:
        data = handle.read(MAX_SHARD_BYTES + 1)
    if len(data) != expanded_bytes:
        raise ValueError(f"Shard expansion mismatch: {relative}")
    connection = sqlite3.connect(":memory:")
    try:
        connection.deserialize(data)
        connection.execute("PRAGMA query_only=ON")
        connection.execute("PRAGMA trusted_schema=OFF")
        arrays = connection.execute("""SELECT a.array_pk, x.accession, a.category,
            a.certainty_score, a.repeat_count, a.spacer_count, a.repeat_consensus,
            a.tool, a.tool_version
            FROM arrays a JOIN sequences s USING(sequence_pk)
            JOIN assemblies x USING(assembly_pk) ORDER BY a.array_pk""").fetchall()
        stats = Counter(shards_verified=1, arrays_seen=len(arrays),
                        compressed_bytes=len(compressed), expanded_bytes=len(data))
        repeats = defaultdict(set)
        for array_pk, sequence in connection.execute("SELECT array_pk,sequence FROM repeats"):
            stats["repeat_occurrences_seen"] += 1
            if eligible(sequence, 20, 55):
                repeats[array_pk].add(sequence)
            else:
                stats["repeat_occurrences_excluded_sequence"] += 1
        selected = {}
        for array_pk, accession, category, confidence, nr, ns, consensus, tool, version in arrays:
            if tool != "CRISPRidentify" or version != "2.3.2":
                raise ValueError(f"Unexpected source tool/version in {relative}: {tool}/{version}")
            if eligible(consensus, 20, 55):
                repeats[array_pk].add(consensus)
            else:
                stats["consensuses_excluded_sequence"] += 1
            if category not in {"Bona-fide", "Possible"}:
                stats["arrays_excluded_category"] += 1
                continue
            if confidence is None or not math.isfinite(confidence) or not 0 <= confidence <= 1:
                stats["arrays_excluded_confidence"] += 1
                continue
            if not any(r[:20] in _SEEDS for r in repeats[array_pk]):
                continue
            if any(c in accession for c in "\t\n\r,"):
                raise ValueError(f"Invalid accession in {relative}")
            selected[array_pk] = (accession, category, confidence, nr, ns)
        spacers = defaultdict(set)
        for array_pk, kind, sequence in connection.execute("SELECT array_pk,kind,sequence FROM spacers"):
            stats["spacer_occurrences_seen"] += 1
            if kind != "observed":
                stats["spacer_occurrences_excluded_kind"] += 1
            elif not eligible(sequence, 20, 70):
                stats["spacer_occurrences_excluded_sequence"] += 1
            elif array_pk in selected:
                spacers[array_pk].add(sequence)
        expected_counts = json.loads(expected_counts_text)
        for source_key, stat_key in (("arrays", "arrays_seen"), ("repeats", "repeat_occurrences_seen"),
                                     ("spacers", "spacer_occurrences_seen")):
            if stats[stat_key] != expected_counts[source_key]:
                raise ValueError(f"Source count mismatch for {source_key}: {relative}")
        result = []
        for array_pk in sorted(selected):
            accession, category, confidence, nr, ns = selected[array_pk]
            if not spacers[array_pk]:
                stats["seed_selected_arrays_without_eligible_spacers"] += 1
                continue
            signature = canonical_signature(category, repeats[array_pk], spacers[array_pk])
            result.append((signature, (accession, confidence, nr or 0, ns or 0)))
            stats["seed_selected_arrays"] += 1
        return result, dict(stats)
    finally:
        connection.close()


def update_representative(target, signature, representative):
    """One best-confidence observed representative; redundancy cannot add weight."""
    previous = target.get(signature)
    rank = lambda row: (-row[1], row[0], row[2], row[3])
    if previous is None or rank(representative) < rank(previous):
        target[signature] = representative


def manifest_path(path):
    return path.with_suffix(path.suffix + ".manifest.json")


def validate_output_paths(outputs, protected):
    """Reject collisions among indices, manifests, staging files, and inputs."""
    planned = [artifact for path in outputs
               for artifact in (path, manifest_path(path), path.with_name(path.name + ".building"))]
    resolved = [path.resolve() for path in planned]
    if len(set(resolved)) != len(planned):
        raise ValueError("Output, manifest, and staging paths must be distinct")
    for path, destination in zip(planned, resolved):
        # exists() alone misses dangling symlinks, which are still user files.
        if os.path.lexists(path) or destination in protected:
            raise FileExistsError(f"Refusing to overwrite existing input/output: {path}")


def write_index(path, records, maximum_bytes):
    """Publish a complete bounded index without replacing an existing file."""
    temporary = path.with_name(path.name + ".building")
    # Open outside the try: a pre-existing staging file is not ours to remove.
    handle = temporary.open("x")
    try:
        with handle:
            handle.write(HEADER)
            if handle.tell() > maximum_bytes:
                raise ValueError(f"Index exceeds configured output bound: {path}")
            for signature in sorted(records):
                category, repeats, spacers = signature
                accession, confidence, nr, ns = records[signature]
                array_id = hashlib.sha256(json.dumps(signature, separators=(",", ":")).encode()).hexdigest()
                handle.write("\t".join((array_id, accession, category, repr(confidence),
                                         str(nr), str(ns), ",".join(repeats), ",".join(spacers))) + "\n")
                if handle.tell() > maximum_bytes:
                    raise ValueError(f"Index exceeds configured output bound: {path}")
        # Both files share a directory/filesystem. A hard link publishes the
        # closed file atomically and fails if a competing writer created path.
        os.link(temporary, path)
    finally:
        temporary.unlink()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot", type=Path, required=True)
    parser.add_argument("--cycles", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--exclusions", type=Path, help="JSON containing exclusion_accession_roots")
    parser.add_argument("--excluded-output", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--maximum-output-mib", type=int, default=1024)
    args = parser.parse_args()
    if bool(args.exclusions) != bool(args.excluded_output):
        parser.error("--exclusions and --excluded-output must be supplied together")
    if not 1 <= args.workers <= 8 or args.maximum_output_mib < 1:
        parser.error("workers must be 1–8; output limit must be positive")
    outputs = [args.output] + ([args.excluded_output] if args.excluded_output else [])
    protected = {args.cycles.resolve(), (args.snapshot / "catalog.sqlite3").resolve(),
                 (args.snapshot / "manifest.json").resolve()}
    if args.exclusions:
        protected.add(args.exclusions.resolve())
    validate_output_paths(outputs, protected)
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    seeds, cycle_records = cycle_seeds(args.cycles)
    print(json.dumps({"stage": "seeds", "cycle_records": cycle_records, "unique_20mers": len(seeds)}), flush=True)
    registry_path = args.snapshot / "catalog.sqlite3"
    registry = sqlite3.connect(registry_path.resolve().as_uri() + "?mode=ro&immutable=1", uri=True)
    try:
        registry.execute("PRAGMA query_only=ON")
        registry.execute("PRAGMA trusted_schema=OFF")
        summary = json.loads(registry.execute("SELECT value FROM export_metadata WHERE key='summary'").fetchone()[0])
        shards = registry.execute("SELECT id,relative_path,sha256,uncompressed_bytes,counts_json FROM shards ORDER BY id").fetchall()
    finally:
        registry.close()
    inventory_hash = hashlib.sha256()
    for shard_id, _, digest, _, _ in shards:
        inventory_hash.update(f"{shard_id}:{digest}\n".encode())
    if inventory_hash.hexdigest() != summary["provenance"]["shards_sha256"]:
        raise ValueError("Catalogue registry shard inventory does not match its publication")
    exclusion = json.loads(args.exclusions.read_text()) if args.exclusions else {}
    excluded_roots = set(exclusion.get("exclusion_accession_roots", []))
    if args.exclusions and not excluded_roots:
        raise ValueError("Exclusion manifest has no accession roots")
    records, excluded_records = {}, {}
    totals = Counter()
    initialize_worker(args.snapshot.resolve(), seeds)
    with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize_worker,
                             initargs=(args.snapshot.resolve(), seeds)) as executor:
        for rows, stats in executor.map(scan_shard, shards, chunksize=8):
            totals.update(stats)
            for signature, representative in rows:
                update_representative(records, signature, representative)
                if representative[0].split(".")[0] in excluded_roots:
                    totals["excluded_assembly_array_occurrences"] += 1
                elif args.excluded_output:
                    update_representative(excluded_records, signature, representative)
            if totals["shards_verified"] % 1000 == 0:
                print(json.dumps({"stage": "scan", "shards": totals["shards_verified"],
                                  "total_shards": len(shards), "retained_signatures": len(records),
                                  "seconds": round(time.monotonic() - started, 2)}), flush=True)
    for source_key, stat_key in (("arrays", "arrays_seen"), ("repeats", "repeat_occurrences_seen"),
                                 ("spacers", "spacer_occurrences_seen")):
        if totals[stat_key] != summary["counts"][source_key]:
            raise ValueError(f"Global source {source_key} count mismatch")
    maximum = args.maximum_output_mib * 1024 * 1024
    write_index(args.output, records, maximum)
    if args.excluded_output:
        write_index(args.excluded_output, excluded_records, maximum)
    manifest = {
        "schema": "zuck.array-database.projection.v1", "created_at": datetime.now(timezone.utc).isoformat(),
        "snapshot": summary, "registry": identity(registry_path), "cycles": identity(args.cycles),
        "builder": identity(Path(__file__)), "source_shards_sha256": inventory_hash.hexdigest(),
        "source_validation": "Every compressed shard hash, expanded size, and array/repeat/spacer count verified",
        "selection": {"kind": "circular-cycle-seed-projection", "seed_length": 20,
                      "both_strands": True, "cycle_records": cycle_records, "unique_seeds": len(seeds),
                      "repeat_lengths": [20, 55], "spacer_lengths": [20, 70], "alphabet": "ACGT",
                      "no_reference_annotations": True},
        "deduplication": "Joint reverse-complement array signature; category and distinct repeat/spacer sets; highest source confidence; lexical representative tie-break",
        "confidence_interpretation": "CRISPRidentify source certainty scores are heuristic, not calibrated probabilities; duplicate assemblies do not multiply support",
        "source_counts": dict(totals), "retained_signatures": len(records),
        "excluded_retained_signatures": len(excluded_records) if args.excluded_output else None,
        "exclusions": identity(args.exclusions) if args.exclusions else None,
        "exclusion_accession_roots": sorted(excluded_roots), "outputs": [identity(p) for p in outputs],
        "elapsed_seconds": time.monotonic() - started, "workers": args.workers,
        "limitations": ["This is the September region-screen snapshot, not the later October full run",
                        "Only accepted Bona-fide/Possible source arrays are available",
                        "Assembly exclusions do not remove related genomes or shared sequences",
                        "Seed projection may retain irrelevant arrays; native inference must verify full matches",
                        "Selected arrays retain eligible spacer content; DB-only spacers must not be emitted"]}
    for path in outputs:
        document = dict(manifest, output_mode="assembly-excluded" if path == args.excluded_output else "full-catalogue-projection")
        with manifest_path(path).open("x") as handle:
            handle.write(json.dumps(document, indent=2) + "\n")
    print(json.dumps({"stage": "complete", "retained_signatures": len(records),
                      "excluded_retained_signatures": len(excluded_records),
                      "seconds": manifest["elapsed_seconds"], "outputs": manifest["outputs"]}), flush=True)


if __name__ == "__main__":
    main()
