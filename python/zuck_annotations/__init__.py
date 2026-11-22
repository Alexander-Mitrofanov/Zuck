"""Source-linked annotations for observed Zuck sequences (no inference changes)."""

import csv
import gzip
import hashlib
import io
import json
import math
import re
import signal
import sqlite3
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote, urlencode

VERSION = "1.0.0"
DNA = frozenset("ACGT")
QUERY_FIELDS = ("query_id", "kind", "sequence", "database_state", "database_count",
                "protospacer_state", "protospacer_count", "skip_reason")
MAX_QUERIES = 100_000
MAX_MATCHES = 1_000_000


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def query_id(kind, sequence):
    return hashlib.sha256(f"{kind}:{sequence}".encode("ascii")).hexdigest()


def source_link(accession, start=None, end=None):
    if not accession:
        return None
    database = "assembly" if accession.startswith(("GCA_", "GCF_")) else "nuccore"
    url = f"https://www.ncbi.nlm.nih.gov/{database}/{quote(accession, safe='')}"
    if database == "nuccore" and start is not None and end is not None:
        url += "?" + urlencode({"from": min(start, end), "to": max(start, end)})
    return url


def read_queries(directory):
    """Include group consensuses and every row, independently of inference tier."""
    queries = {}
    for name, columns in (("arrays.tsv", (("repeat_consensus", "repeat"),)),
                          ("entries.tsv", (("repeat_variant", "repeat"), ("spacer", "spacer")))):
        with (directory / name).open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                for column, kind in columns:
                    sequence = row[column].upper()
                    if not sequence or len(sequence) > 20_000 or not re.fullmatch("[ACGTRYSWKMBDHVN]+", sequence):
                        raise ValueError("Invalid result sequence for annotation")
                    key = query_id(kind, sequence)
                    queries[key] = dict(query_id=key, kind=kind, sequence=sequence,
                                        database_state="not_configured", database_count=0,
                                        protospacer_state="not_applicable" if kind == "repeat" else "not_configured",
                                        protospacer_count=0, skip_reason="")
                    if len(queries) > MAX_QUERIES:
                        raise ValueError("Too many unique sequences to annotate; no completed annotation exists")
    return queries


def readonly(path):
    db = sqlite3.connect(Path(path).resolve().as_uri() + "?mode=ro&immutable=1", uri=True)
    db.row_factory = sqlite3.Row
    db.execute("PRAGMA query_only=ON")
    db.execute("PRAGMA trusted_schema=OFF")
    db.execute("PRAGMA cache_size=-4096")
    return db


def catalogue_matches(root, queries, emit, checkpoint):
    """Use CRISPRloci's exact index to nominate shards, then verify full sequences.

    The index contains 128-bit hash prefixes, so a candidate alone is never a hit.
    All nominated shards are scanned; neither source rows nor queries are capped
    silently. Repeat consensuses and individual occurrences are both included.
    """
    registry_path = root / "catalog.sqlite3"
    index_path = root / "exact-search.sqlite3"
    registry_digest = sha256(registry_path)
    registry, index = readonly(registry_path), readonly(index_path)
    try:
        summary = json.loads(registry.execute("SELECT value FROM export_metadata WHERE key='summary'").fetchone()[0])
        metadata = dict(index.execute("SELECT key,value FROM metadata"))
        expected = {"schema": "crisprloci.exact-search.v1", "term_version": "source-label-aliases-v1",
                    "status": "complete", "snapshot_id": summary["snapshot_id"],
                    "registry_sha256": registry_digest, "shards_sha256": summary["provenance"]["shards_sha256"]}
        if any(metadata.get(k) != v for k, v in expected.items()):
            raise ValueError("Catalogue index does not match its snapshot")
        nominations = defaultdict(set)
        by_sequence = defaultdict(list)
        for key, query in queries.items():
            checkpoint()
            sequence, kind = query["sequence"], query["kind"]
            if not set(sequence) <= DNA:
                query["database_state"] = "ineligible"
                continue
            query["database_state"] = "complete"
            for oriented in sorted({sequence, reverse_complement(sequence)}):
                by_sequence[(kind, oriented)].append(key)
                for term_kind in ((1, 2) if kind == "repeat" else (3,)):
                    digest = hashlib.sha256(oriented.encode()).digest()[:16]
                    for (shard_id,) in index.execute("SELECT shard_id FROM hits WHERE kind=? AND key_hash=? ORDER BY shard_id", (term_kind, digest)):
                        nominations[shard_id].add(key)
        for shard_id in sorted(nominations):
            checkpoint()
            row = registry.execute("SELECT relative_path,sha256,uncompressed_bytes FROM shards WHERE id=?", (shard_id,)).fetchone()
            if row is None or row[0] != f"shards/{shard_id:08d}.sqlite3.gz" or not 0 < row[2] <= 128 * 1024**2:
                raise ValueError("Invalid catalogue shard reference")
            path = root / row[0]
            if path.is_symlink() or path.stat().st_size > 64 * 1024**2:
                raise ValueError("Invalid catalogue shard")
            with path.open('rb') as stream:
                compressed = stream.read(64 * 1024**2 + 1)
            if len(compressed) > 64 * 1024**2 or hashlib.sha256(compressed).hexdigest() != row[1]:
                raise ValueError("Catalogue shard checksum mismatch")
            with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as stream:
                data = stream.read(row[2] + 1)
            if len(data) != row[2]:
                raise ValueError("Catalogue shard expansion mismatch")
            shard = sqlite3.connect(":memory:")
            try:
                shard.deserialize(data)
                shard.row_factory = sqlite3.Row
                shard.execute("PRAGMA query_only=ON")
                shard.execute("PRAGMA trusted_schema=OFF")
                base = """ FROM arrays a JOIN sequences s USING(sequence_pk)
                    JOIN assemblies x USING(assembly_pk) """
                fields = """a.array_pk, a.call_id, a.category, a.certainty_score, a.tool, a.tool_version,
                    s.sequence_id AS accession, x.accession AS assembly_accession"""
                scans = (
                    ("repeat", "consensus", "SELECT " + fields + ", a.array_pk AS source_id, a.repeat_consensus AS sequence, a.start, a.end" + base),
                    ("repeat", "repeat", "SELECT " + fields + ", e.repeat_pk AS source_id, e.sequence, e.start, e.end" + base + "JOIN repeats e USING(array_pk)"),
                    ("spacer", "spacer", "SELECT " + fields + ", e.spacer_pk AS source_id, e.sequence, e.start, e.end" + base + "JOIN spacers e USING(array_pk) WHERE e.kind='observed'"),
                )
                for kind, source_kind, sql in scans:
                    for n, record in enumerate(shard.execute(sql)):
                        if n % 1000 == 0:
                            checkpoint()
                        sequence = (record["sequence"] or "").upper()
                        for key in by_sequence.get((kind, sequence), ()):
                            if key not in nominations[shard_id]:
                                continue
                            match = dict(record)
                            match.update(source_kind=source_kind, snapshot_id=summary["snapshot_id"],
                                         relative_strand="+" if sequence == queries[key]["sequence"] else "-",
                                         url=source_link(record["accession"], record["start"], record["end"]),
                                         assembly_url=source_link(record["assembly_accession"]))
                            emit(key, "database", match)
            finally:
                shard.close()
        if sha256(registry_path) != registry_digest:
            raise ValueError("Catalogue changed during annotation")
        return {"state": "complete", "snapshot_id": summary["snapshot_id"],
                "registry_sha256": registry_digest, "shards_sha256": expected["shards_sha256"],
                "method": "full-length exact sequence, both strands", "candidate_shards_checked": len(nominations)}
    finally:
        registry.close()
        index.close()


def projection_matches(path, queries, emit, checkpoint):
    """Also support the native --array-database TSV, with its limited scope explicit."""
    digest = sha256(path)
    by_sequence = defaultdict(list)
    for key, query in queries.items():
        if set(query["sequence"]) <= DNA:
            query["database_state"] = "complete"
            for sequence in {query["sequence"], reverse_complement(query["sequence"])}:
                by_sequence[(query["kind"], sequence)].append(key)
        else:
            query["database_state"] = "ineligible"
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != ["array_id", "accession", "category", "confidence", "repeat_count", "spacer_count", "repeats", "spacers"]:
            raise ValueError("Invalid array database header")
        for row in reader:
            checkpoint()
            confidence = float(row["confidence"])
            if not math.isfinite(confidence) or not 0 <= confidence <= 1:
                raise ValueError("Invalid catalogue confidence")
            for kind in ("repeat", "spacer"):
                for sequence in sorted(set(row[kind + "s"].split(","))):
                    for key in by_sequence.get((kind, sequence), ()):
                        emit(key, "database", {"source_id": row["array_id"], "source_kind": "array_projection",
                            "accession": row["accession"], "category": row["category"], "certainty_score": confidence,
                            "sequence": sequence, "relative_strand": "+" if sequence == queries[key]["sequence"] else "-",
                            "url": source_link(row["accession"])})
    if sha256(path) != digest:
        raise ValueError("Array database changed during annotation")
    return {"state": "complete", "sha256": digest, "scope": "supplied array projection only",
            "method": "full-length exact sequence, both strands"}


def viral_matches(work, queries, emit, *, reference, expected_sha256, executable, blastn,
                  mismatches, threads, timeout):
    if not re.fullmatch("[a-f0-9]{64}", expected_sha256 or ""):
        raise ValueError("Protospacer search requires a pinned reference SHA-256")
    if sha256(reference / "reference.json") != expected_sha256:
        raise ValueError("Protospacer reference identity mismatch")
    spacers = {key: q for key, q in queries.items() if q["kind"] == "spacer"}
    if not spacers:
        return {"state": "complete", "outcome": "no_queries", "reference_sha256": expected_sha256}
    fasta = work / "spacers.fasta"
    fasta.write_text("".join(f">{key}\n{q['sequence']}\n" for key, q in spacers.items()))
    native = work / "native-protospacers"
    commands = [
        [executable, "viral-search", "--spacers", str(fasta), "--reference", str(reference),
         "--expected-reference-sha256", expected_sha256, "--blastn", blastn,
         "--max-mismatches", str(mismatches), "--threads", str(threads), "--timeout", str(timeout), "--output", str(native)],
        [executable, "validate-viral-result", str(native), "--reference", str(reference),
         "--expected-reference-sha256", expected_sha256],
    ]
    for command in commands:
        process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        try:
            if process.wait(timeout=timeout):
                raise ValueError('CRISPRspacer search or evidence replay failed')
        except BaseException:
            # CRISPRspacer owns a separate BLAST process group. SIGINT invokes
            # its finally cleanup; killing its Python parent first leaks BLAST.
            if process.poll() is None:
                process.send_signal(signal.SIGINT)
                try:
                    process.wait(timeout=4)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()
            raise
    report = json.loads((native / "result.json").read_text())
    inputs = json.loads((native / "inputs.json").read_text())
    if (report.get("schema") != {"name": "crisprspacer.viral-search", "version": "1.0.0"}
        or report.get("execution_completed") is not True or report.get("retrieval_completeness") != "heuristic"
        or report["reference"]["reference_sha256"] != expected_sha256
        or report["policy"]["max_substitutions"] != mismatches
        or inputs["spacer_fasta_sha256"] != sha256(fasta)
        or [(q["original_id"], q["sequence"]) for q in inputs["queries"]] != [(k, q["sequence"]) for k, q in spacers.items()]):
        raise ValueError("Protospacer result does not bind to Zuck's observed spacers and reference")
    seen = set()
    for query in report["queries"]:
        key = query["original_id"]
        if key not in spacers or key in seen or query["sha256"] != hashlib.sha256(spacers[key]["sequence"].encode()).hexdigest():
            raise ValueError("Invalid protospacer query identity")
        seen.add(key)
        spacers[key]["protospacer_state"] = "complete" if query["eligible"] else "ineligible"
        spacers[key]["skip_reason"] = ", ".join(query["skip_reasons"])
    if seen != spacers.keys():
        raise ValueError("Missing protospacer queries")
    for match in report["matches"]:
        key = match["query_id"]
        if key not in spacers or match["query_sequence"] != spacers[key]["sequence"]:
            raise ValueError("Invalid protospacer match binding")
        if (match["query_coverage"] != 1 or match["gaps"] != 0
            or not 0 <= match["substitutions"] <= mismatches
            or len(match["subject_alignment"]) != len(spacers[key]["sequence"])):
            raise ValueError("Invalid protospacer match policy")
        emit(key, "protospacer", {**match, "accession": match["subject_accession"],
             "url": source_link(match["subject_accession"], match["subject_start"], match["subject_end"])})
    return {"state": "complete", "method": "CRISPRspacer viral-search (CRISPRloci v2)",
            **{k: report[k] for k in ("tool_version", "reference", "policy", "coordinate_system",
                                      "retrieval_completeness", "counts", "limitations")},
            "blast": {"settings": report["blast"]["settings"],
                      "executable": {k: v for k, v in report["blast"]["executable"].items() if k != "path"}}}


def annotate(results, destination=None, *, catalogue=None, array_database=None,
             viral_reference=None, viral_sha256=None, crisprspacer="crisprspacer", blastn="blastn",
             mismatches=2, threads=2, timeout=300, checkpoint=lambda: None):
    results = Path(results).resolve()
    destination = Path(destination).resolve() if destination else results / "annotations"
    if destination.exists():
        raise FileExistsError("Annotation destination already exists; choose a new --output")
    if catalogue and array_database:
        raise ValueError("Choose a full catalogue or an array projection")
    if bool(viral_reference) != bool(viral_sha256):
        raise ValueError("Supply both the protospacer reference and its SHA-256")
    queries = read_queries(results)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".annotations-", dir=destination.parent) as temporary:
        work = Path(temporary) / "result"
        work.mkdir()
        counts = defaultdict(int)
        with (work / "matches.jsonl").open("w") as matches:
            def emit(key, kind, evidence):
                checkpoint()
                counts[kind] += 1
                if sum(counts.values()) > MAX_MATCHES:
                    raise ValueError("Annotation match budget exceeded; no completed annotation exists")
                queries[key][kind + "_count"] += 1
                matches.write(json.dumps({"query_id": key, "kind": kind, "evidence": evidence}, allow_nan=False) + "\n")
            database = {"state": "not_configured"}
            if catalogue:
                database = catalogue_matches(Path(catalogue), queries, emit, checkpoint)
            elif array_database:
                database = projection_matches(Path(array_database), queries, emit, checkpoint)
            protospacer = {"state": "not_configured"}
            if viral_reference:
                protospacer = viral_matches(work, queries, emit, reference=Path(viral_reference).resolve(),
                    expected_sha256=viral_sha256, executable=str(crisprspacer), blastn=str(blastn),
                    mismatches=mismatches, threads=threads, timeout=timeout)
        with (work / "queries.tsv").open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=QUERY_FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(queries.values())
        info = {"schema_version": 1, "annotator_version": VERSION, "database": database,
                "protospacer": protospacer, "counts": dict(counts), "query_count": len(queries),
                "input_sha256": {name: sha256(results / name) for name in ("arrays.tsv", "entries.tsv")},
                "files": {name: sha256(work / name) for name in ("queries.tsv", "matches.jsonl")}}
        (work / "manifest.json").write_text(json.dumps(info, indent=2, allow_nan=False) + "\n")
        checkpoint()
        work.rename(destination)
    return info
