#!/usr/bin/env python3
"""Independent negative controls for the catalogue projection preparer."""
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/build_array_database.py"
SPEC = importlib.util.spec_from_file_location("build_array_database", SCRIPT)
DB = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(DB)
REPEAT = "ACGTTGCAACGATCGATTCGAGCTA"
SPACER1 = "GGCATACCTTAGGCTATGACCGATGCA"
SPACER2 = "TCAGGCTAGACGTACCTTGGATCGAA"


def make_snapshot(path):
    (path / "shards").mkdir(parents=True)
    connection = sqlite3.connect(":memory:")
    connection.executescript("""
      CREATE TABLE assemblies(assembly_pk INTEGER,accession TEXT);
      CREATE TABLE sequences(sequence_pk INTEGER,assembly_pk INTEGER);
      CREATE TABLE arrays(array_pk INTEGER,sequence_pk INTEGER,category TEXT,certainty_score REAL,
        repeat_count INTEGER,spacer_count INTEGER,repeat_consensus TEXT,tool TEXT,tool_version TEXT);
      CREATE TABLE repeats(array_pk INTEGER,sequence TEXT);
      CREATE TABLE spacers(array_pk INTEGER,kind TEXT,sequence TEXT);
    """)
    for pk, accession, score in ((1, "GCA_000000001.1", .95), (2, "GCA_000000002.1", .82)):
        connection.execute("INSERT INTO assemblies VALUES (?,?)", (pk, accession))
        connection.execute("INSERT INTO sequences VALUES (?,?)", (pk, pk))
        connection.execute("INSERT INTO arrays VALUES (?,?,?,?,?,?,?,?,?)",
                           (pk, pk, "Bona-fide", score, 3, 2, REPEAT, "CRISPRidentify", "2.3.2"))
        connection.executemany("INSERT INTO repeats VALUES (?,?)", [(pk, REPEAT)] * 3)
        connection.executemany("INSERT INTO spacers VALUES (?,?,?)",
                               [(pk, "observed", SPACER1), (pk, "observed", SPACER2)])
    raw = connection.serialize()
    compressed = gzip.compress(raw, mtime=0)
    digest = hashlib.sha256(compressed).hexdigest()
    (path / "shards/00000001.sqlite3.gz").write_bytes(compressed)
    counts = {"arrays": 2, "repeats": 6, "spacers": 4}
    inventory_hash = hashlib.sha256(f"1:{digest}\n".encode()).hexdigest()
    summary = {"snapshot_id": "synthetic", "counts": counts,
               "provenance": {"shards_sha256": inventory_hash}}
    registry = sqlite3.connect(path / "catalog.sqlite3")
    registry.executescript("CREATE TABLE export_metadata(key TEXT,value TEXT); CREATE TABLE shards(id INTEGER,relative_path TEXT,sha256 TEXT,uncompressed_bytes INTEGER,counts_json TEXT);")
    registry.execute("INSERT INTO export_metadata VALUES ('summary',?)", (json.dumps(summary),))
    registry.execute("INSERT INTO shards VALUES (?,?,?,?,?)",
                     (1, "shards/00000001.sqlite3.gz", digest, len(raw), json.dumps(counts)))
    registry.commit()
    registry.close()
    return (1, "shards/00000001.sqlite3.gz", digest, len(raw), json.dumps(counts))


class BuilderTests(unittest.TestCase):
    def test_circular_and_reverse_seed_invariance(self):
        with tempfile.TemporaryDirectory() as directory:
            core = REPEAT + SPACER1
            rotated = core[13:] + core[:13]
            paths = []
            for index, sequence in enumerate((core, rotated, DB.reverse_complement(core))):
                path = Path(directory) / str(index)
                path.write_text(sequence + sequence[:22] + "\n")
                paths.append(path)
            self.assertEqual(DB.cycle_seeds(paths[0]), DB.cycle_seeds(paths[1]))
            self.assertEqual(DB.cycle_seeds(paths[0]), DB.cycle_seeds(paths[2]))

    def test_missing_closure_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cycles"
            path.write_text(REPEAT + SPACER1 + "A" * 22 + "\n")
            with self.assertRaisesRegex(ValueError, "closure"):
                DB.cycle_seeds(path)

    def test_legacy_kmer_gluing(self):
        self.assertEqual(DB.glue_cycle("ACGT CGTT GTTA"), "ACGTTA")

    def test_signature_joint_orientation_and_duplicate_invariance(self):
        signature = DB.canonical_signature("Bona-fide", [REPEAT, REPEAT], [SPACER1, SPACER2])
        reverse = DB.canonical_signature("Bona-fide", [DB.reverse_complement(REPEAT)],
                                         [DB.reverse_complement(SPACER2), DB.reverse_complement(SPACER1)])
        self.assertEqual(signature, reverse)

    def test_duplicate_confidence_max_not_sum(self):
        records = {}
        for value in [("GCA_2.1", .9, 3, 2), ("GCA_1.1", .8, 3, 2), ("GCA_0.1", .9, 3, 2)]:
            DB.update_representative(records, ("signature",), value)
        self.assertEqual(records[("signature",)], ("GCA_0.1", .9, 3, 2))

    def test_ambiguous_or_out_of_range_sequence_rejected(self):
        self.assertFalse(DB.eligible("N" * 25, 20, 55))
        self.assertFalse(DB.eligible("A" * 19, 20, 55))
        self.assertFalse(DB.eligible("A" * 56, 20, 55))

    def test_seed_miss_cannot_select_array(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            row = make_snapshot(root)
            DB.initialize_worker(root, frozenset({"T" * 20}))
            rows, stats = DB.scan_shard(row)
            self.assertEqual(rows, [])
            self.assertEqual(stats["arrays_seen"], 2)

    def test_corrupt_source_hash_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            row = make_snapshot(root)
            DB.initialize_worker(root, frozenset({REPEAT[:20]}))
            with (root / row[1]).open("ab") as handle:
                handle.write(b"corrupt")
            with self.assertRaisesRegex(ValueError, "checksum"):
                DB.scan_shard(row)

    def test_wrong_source_count_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            row = list(make_snapshot(root))
            row[-1] = json.dumps({"arrays": 99, "repeats": 6, "spacers": 4})
            DB.initialize_worker(root, frozenset({REPEAT[:20]}))
            with self.assertRaisesRegex(ValueError, "count mismatch"):
                DB.scan_shard(row)

    def test_output_artifact_collisions_fail_before_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "full.tsv"
            for collision in (DB.manifest_path(output), output.with_name("full.tsv.building")):
                with self.assertRaisesRegex(ValueError, "distinct"):
                    DB.validate_output_paths([output, collision], set())
            self.assertEqual(list(Path(directory).iterdir()), [])

    def test_dangling_output_symlink_is_not_replaced(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "full.tsv"
            output.symlink_to(Path(directory) / "absent")
            with self.assertRaises(FileExistsError):
                DB.validate_output_paths([output], set())
            self.assertTrue(output.is_symlink())

    def test_failed_index_write_removes_owned_staging_file(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "full.tsv"
            records = {DB.canonical_signature("Bona-fide", [REPEAT], [SPACER1]):
                       ("GCA_1.1", 0.9, 2, 1)}
            with self.assertRaisesRegex(ValueError, "output bound"):
                DB.write_index(output, records, len(DB.HEADER))
            self.assertEqual(list(Path(directory).iterdir()), [])
            DB.write_index(output, records, 4096)
            before = output.read_bytes()
            with self.assertRaises(FileExistsError):
                DB.write_index(output, {}, 4096)
            self.assertEqual(output.read_bytes(), before)
            self.assertFalse(output.with_name("full.tsv.building").exists())

    def test_preexisting_staging_file_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "full.tsv"
            staging = output.with_name("full.tsv.building")
            staging.write_text("other writer")
            with self.assertRaises(FileExistsError):
                DB.write_index(output, {}, 4096)
            self.assertEqual(staging.read_text(), "other writer")

    def test_end_to_end_exclusion_recomputes_best_confidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            snapshot = root / "snapshot"
            make_snapshot(snapshot)
            cycles = root / "cycles"
            core = REPEAT + SPACER1
            cycles.write_text(core + core[:22] + "\n")
            exclusion = root / "exclude.json"
            exclusion.write_text(json.dumps({"exclusion_accession_roots": ["GCA_000000001"]}))
            command = [sys.executable, str(SCRIPT), "--snapshot", str(snapshot), "--cycles", str(cycles),
                       "--output", str(root / "full.tsv"), "--exclusions", str(exclusion),
                       "--excluded-output", str(root / "excluded.tsv"), "--workers", "1"]
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            full = (root / "full.tsv").read_text().splitlines()
            excluded = (root / "excluded.tsv").read_text().splitlines()
            self.assertEqual(len(full), 2)
            self.assertEqual(len(excluded), 2)
            self.assertEqual(full[1].split("\t")[3], "0.95")
            self.assertEqual(excluded[1].split("\t")[3], "0.82")
            self.assertEqual(full[1].split("\t")[0], excluded[1].split("\t")[0])
            self.assertNotEqual(subprocess.run(command, capture_output=True).returncode, 0)


if __name__ == "__main__":
    unittest.main()
