"""Annotation controls with synthetic catalogue sources and native BLAST."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'python'))
from zuck_annotations import annotate, query_id, reverse_complement, sha256

REPEAT = 'GTTCACTGCCGTACAGGCAGCTTAGAAA'
SPACER = 'ACGTCGATGCTAGCTACCGATCGTACGATCGT'
OTHER = 'TTCAGGCTAACCTGAGTCGATGCCAT'


def native_results(root, repeat=REPEAT, spacer=SPACER):
    root.mkdir(parents=True, exist_ok=True)
    (root/'arrays.tsv').write_text('array_id\tkind\trepeat_consensus\tentry_count\nArray_1\tprimary\t'+repeat+'\t1\n')
    (root/'entries.tsv').write_text('array_id\tentry_index\trepeat_variant\tspacer\nArray_1\t1\t'+repeat+'\t'+spacer+'\n')
    (root/'run_manifest.json').write_text(json.dumps({'schema_version':1,'tool':'Zuck','counts':{'primary_arrays':1,'primary_entries':1}}))


def catalogue(root):
    (root/'shards').mkdir(parents=True)
    db=sqlite3.connect(':memory:')
    db.executescript('''
      CREATE TABLE assemblies(assembly_pk INTEGER PRIMARY KEY,accession TEXT);
      CREATE TABLE sequences(sequence_pk INTEGER PRIMARY KEY,assembly_pk INTEGER,sequence_id TEXT);
      CREATE TABLE arrays(array_pk INTEGER PRIMARY KEY,sequence_pk INTEGER,call_id TEXT,category TEXT,certainty_score REAL,tool TEXT,tool_version TEXT,repeat_consensus TEXT,start INTEGER,end INTEGER);
      CREATE TABLE repeats(repeat_pk INTEGER PRIMARY KEY,array_pk INTEGER,sequence TEXT,start INTEGER,end INTEGER);
      CREATE TABLE spacers(spacer_pk INTEGER PRIMARY KEY,array_pk INTEGER,sequence TEXT,start INTEGER,end INTEGER,kind TEXT);
      INSERT INTO assemblies VALUES(1,'GCF_000000001.1');
      INSERT INTO sequences VALUES(1,1,'NC_000001.1');
    ''')
    # Spacer independently matches an array with a different repeat; a repeat
    # family/group threshold must not suppress its sequence evidence.
    db.execute('INSERT INTO arrays VALUES(1,1,?,?,?,?,?,?,100,200)',('call1','Possible',0.6,'CRISPRidentify','2.3.2',REPEAT))
    db.execute('INSERT INTO arrays VALUES(2,1,?,?,?,?,?,?,300,400)',('call2','Bona-fide',0.99,'CRISPRidentify','2.3.2',OTHER))
    db.execute('INSERT INTO repeats VALUES(1,1,?,100,126)',(reverse_complement(REPEAT),))
    db.execute('INSERT INTO spacers VALUES(1,2,?,330,361,\'observed\')',(reverse_complement(SPACER),))
    db.execute('INSERT INTO spacers VALUES(2,2,?,370,401,\'inferred\')',(SPACER,))
    db.commit()
    data=db.serialize();db.close()
    compressed=gzip.compress(data)
    path=root/'shards/00000001.sqlite3.gz';path.write_bytes(compressed)
    shard_hash=sha256(path)
    inventory=hashlib.sha256(f'1:{shard_hash}\n'.encode()).hexdigest()
    summary={'snapshot_id':'synthetic-test','provenance':{'shards_sha256':inventory}}
    db=sqlite3.connect(root/'catalog.sqlite3')
    db.executescript('CREATE TABLE export_metadata(key TEXT,value TEXT); CREATE TABLE shards(id INTEGER PRIMARY KEY,relative_path TEXT,sha256 TEXT,uncompressed_bytes INTEGER);')
    db.execute('INSERT INTO export_metadata VALUES(?,?)',('summary',json.dumps(summary)))
    db.execute('INSERT INTO shards VALUES(1,?,?,?)',('shards/00000001.sqlite3.gz',shard_hash,len(data)))
    db.commit();db.close()
    db=sqlite3.connect(root/'exact-search.sqlite3')
    db.executescript('CREATE TABLE metadata(key TEXT,value TEXT); CREATE TABLE hits(kind INTEGER,key_hash BLOB,shard_id INTEGER,PRIMARY KEY(kind,key_hash,shard_id));')
    meta={'schema':'crisprloci.exact-search.v1','term_version':'source-label-aliases-v1','status':'complete',
          'snapshot_id':'synthetic-test','registry_sha256':sha256(root/'catalog.sqlite3'),'shards_sha256':inventory}
    db.executemany('INSERT INTO metadata VALUES(?,?)',meta.items())
    for kind,sequence in ((1,REPEAT),(2,reverse_complement(REPEAT)),(3,reverse_complement(SPACER)),(3,SPACER)):
        db.execute('INSERT INTO hits VALUES(?,?,1)',(kind,hashlib.sha256(sequence.encode()).digest()[:16]))
    db.commit();db.close()


def projection(path, repeat=REPEAT, spacer=SPACER, copies=1):
    path.write_text('array_id\taccession\tcategory\tconfidence\trepeat_count\tspacer_count\trepeats\tspacers\n'+
        ''.join(f'source{i}\tGCF_{i+1:09}.1\tPossible\t0.6\t3\t2\t{repeat}\t{spacer}\n' for i in range(copies)))


def viral_reference(root, executable, spacer=SPACER):
    source=root/'viruses.fa'
    # Forward exact, reverse exact, internal mismatch, negative target.
    mutant=spacer[:15]+('A' if spacer[15]!='A' else 'C')+spacer[16:]
    source.write_text(f'>NC_000001.1 Synthetic exact\nTTTTT{spacer}AAAAA\n'
                      f'>NC_000002.1 Synthetic reverse\nGGGGG{reverse_complement(spacer)}CCCCC\n'
                      f'>NC_000003.1 Synthetic substitution\nTTTTT{mutant}AAAAA\n')
    reference=root/'viral-reference'
    subprocess.run([executable,'viral-index','--sources',str(source),'--name','synthetic',
        '--source-release','test','--source-url','https://example.org/synthetic',
        '--output',str(reference)],check=True,stdout=subprocess.DEVNULL)
    return reference


class AnnotationTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name)
        self.result=self.root/'output'
        native_results(self.result)

    def test_every_sequence_both_strands_and_source_links(self):
        catalogue(self.root/'catalogue')
        info=annotate(self.result,catalogue=self.root/'catalogue')
        hits=[json.loads(line) for line in (self.result/'annotations/matches.jsonl').read_text().splitlines()]
        self.assertEqual(info['counts']['database'],3)
        self.assertEqual({h['evidence']['source_kind'] for h in hits},{'consensus','repeat','spacer'})
        spacer=next(h for h in hits if h['query_id']==query_id('spacer',SPACER))
        self.assertEqual(spacer['evidence']['relative_strand'],'-')
        self.assertEqual(spacer['evidence']['url'],'https://www.ncbi.nlm.nih.gov/nuccore/NC_000001.1?from=330&to=361')
        self.assertEqual(spacer['evidence']['array_pk'],2)
        self.assertEqual(info['database']['candidate_shards_checked'],1)

    def test_corrupt_shard_fails_atomically(self):
        catalogue(self.root/'catalogue')
        (self.root/'catalogue/shards/00000001.sqlite3.gz').write_bytes(b'corrupt')
        with self.assertRaisesRegex(ValueError,'checksum'):
            annotate(self.result,catalogue=self.root/'catalogue')
        self.assertFalse((self.result/'annotations').exists())

    def test_index_identity_and_false_candidates(self):
        catalogue(self.root/'catalogue')
        db=sqlite3.connect(self.root/'catalogue/exact-search.sqlite3')
        db.execute('INSERT INTO hits VALUES(?,?,1)',(3,hashlib.sha256(('A'*30).encode()).digest()[:16]))
        db.commit();db.close()
        native_results(self.result,spacer='A'*30)
        annotate(self.result,catalogue=self.root/'catalogue')
        rows=list(csv.DictReader((self.result/'annotations/queries.tsv').read_text().splitlines(),delimiter='\t'))
        self.assertEqual(next(r for r in rows if r['kind']=='spacer')['database_count'],'0')

    def test_projection_all_sources_and_unknown_states(self):
        projection(self.root/'arrays.tsv',copies=8)
        info=annotate(self.result,array_database=self.root/'arrays.tsv')
        self.assertEqual(info['counts']['database'],16)
        rows=list(csv.DictReader((self.result/'annotations/queries.tsv').read_text().splitlines(),delimiter='\t'))
        self.assertEqual(next(r for r in rows if r['kind']=='spacer')['protospacer_state'],'not_configured')
        with self.assertRaises(FileExistsError):
            annotate(self.result)

    def test_unconfigured_is_not_a_negative_search_and_empty_results(self):
        info=annotate(self.result)
        self.assertEqual(info['database']['state'],'not_configured')
        self.assertEqual(info['protospacer']['state'],'not_configured')
        empty=self.root/'empty';empty.mkdir()
        (empty/'arrays.tsv').write_text('array_id\trepeat_consensus\tentry_count\n')
        (empty/'entries.tsv').write_text('array_id\tentry_index\trepeat_variant\tspacer\n')
        self.assertEqual(annotate(empty)['query_count'],0)

    @unittest.skipUnless(os.getenv('ZUCK_CRISPRSPACER_TEST') and shutil.which('blastn'),'set ZUCK_CRISPRSPACER_TEST to test native BLAST')
    def test_native_protospacers_replay_both_strands_and_mismatch(self):
        executable=os.environ['ZUCK_CRISPRSPACER_TEST']
        reference=viral_reference(self.root,executable)
        info=annotate(self.result,viral_reference=reference,viral_sha256=sha256(reference/'reference.json'),crisprspacer=executable)
        hits=[json.loads(line)['evidence'] for line in (self.result/'annotations/matches.jsonl').read_text().splitlines()]
        self.assertEqual({h['accession'] for h in hits},{'NC_000001.1','NC_000002.1','NC_000003.1'})
        self.assertEqual({h['relative_strand'] for h in hits},{'+','-'})
        self.assertEqual({h['substitutions'] for h in hits},{0,1})
        self.assertEqual(info['protospacer']['retrieval_completeness'],'heuristic')
        self.assertTrue(all(h['query_coverage']==1 and h['gaps']==0 for h in hits))
        # A negative and an ineligible query remain distinct; clipped matches
        # from the positive query must not create a partial accepted target.
        negative=self.root/'negative'
        native_results(negative,spacer='A'*32)
        annotate(negative,viral_reference=reference,viral_sha256=sha256(reference/'reference.json'),crisprspacer=executable)
        rows=list(csv.DictReader((negative/'annotations/queries.tsv').read_text().splitlines(),delimiter='\t'))
        row=next(r for r in rows if r['kind']=='spacer')
        self.assertEqual((row['protospacer_state'],row['protospacer_count']),('complete','0'))
        ambiguous=self.root/'ambiguous'
        native_results(ambiguous,spacer='N'*32)
        annotate(ambiguous,viral_reference=reference,viral_sha256=sha256(reference/'reference.json'),crisprspacer=executable)
        rows=list(csv.DictReader((ambiguous/'annotations/queries.tsv').read_text().splitlines(),delimiter='\t'))
        row=next(r for r in rows if r['kind']=='spacer')
        self.assertEqual(row['protospacer_state'],'ineligible')
        self.assertIn('ambiguous_query',row['skip_reason'])


if __name__=='__main__':
    unittest.main()
