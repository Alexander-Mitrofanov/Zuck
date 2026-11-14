import test from 'node:test';
import assert from 'node:assert/strict';
import { evidenceUrl, evidenceLabel } from '../matches.js';

test('source links accept NCBI records and reject active/foreign URLs', () => {
  assert.equal(evidenceUrl('https://www.ncbi.nlm.nih.gov/nuccore/NC_1.1?from=1&to=30'), 'https://www.ncbi.nlm.nih.gov/nuccore/NC_1.1?from=1&to=30');
  for (const value of ['javascript:alert(1)', 'https://evil.test/nuccore/x', 'https://www.ncbi.nlm.nih.gov.evil.test/nuccore/x', 'http://www.ncbi.nlm.nih.gov/nuccore/x', 'https://name:pass@www.ncbi.nlm.nih.gov/nuccore/x']) assert.equal(evidenceUrl(value), null);
});

test('protospacer labels retain strand, coordinates and alignment evidence', () => {
  assert.equal(evidenceLabel({accession:'NC_1.1',subject_start:12,subject_end:42,relative_strand:'-',substitutions:1,query_coverage:1},'protospacer'), 'NC_1.1:12–42 (-) · 1 substitution · 100% coverage');
});
