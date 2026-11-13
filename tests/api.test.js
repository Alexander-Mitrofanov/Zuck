import test from 'node:test';
import assert from 'node:assert/strict';
import {parseRecovery, validateFiles, makeApi} from '../api.js';
test('recovery requires an exact job and capability fragment', () => {
 const valid = '#job='+'a'.repeat(32)+'.'+'B'.repeat(43);
 assert.deepEqual(parseRecovery(valid), {id:'a'.repeat(32), token:'B'.repeat(43)});
 for(const bad of [valid+'&extra=x',valid+'/',valid.replace('#','?'),'#job=../test',valid.replace('job=','token=')]) assert.equal(parseRecovery(bad),null);
});
test('uploads enforce mode, nonempty files, and combined size', () => {
 assert.doesNotThrow(()=>validateFiles([{size:2},{size:3}],'reads',5));
 for(const args of [[[],'reads',5],[[{size:0}],'reads',5],[[{size:3},{size:3}],'reads',5],[[{size:1},{size:1}],'cycles',5]]) assert.throws(()=>validateFiles(...args));
});
test('only HTTPS production origins and explicit localhost development are allowed',()=>{
 assert.doesNotThrow(()=>makeApi('https://example.org/zuck'));
 assert.doesNotThrow(()=>makeApi('http://127.0.0.1:8001'));
 for(const url of ['http://example.org','https://name:secret@example.org','https://example.org?token=x','javascript:alert(1)']) assert.throws(()=>makeApi(url));
});
