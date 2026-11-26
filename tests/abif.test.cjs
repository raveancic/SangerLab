'use strict';
// Scientific parser regression tests. All fixtures below are synthetic.
// Run: node --test tests/abif.test.cjs
const test = require('node:test');
const assert = require('node:assert/strict');
const {parseAb1, reverseComplementTrace} = require('../src/abif.js');

function bytes(values, width = 1, signed = false) {
  const b = Buffer.alloc(values.length * width);
  values.forEach((x, i) => {
    if (width === 1) b.writeUInt8(x, i);
    else if (width === 2 && signed) b.writeInt16BE(x, i * 2);
    else if (width === 2) b.writeUInt16BE(x, i * 2);
    else if (width === 4 && signed) b.writeInt32BE(x, i * 4);
    else b.writeUInt32BE(x, i * 4);
  });
  return b;
}

function tag(name, number, type, size, data, count = data.length / size) {
  return {name, number, type, size, data, count};
}

const BASES = 'ACGTRN';
const QUALITIES = [5, 10, 20, 30, 40, 50];
const PEAKS = [0, 2, 4, 6, 8, 10];
const ORDER = 'GATC';
const CHANNELS = Object.fromEntries([...ORDER].map((b, i) => [b,
  [-32768 + i, -5 - i, i, 1 + i, 513 + i, 100 + i, 20 + i,
    300 + i, 1500 + i, 20 + i, 10 + i, 40 + i, 100 + i, 32767 - i]
]));

function defaultTags(generation = 2) {
  return [
    tag('FWO_', 1, 2, 1, Buffer.from(ORDER)),
    ...[...ORDER].map((b, i) => tag('DATA', 9 + i, 4, 2, bytes(CHANNELS[b], 2, true))),
    tag('PBAS', generation, 2, 1, Buffer.from(BASES)),
    tag('PCON', generation, 2, 1, bytes(QUALITIES)),
    tag('PLOC', generation, 4, 2, bytes(PEAKS, 2, true))
  ];
}

function fixture(tags = defaultTags(), version = 101) {
  const directory = 128;
  const dataStart = directory + tags.length * 28;
  const dataBytes = tags.reduce((n, t) => n + (t.data.length > 4 ? t.data.length : 0), 0);
  const buf = Buffer.alloc(dataStart + dataBytes);
  buf.write('ABIF', 0, 'ascii');
  buf.writeUInt16BE(version, 4);
  buf.write('tdir', 6, 'ascii');
  buf.writeUInt32BE(1, 10);
  buf.writeUInt16BE(1023, 14);
  buf.writeUInt16BE(28, 16);
  buf.writeUInt32BE(tags.length, 18);
  buf.writeUInt32BE(tags.length * 28, 22);
  buf.writeUInt32BE(directory, 26);
  let at = dataStart;
  tags.forEach((t, i) => {
    const d = directory + i * 28;
    buf.write(t.name, d, 'ascii');
    buf.writeUInt32BE(t.number, d + 4);
    buf.writeUInt16BE(t.type, d + 8);
    buf.writeUInt16BE(t.size, d + 10);
    buf.writeUInt32BE(t.count, d + 12);
    buf.writeUInt32BE(t.data.length, d + 16);
    if (t.data.length <= 4) t.data.copy(buf, d + 20);
    else {
      buf.writeUInt32BE(at, d + 20);
      t.data.copy(buf, at);
      at += t.data.length;
    }
  });
  return buf;
}

function replacing(tags, name, number, replacement) {
  return tags.map(t => t.name === name && t.number === number ? replacement : t);
}

test('parses analyzed trace, original calls, quality and signed big-endian samples', () => {
  const r = parseAb1(fixture(), 'synthetic.ab1');
  assert.equal(r.sequence, BASES);
  assert.deepEqual(r.qualities, QUALITIES);
  assert.deepEqual(r.peakPositions, PEAKS);
  assert.equal(r.channelOrder, ORDER);
  assert.equal(r.sampleCount, 14);
  for (const b of ORDER) assert.deepEqual(Array.from(r.channels[b]), CHANNELS[b]);
  assert.deepEqual(r.tags, {baseCalls: 'PBAS2', peaks: 'PLOC2', qualities: 'PCON2'});
  assert.deepEqual(r.warnings, []);
  assert.equal(r.abiVersion, 101);
});

test('also accepts ArrayBuffer and an offset byte view without adjacent bytes', () => {
  const b = fixture();
  const large = Buffer.concat([Buffer.from('prefix'), b, Buffer.from('suffix')]);
  assert.equal(parseAb1(large.subarray(6, 6 + b.length)).sequence, BASES);
  assert.equal(parseAb1(b.buffer.slice(b.byteOffset, b.byteOffset + b.length)).sequence, BASES);
});

test('inline fields are data bytes, not data offsets (one-base calls/quality/peak)', () => {
  let tags = defaultTags();
  tags = replacing(tags, 'PBAS', 2, tag('PBAS', 2, 2, 1, Buffer.from('T')));
  tags = replacing(tags, 'PCON', 2, tag('PCON', 2, 2, 1, bytes([37])));
  tags = replacing(tags, 'PLOC', 2, tag('PLOC', 2, 4, 2, bytes([5], 2, true)));
  const r = parseAb1(fixture(tags));
  assert.equal(r.sequence, 'T');
  assert.deepEqual(r.qualities, [37]);
  assert.deepEqual(r.peakPositions, [5]);
});

test('absence of quality is explicit; generation-1 quality is not substituted', () => {
  const tags = defaultTags().filter(t => t.name !== 'PCON');
  tags.push(tag('PCON', 1, 2, 1, bytes([99, 99, 99, 99, 99, 99])));
  const r = parseAb1(fixture(tags));
  assert.equal(r.qualities, null);
  assert.equal(r.tags.qualities, null);
  assert.match(r.warnings.join(' '), /PCON2 assente/);
});

test('generation-1 fallback is coherent and reported', () => {
  const r = parseAb1(fixture(defaultTags(1)));
  assert.equal(r.tags.baseCalls, 'PBAS1');
  assert.equal(r.tags.peaks, 'PLOC1');
  assert.equal(r.tags.qualities, 'PCON1');
  assert.equal(r.warnings.length, 1);
});

test('complete generation 2 takes precedence over edited generation 1', () => {
  const tags = defaultTags();
  tags.push(tag('PBAS', 1, 2, 1, Buffer.from('TTTTTT')));
  tags.push(tag('PCON', 1, 2, 1, bytes([1, 1, 1, 1, 1, 1])));
  tags.push(tag('PLOC', 1, 4, 2, bytes([1, 3, 5, 7, 9, 11], 2, true)));
  const r = parseAb1(fixture(tags));
  assert.equal(r.sequence, BASES);
  assert.deepEqual(r.qualities, QUALITIES);
  assert.deepEqual(r.peakPositions, PEAKS);
});

test('partial generation 2 does not silently fall back to generation 1', () => {
  const tags = defaultTags(1);
  tags.push(tag('PBAS', 2, 2, 1, Buffer.from('TTTTTT')));
  assert.throws(() => parseAb1(fixture(tags)), /Coppia coerente/);
});

test('reverse complement transforms all coordinates, channels and quality together', () => {
  const r = parseAb1(fixture());
  const rc = reverseComplementTrace(r);
  assert.equal(rc.sequence, 'NYACGT');
  assert.deepEqual(rc.qualities, [...QUALITIES].reverse());
  assert.deepEqual(rc.peakPositions, [3, 5, 7, 9, 11, 13]);
  assert.deepEqual(rc.sourceIndices, [5, 4, 3, 2, 1, 0]);
  assert.deepEqual(Array.from(rc.channels.A), [...CHANNELS.T].reverse());
  assert.deepEqual(Array.from(rc.channels.C), [...CHANNELS.G].reverse());
  assert.deepEqual(reverseComplementTrace(rc), r);
  assert.equal(r.sequence, BASES, 'source record must remain unchanged');
});

test('reverse complement preserves missing quality', () => {
  const r = parseAb1(fixture(defaultTags().filter(t => t.name !== 'PCON')));
  assert.equal(reverseComplementTrace(r).qualities, null);
});

test('invalid signature, unsupported version and materially truncated file are rejected', () => {
  const b = fixture();
  for (const length of [0, 10, 33, 100, 129, b.length - 1]) {
    assert.throws(() => parseAb1(b.subarray(0, length)), undefined, `truncated at ${length}`);
  }
  const badHeader = Buffer.from(b); badHeader.write('NOPE', 0);
  assert.throws(() => parseAb1(badHeader), /Firma ABIF/);
  assert.throws(() => parseAb1(fixture(defaultTags(), 999)), /Versione ABIF/);
});

test('missing or duplicate channel order, and invalid nucleotide order are rejected', () => {
  assert.throws(() => parseAb1(fixture(defaultTags().filter(t => t.name !== 'FWO_'))), /FWO_1 assente/);
  for (const order of ['AACG', 'XYZW', 'GAT', 'GATCN']) {
    const tags = replacing(defaultTags(), 'FWO_', 1, tag('FWO_', 1, 2, 1, Buffer.from(order)));
    assert.throws(() => parseAb1(fixture(tags)), /Ordine/);
  }
  const tags = defaultTags(); tags.push(tags[0]);
  assert.throws(() => parseAb1(fixture(tags)), /duplicato/);
});

test('out-of-range, negative and descending peaks are rejected', () => {
  for (const peaks of [[0, 2, 4, 6, 8, 14], [-1, 2, 4, 6, 8, 10], [0, 2, 4, 3, 8, 10]]) {
    const tags = replacing(defaultTags(), 'PLOC', 2, tag('PLOC', 2, 4, 2, bytes(peaks, 2, true)));
    assert.throws(() => parseAb1(fixture(tags)), /Coordinate/);
  }
});

test('inconsistent base, quality and trace lengths are rejected', () => {
  const shortQ = replacing(defaultTags(), 'PCON', 2, tag('PCON', 2, 2, 1, bytes([20])));
  assert.throws(() => parseAb1(fixture(shortQ)), /Numero di qualità/);
  const shortP = replacing(defaultTags(), 'PLOC', 2, tag('PLOC', 2, 4, 2, bytes([1], 2, true)));
  assert.throws(() => parseAb1(fixture(shortP)), /Numero di basi/);
  const shortT = replacing(defaultTags(), 'DATA', 9, tag('DATA', 9, 4, 2, bytes([1, 2], 2, true)));
  assert.throws(() => parseAb1(fixture(shortT)), /lunghezze diverse/);
});

test('invalid sequence characters and missing analyzed traces are rejected', () => {
  const invalid = replacing(defaultTags(), 'PBAS', 2, tag('PBAS', 2, 2, 1, Buffer.from('ACGT-X')));
  assert.throws(() => parseAb1(fixture(invalid)), /non IUPAC/);
  assert.throws(() => parseAb1(fixture(defaultTags().filter(t => !(t.name === 'DATA' && t.number === 9)))), /DATA9 assente/);
});

test('corrupted file offsets and impossible directory counts are rejected', () => {
  const b = fixture();
  const corruptRoot = Buffer.from(b); corruptRoot.writeUInt32BE(0xffffffff, 18);
  assert.throws(() => parseAb1(corruptRoot), /Numero di elementi/);
  const corruptData = Buffer.from(b);
  // Second directory entry is DATA9 (non-inline).
  corruptData.writeUInt32BE(b.length + 100, 128 + 28 + 20);
  assert.throws(() => parseAb1(corruptData), /incompleto o corrotto/);
});

test('structural ABIF300 profile is accepted (synthetic; no physical instrument claim)', () => {
  assert.equal(parseAb1(fixture(defaultTags(), 300)).sequence, BASES);
});
