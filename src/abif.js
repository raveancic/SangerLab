/*
 * ABIF parser adapted and hardened from Sanger Viewer:
 * https://github.com/animeshkundu/sanger-viewer/blob/
 * 3ceda33e18dde1700f6a53414f58ad4b1d8989be/src/parsers/abif.ts
 *
 * MIT License
 * Copyright (c) 2026 Animesh
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * SangerLab changes: bounds checks, typed signed traces, matched PBAS/PLOC/PCON
 * generations, no assumed FWO_ channel order, length consistency checks,
 * no extraction of personal instrument/sample metadata, reverse-complement
 * of sequence/quality/peaks/channels as one coherent record.
 */

(function (global) {
"use strict";

const ABIF_COMPLEMENT = Object.freeze({
  A: 'T', T: 'A', C: 'G', G: 'C', R: 'Y', Y: 'R', S: 'S', W: 'W',
  K: 'M', M: 'K', B: 'V', V: 'B', D: 'H', H: 'D', N: 'N'
});

/** Parse ABI analyzed traces and the original basecaller calls (PBAS2).
 * Falls back to generation 1 only when the whole PBAS2/PLOC2 pair is absent.
 * A generation's quality scores are never mixed with another generation.
 * All per-base and per-trace coordinates are zero based internally.
 */
function parseAb1(buffer, fileName = 'sequence.ab1') {
  if (ArrayBuffer.isView(buffer)) {
    buffer = buffer.buffer.slice(buffer.byteOffset, buffer.byteOffset + buffer.byteLength);
  }
  if (!(buffer instanceof ArrayBuffer) || buffer.byteLength < 34) {
    throw new Error('File AB1 vuoto, troppo corto o non valido.');
  }
  if (buffer.byteLength > 64 * 1024 * 1024) {
    throw new Error('File AB1 oltre il limite di 64 MiB di questa versione.');
  }
  const view = new DataView(buffer);
  const warnings = [];
  function check(offset, length, what) {
    if (!Number.isSafeInteger(offset) || !Number.isSafeInteger(length) ||
        offset < 0 || length < 0 || offset + length > buffer.byteLength) {
      throw new Error(`File AB1 incompleto o corrotto: ${what}.`);
    }
  }
  function ascii(offset, length) {
    check(offset, length, 'testo oltre la fine del file');
    let result = '';
    for (let i = 0; i < length; i++) result += String.fromCharCode(view.getUint8(offset + i));
    return result;
  }
  function entry(offset) {
    check(offset, 28, 'directory');
    return {
      tag: ascii(offset, 4), number: view.getUint32(offset + 4, false),
      type: view.getUint16(offset + 8, false), size: view.getUint16(offset + 10, false),
      count: view.getUint32(offset + 12, false), bytes: view.getUint32(offset + 16, false),
      dataOffset: view.getUint32(offset + 20, false), inlineOffset: offset + 20
    };
  }
  function location(e) {
    const at = e.bytes <= 4 ? e.inlineOffset : e.dataOffset;
    check(at, e.bytes, `${e.tag}${e.number}`);
    return at;
  }
  function string(e) {
    let at = location(e), size = e.bytes;
    if (e.type === 18) {
      if (!size) return '';
      const length = view.getUint8(at);
      if (length + 1 > size) throw new Error('Stringa AB1 Pascal incompleta.');
      at++; size = length;
    }
    return ascii(at, size).replace(/\0+$/g, '');
  }
  function numbers(e, byteQuality = false) {
    const at = location(e);
    const supportedSizes = {1: 1, 2: 1, 3: 2, 4: 2, 5: 4, 7: 4, 8: 8};
    const size = supportedSizes[e.type];
    if (!size || e.size !== size || e.count * size > e.bytes) {
      throw new Error(`Formato numerico non supportato o incompleto: ${e.tag}${e.number}.`);
    }
    if (byteQuality && size !== 1) throw new Error('Formato PCON non supportato: attesi byte.');
    const out = new Array(e.count);
    for (let i = 0; i < e.count; i++) {
      const p = at + i * size;
      if (size === 1) out[i] = view.getUint8(p);
      else if (e.type === 3) out[i] = view.getUint16(p, false);
      else if (e.type === 4) out[i] = view.getInt16(p, false);
      else if (e.type === 5) out[i] = view.getInt32(p, false);
      else if (e.type === 7) out[i] = view.getFloat32(p, false);
      else if (e.type === 8) out[i] = view.getFloat64(p, false);
      if (!Number.isFinite(out[i])) throw new Error('Dato numerico AB1 non finito.');
    }
    return out;
  }
  if (ascii(0, 4) !== 'ABIF') throw new Error('Firma ABIF assente: selezionare un file .ab1 valido.');
  const version = view.getUint16(4, false);
  if (![1, 3].includes(Math.floor(version / 100))) {
    throw new Error(`Versione ABIF ${version} non supportata.`);
  }
  const root = entry(6);
  if (!root.count || root.count > 100000) throw new Error('Numero di elementi ABIF non valido.');
  check(root.dataOffset, root.count * 28, 'tabella delle directory');
  const entries = new Map();
  for (let i = 0; i < root.count; i++) {
    const e = entry(root.dataOffset + i * 28), key = `${e.tag}${e.number}`;
    if (entries.has(key)) throw new Error(`Tag ABIF duplicato: ${key}.`);
    entries.set(key, e);
  }
  const fwo = entries.get('FWO_1');
  if (!fwo) throw new Error('FWO_1 assente: ordine dei canali non determinabile.');
  const channelOrder = string(fwo).toUpperCase();
  if (channelOrder.length !== 4 || [...channelOrder].sort().join('') !== 'ACGT') {
    throw new Error('Ordine dei quattro canali FWO_1 non valido.');
  }
  const channels = {};
  for (let i = 0; i < 4; i++) {
    const e = entries.get(`DATA${i + 9}`);
    if (!e) throw new Error(`Traccia analizzata DATA${i + 9} assente.`);
    channels[channelOrder[i]] = Float32Array.from(numbers(e));
  }
  const sampleCount = channels.A.length;
  if (!sampleCount || !['C', 'G', 'T'].every(b => channels[b].length === sampleCount)) {
    throw new Error('I quattro canali del cromatogramma hanno lunghezze diverse o sono vuoti.');
  }
  let generation;
  if (entries.has('PBAS2') && entries.has('PLOC2')) generation = 2;
  else if (!entries.has('PBAS2') && !entries.has('PLOC2') && entries.has('PBAS1') && entries.has('PLOC1')) {
    generation = 1;
    warnings.push('PBAS2/PLOC2 assenti: letti PBAS1/PLOC1, che possono contenere modifiche precedenti.');
  } else throw new Error('Coppia coerente PBAS/PLOC assente: file AB1 non analizzato o incompleto.');
  const sequence = string(entries.get(`PBAS${generation}`)).toUpperCase();
  if (!sequence.length || /[^ACGTRYSWKMBDHVN]/.test(sequence)) {
    throw new Error('Sequenza ABI vuota o contenente caratteri non IUPAC.');
  }
  const peakPositions = numbers(entries.get(`PLOC${generation}`));
  if (peakPositions.length !== sequence.length) throw new Error('Numero di basi e coordinate dei picchi non coincidente.');
  if (peakPositions.some((x, i) => !Number.isInteger(x) || x < 0 || x >= sampleCount || (i > 0 && x < peakPositions[i - 1]))) {
    throw new Error('Coordinate dei picchi fuori intervallo o non ordinate.');
  }
  const qualityEntry = entries.get(`PCON${generation}`);
  const qualities = qualityEntry ? numbers(qualityEntry, true) : null;
  if (qualities && qualities.length !== sequence.length) throw new Error('Numero di qualità PCON e di basi non coincidente.');
  if (!qualities) warnings.push(`PCON${generation} assente: qualità e trimming automatico non disponibili.`);
  return {
    format: 'ab1', fileName, abiVersion: version, sampleCount, channelOrder,
    channels, sequence, baseCalls: sequence.split(''), peakPositions, qualities,
    sourceIndices: Array.from({length: sequence.length}, (_, i) => i),
    orientation: 'original',
    tags: {baseCalls: `PBAS${generation}`, peaks: `PLOC${generation}`, qualities: qualities ? `PCON${generation}` : null},
    warnings
  };
}

/** Reverse the plotted signal AND complement its channel labels.
 * Merely reversing a sequence string would leave a scientifically wrong trace.
 * This operation is an involution: applying it twice restores all data.
 */
function reverseComplementTrace(record) {
  const channels = {};
  for (const base of ['A', 'C', 'G', 'T']) channels[base] = Float32Array.from(record.channels[ABIF_COMPLEMENT[base]]).reverse();
  const baseCalls = record.baseCalls.slice().reverse().map(b => ABIF_COMPLEMENT[b]);
  return {
    ...record,
    channels,
    baseCalls,
    sequence: baseCalls.join(''),
    peakPositions: record.peakPositions.slice().reverse().map(x => record.sampleCount - 1 - x),
    qualities: record.qualities ? record.qualities.slice().reverse() : null,
    sourceIndices: record.sourceIndices.slice().reverse(),
    orientation: record.orientation === 'reverse-complement' ? 'original' : 'reverse-complement'
  };
}

const api = Object.freeze({parseAb1, reverseComplementTrace});
global.SangerAbif = api;
if (typeof module !== "undefined" && module.exports) module.exports = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
