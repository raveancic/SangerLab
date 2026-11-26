'use strict';
const assert = require('node:assert/strict');
const {test} = require('node:test');
const C = require('../src/core.js');

function dna(n,seed=137) {
  let x=seed>>>0,s='';
  for(let i=0;i<n;++i){x=(Math.imul(x,1664525)+1013904223)>>>0;s+='ACGT'[(x>>>28)&3];}
  return s;
}
const ref=dna(240);
const opts={circular:false,orientation:'forward',minReadLength:1};
const q=(s,n=40)=>Array(s.length).fill(n);

test('FASTA single record and IUPAC validation',()=>{
  const parsed=C.parseFasta('\uFEFF>rCRS_15990 comment\nacgt\nNNRY\n');
  assert.equal(parsed.sequence,'ACGTNNRY');assert.equal(parsed.name,'rCRS_15990');assert.equal(parsed.warnings.length,1);
  for(const s of ['ACGT','>a\nAC\n>b\nGT','>a\nAC-G','>a\nACU','>a\n']) assert.throws(()=>C.parseFasta(s));
});
test('reverse complement preserves every IUPAC symbol',()=>{
  const s='ACGTRYSWKMBDHVN';assert.equal(C.reverseComplement(C.reverseComplement(s)),s);
  assert.equal(C.reverseComplement('AaTC'),'GATT');
});
test('Mott retains internal low-quality bases and trims poor tails',()=>{
  const quality=[...Array(5).fill(0),...Array(40).fill(40),5,...Array(40).fill(40),...Array(4).fill(0)];
  const t=C.mottTrim(quality,.01);assert.equal(t.start,5);assert.equal(t.end,86);assert.equal(t.length,81);
  assert.equal(C.mottTrim(Array(10).fill(0)).empty,true);
  assert.deepEqual(Object.fromEntries(['start','end'].map(k=>[k,C.mottTrim([20,40,20])[k]])),{start:1,end:2});
  assert.throws(()=>C.mottTrim([40,NaN]));assert.throws(()=>C.mottTrim([40],0));
});
test('forward exact alignment retains all read positions and qualities',()=>{
  const s=ref.slice(35,115),a=C.alignRead(s,q(s),ref,{...opts,trimStart:7,origin:100});
  assert.equal(a.score,160);assert.equal(a.summary.fastaStart,36);assert.equal(a.summary.fastaEnd,115);
  assert.equal(a.summary.rcrsStart,135);assert.equal(a.summary.rcrsEnd,214);assert.equal(a.summary.identity,1);
  assert.equal(a.columns[0].readPos,8);assert.equal(a.columns.at(-1).readPos,87);assert.equal(a.variants.length,0);
  assert.equal(a.summary.accepted,true);assert.equal(a.summary.readCoverage,1);
});
test('automatic reverse-complement reverses quality and original read coordinates',()=>{
  const s=C.reverseComplement(ref.slice(40,120)),quality=Array.from({length:s.length},(_,i)=>30+(i%20));
  const a=C.alignRead(s,quality,ref,{...opts,orientation:'auto',trimStart:10});
  assert.equal(a.orientation,'reverse');assert.equal(a.sequence,ref.slice(40,120));
  assert.equal(a.columns[0].readPos,90);assert.equal(a.columns.at(-1).readPos,11);
  assert.equal(a.columns[0].quality,quality.at(-1));assert.equal(a.columns.at(-1).quality,quality[0]);
});
test('circular target crosses the FASTA end without clipping',()=>{
  const s=ref.slice(220)+ref.slice(0,40),a=C.alignRead(s,q(s),ref,{...opts,circular:true,origin:150});
  assert.equal(a.summary.fastaStart,221);assert.equal(a.summary.fastaEnd,40);assert.equal(a.summary.referenceSpan,60);
  assert.equal(a.summary.wrapsFasta,true);assert.equal(a.summary.identity,1);assert.equal(a.columns.length,s.length);
  assert.equal(a.summary.rcrsStart,130);assert.equal(a.summary.rcrsEnd,189);
});
test('rCRS numbering wraps inside a rotated FASTA',()=>{
  const s=ref.slice(80,130),a=C.alignRead(s,q(s),ref,{...opts,circular:true,origin:150});
  assert.equal(a.summary.rcrsStart,230);assert.equal(a.summary.rcrsEnd,39);assert.equal(a.summary.wrapsRcrs,true);
  assert.equal(a.summary.wrapsFasta,false);
  const coordinates=a.columns.map(c=>c.rcrsPos);assert.deepEqual(coordinates.slice(9,13),[239,240,1,2]);
});
test('SNP quality flag and ambiguous read are separate evidence',()=>{
  const s=ref.slice(25,125).split('');s[30]='ACGT'.replace(s[30],'')[0];s[50]='N';
  const qual=q(s);qual[30]=12;const a=C.alignRead(s.join(''),qual,ref,opts);
  const snp=a.variants.find(v=>v.type==='SNP');assert.equal(snp.fastaPos,56);assert.equal(snp.quality,12);assert.equal(snp.lowQuality,true);
  const amb=a.variants.find(v=>v.type==='AMBIGUOUS');assert.equal(amb.fastaPos,76);assert.equal(amb.alt,'N');
  assert.equal(a.summary.ambiguous,1);
});
test('N in the reference is never called as a definite SNP',()=>{
  const target=ref.slice(0,70)+'N'+ref.slice(71);const s=ref.slice(20,120),a=C.alignRead(s,q(s),target,opts);
  assert.equal(a.variants.length,1);assert.equal(a.variants[0].type,'AMBIGUOUS');assert.equal(a.variants[0].ref,'N');
});
test('insertion has two flanks and does not borrow a reference position',()=>{
  const s=ref.slice(20,70)+'GAG'+ref.slice(70,120),a=C.alignRead(s,q(s),ref,opts);
  const v=a.variants.find(v=>v.type==='INS');assert.ok(v);assert.equal(v.alt.length,3);assert.equal(v.positions.length,0);
  assert.equal(v.leftAnchor.rcrsPos+1,v.rightAnchor.rcrsPos);assert.equal(a.summary.inserted,3);
  assert.equal(a.columns.filter(c=>c.op==='I').every(c=>c.rcrsPos===null),true);
});
test('deletion has reference positions but no invented base quality',()=>{
  const s=ref.slice(20,70)+ref.slice(73,120),a=C.alignRead(s,q(s),ref,opts);
  const v=a.variants.find(v=>v.type==='DEL');assert.ok(v);assert.equal(v.ref.length,3);assert.equal(v.positions.length,3);
  assert.equal(v.quality,null);assert.equal(v.flankQuality,40);assert.equal(a.summary.deleted,3);
});
test('homopolymer indel placement is explicitly ambiguous',()=>{
  const target=dna(45,5)+'CCCCCCCC'+dna(45,7),s=target.slice(10,49)+target.slice(50,88);
  const a=C.alignRead(s,q(s),target,opts),v=a.variants.find(x=>x.type==='DEL');
  assert.ok(v);assert.equal(v.placementAmbiguous,true);assert.match(v.note,/non normalizzato/);
});
test('nonhomologous reads are rejected by identity without local clipping',()=>{
  const s='A'.repeat(80),a=C.alignRead(s,q(s),ref,{...opts,minReadLength:50});
  assert.equal(a.summary.accepted,false);assert.equal(a.columns.filter(c=>c.readBase!=='-').length,s.length);
});
test('input and memory guards fail explicitly',()=>{
  assert.throws(()=>C.alignRead('ACGT',[40],ref,opts));
  assert.throws(()=>C.alignRead('ACGT',q('ACGT'),ref,{...opts,origin:0}));
  assert.throws(()=>C.alignRead('ACGT',q('ACGT'),ref,{...opts,maxCells:20}));
  assert.throws(()=>C.alignRead(ref+ref,q(ref+ref),ref,{circular:true}));
});

// Independent full-matrix affine scoring (no traceback, no rolling rows).
function fullScore(a,b) {
  const states=Array.from({length:a.length+1},()=>Array.from({length:b.length+1},()=>[-Infinity,-Infinity,-Infinity]));
  for(let j=0;j<=b.length;++j) states[0][j][0]=0;
  for(let i=1;i<=a.length;++i) {
    states[i][0][1]=-6-(i-1);
    for(let j=1;j<=b.length;++j) {
      const v=states[i][j];
      v[0]=Math.max(...states[i-1][j-1])+(a[i-1]===b[j-1]?2:-3);
      v[1]=Math.max(states[i-1][j][0]-6,states[i-1][j][1]-1,states[i-1][j][2]-6);
      v[2]=Math.max(states[i][j-1][0]-6,states[i][j-1][1]-6,states[i][j-1][2]-1);
    }
  }
  return Math.max(...states[a.length].flatMap(v=>v.slice(0,2)));
}
test('rolling DP agrees with independently indexed full matrices on varied inputs',()=>{
  for(let seed=1;seed<=80;++seed) {
    const a=dna(1+(seed%18),seed),b=dna(8+(seed%24),seed*13);
    const result=C.alignRead(a,q(a),b,opts);assert.equal(result.score,fullScore(a,b));
    assert.equal(result.alignment.query.replace(/-/g,''),a);
  }
});
