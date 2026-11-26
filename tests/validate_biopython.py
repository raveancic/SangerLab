#!/usr/bin/env python3
"""Independent SangerLab scientific validation against Biopython.

Requires Python >=3.10, Biopython and Node.js. No sample is embedded.
Run synthetic checks:
    python tests/validate_biopython.py
Add your own data, without uploading them:
    python tests/validate_biopython.py --ab1 sample.ab1 \
        --reference rotated.fasta --canonical NC_012920.1.fasta --origin 15990 \
        --output validation.json

Scoring: +2/-3, affine gap 6 + (length-1), all IUPAC ambiguous pairs score 0.
Full read aligns, reference terminal overhangs are free. Circular reference is
represented twice. Mott is tested against exhaustive contiguous sums, not
Biopython's historical private trimming function.
"""
from __future__ import annotations
import argparse
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
from typing import Any

from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
import Bio

DNA='ACGTRYSWKMBDHVN'
ROOT=Path(__file__).resolve().parents[1]
NODE_DRIVER=r"""
const fs=require('node:fs');
const C=require(process.argv[1]);
const A=require(process.argv[2]);
const input=JSON.parse(fs.readFileSync(0,'utf8'));
const results=input.cases.map(c=>({name:c.name,result:C.alignRead(c.sequence,c.qualities,c.reference,c.options)}));
const mott=input.mott.map(c=>C.mottTrim(c.qualities,c.cutoff));
let abi=null;
if(input.abiPath){const b=fs.readFileSync(input.abiPath);const a=A.parseAb1(b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength),'validation.ab1');
 abi={sequence:a.sequence,qualities:a.qualities,peaks:a.peakPositions,channels:Object.fromEntries(Object.entries(a.channels).map(([k,v])=>[k,Array.from(v)]))};}
process.stdout.write(JSON.stringify({results,mott,abi}));
"""


def parser() -> argparse.ArgumentParser:
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--ab1',type=Path,help='Optional real AB1 file; must also provide --reference.')
    p.add_argument('--reference',type=Path,help='One-record FASTA used by the app.')
    p.add_argument('--canonical',type=Path,help='Optional canonical FASTA to verify rotation and invariant SNP coordinates.')
    p.add_argument('--origin',type=int,default=1,help='Canonical 1-based coordinate of first reference base (default: 1).')
    p.add_argument('--node',default=shutil.which('node') or os.environ.get('CODEX_PRIMARY_RUNTIME_NODE','node'))
    p.add_argument('--output',type=Path,help='Optional JSON report. May contain user sequence results; keep private.')
    return p


def mott_exhaustive(q: list[int],cutoff: float) -> tuple[int,int,float]:
    """O(n^2) oracle, earlier-finished ties then later start after zero prefix."""
    scores=[cutoff-10**(-value/10) for value in q]
    best=0.;start=end=0
    # End-order supplies earlier-finished ties. Later start removes zero prefixes.
    for hi in range(1,len(q)+1):
        running=0.
        for lo in range(hi-1,-1,-1):
            running+=scores[lo]
            if running>best+1e-12:
                best=running;start=lo;end=hi
    return start,end,best


def mott_linear(q: list[int],cutoff: float) -> tuple[int,int]:
    total=best=0.;start=lo=hi=0
    for i,value in enumerate(q):
        total+=cutoff-10**(-value/10)
        if total<=1e-12:total=0.;start=i+1
        elif total>best+1e-12:best=total;lo=start;hi=i+1
    return lo,hi


def make_aligner() -> Align.PairwiseAligner:
    matrix=substitution_matrices.Array(alphabet=DNA,dims=2)
    for a in DNA:
        for b in DNA:matrix[a,b]=(2 if a==b else -3) if a in 'ACGT' and b in 'ACGT' else 0
    x=Align.PairwiseAligner(mode='global',open_gap_score=-6,extend_gap_score=-1)
    x.substitution_matrix=matrix
    # Gaps in query are reference overhangs. Gaps in target remain penalized.
    if 'left_deletion_score' in dir(x):
        x.left_deletion_score=0;x.right_deletion_score=0
    else:
        x.query_left_gap_score=0;x.query_right_gap_score=0
    return x


def ensure(condition: bool,message: str) -> None:
    if not condition:raise AssertionError(message)


def score_alignment(ref: str,query: str) -> float:
    score=0;previous_gap=None
    for a,b in zip(ref,query):
        if a=='-' or b=='-':
            ensure(not(a=='-' and b=='-'),'Double-gap column')
            gap='I' if a=='-' else 'D'
            score-=1 if gap==previous_gap else 6;previous_gap=gap
        else:
            score+=(2 if a==b else -3) if a in 'ACGT' and b in 'ACGT' else 0
            previous_gap=None
    return score


def main() -> int:
    args=parser().parse_args()
    if bool(args.ab1)!=bool(args.reference):parser().error('Use --ab1 and --reference together.')
    if args.canonical and not args.reference:parser().error('--canonical requires --reference.')
    rng=random.Random(137)
    dna=lambda n:''.join(rng.choice('ACGT') for _ in range(n))
    cases=[]
    def add(name: str,s: str,r: str,q: list[int]|None=None,**options: Any) -> None:
        cases.append({'name':name,'sequence':s,'reference':r,'qualities':q or [40]*len(s),
          'options':{'circular':False,'orientation':'auto','minReadLength':1,**options}})
    ref=dna(270)
    add('synthetic_exact',ref[35:145],ref)
    add('synthetic_reverse',str(Seq(ref[35:145]).reverse_complement()),ref,list(range(20,130)),trimStart=11)
    add('synthetic_fasta_join',ref[-65:]+ref[:70],ref,circular=True,origin=160)
    add('synthetic_coordinate_join',ref[90:160],ref,circular=True,origin=160)
    insert=ref[25:80]+'ATG'+ref[80:170]
    delete=ref[25:80]+ref[83:170]
    add('synthetic_insertion',insert,ref)
    add('synthetic_deletion',delete,ref)
    h=dna(50)+'CCCCCCCC'+dna(50)
    add('synthetic_homopolymer',h[10:53]+h[54:98],h)
    mutated=list(ref[30:170]);mutated[35]=next(x for x in 'ACGT' if x!=mutated[35]);mutated[55]='N';mutated[60]='R'
    quality=[40]*len(mutated);quality[35]=12
    add('synthetic_snp_iupac',''.join(mutated),ref,quality,qualityThreshold=30)
    add('synthetic_reference_N',ref[30:170],ref[:85]+'N'+ref[86:])
    add('synthetic_nonhomologous','A'*90,ref,minIdentity=0.80,minReadLength=50)
    for i in range(30):
        add('random_'+str(i),''.join(rng.choice(DNA if i%4==0 else 'ACGT') for _ in range(rng.randint(3,40))),dna(rng.randint(45,80)))
    mott=[{'qualities':q,'cutoff':c} for q in ([0]*10,[40]*10,[20,40,20],[10,0,50,0,50],[],[40,20,40]) for c in (.01,.05)]
    for i in range(120):mott.append({'qualities':[rng.choice((0,10,20,30,40,55)) for _ in range(rng.randrange(2,30))],'cutoff':.01 if i%2 else .05})
    record=None;canonical=None;reference=None;trim_report=[]
    if args.ab1:
        record=SeqIO.read(args.ab1,'abi');reference=str(SeqIO.read(args.reference,'fasta').seq).upper()
        ensure(1<=args.origin<=len(reference),'Origin outside reference')
        if args.canonical:
            canonical=str(SeqIO.read(args.canonical,'fasta').seq).upper()
            ensure(len(canonical)==len(reference),'Canonical and reference lengths differ')
            offset=args.origin-1
            ensure(reference==canonical[offset:]+canonical[:offset],'Reference is not exact rotation at supplied origin')
        seq=str(record.seq).upper();quality=record.letter_annotations['phred_quality']
        add('real_untrimmed',seq,reference,quality,circular=True,origin=args.origin)
        for cutoff in (.01,.05):
            lo,hi=mott_linear(quality,cutoff);ensure(hi>lo,'Real read has no positive Mott interval')
            trim_report.append({'cutoff':cutoff,'start_0':lo,'end_0_exclusive':hi,'length':hi-lo})
            prefix=f'real_p{cutoff:g}'
            add(prefix+'_rotated',seq[lo:hi],reference,quality[lo:hi],circular=True,origin=args.origin,trimStart=lo,qualityThreshold=30)
            add(prefix+'_reverse',str(Seq(seq[lo:hi]).reverse_complement()),reference,quality[lo:hi][::-1],circular=True,origin=args.origin,qualityThreshold=30)
            if canonical:
                add(prefix+'_canonical',seq[lo:hi],canonical,quality[lo:hi],circular=True,origin=1,trimStart=lo,qualityThreshold=30)
        for name,r,origin in [('rotated',reference,args.origin)]+([('canonical',canonical,1)] if canonical else []):
            add('real_reference_'+name+'_join',r[-65:]+r[:70],r,circular=True,origin=origin)
    payload={'cases':cases,'mott':mott,'abiPath':str(args.ab1.resolve()) if args.ab1 else None}
    proc=subprocess.run([args.node,'-e',NODE_DRIVER,str(ROOT/'src/core.js'),str(ROOT/'src/abif.js')],input=json.dumps(payload),text=True,capture_output=True,check=True)
    data=json.loads(proc.stdout)
    for inp,result in zip(mott,data['mott']):
        lo,hi,best=mott_exhaustive(inp['qualities'],inp['cutoff'])
        ensure(result['start']==lo and result['end']==hi,f'Mott interval mismatch: {inp}: JS{result}, oracle{(lo,hi)}')
        ensure(abs(result['score']-best)<1e-9,'Mott score mismatch')
    aligner=make_aligner();checks=[];by_name={}
    for inp,row in zip(cases,data['results']):
        result=row['result'];by_name[inp['name']]=result
        r=inp['reference'];s=inp['sequence'];opts=inp['options'];circular=opts.get('circular',False)
        target=r+r if circular else r
        expected_scores={'forward':aligner.score(target,s),'reverse':aligner.score(target,str(Seq(s).reverse_complement()))}
        ensure(result['score']==max(expected_scores.values()),f'{inp["name"]}: DP score mismatch {result["score"]}, {expected_scores}')
        ensure(result['score']==expected_scores[result['orientation']],inp['name']+': wrong orientation')
        expected_seq=s if result['orientation']=='forward' else str(Seq(s).reverse_complement())
        ensure(result['alignment']['query'].replace('-','')==expected_seq,inp['name']+': lost query bases')
        ensure(score_alignment(result['alignment']['reference'],result['alignment']['query'])==result['score'],inp['name']+': traceback score mismatch')
        seen=[];last=None;origin=opts.get('origin',1);trim=opts.get('trimStart',0)
        for col in result['columns']:
            if col['queryIndex0'] is not None:
                qi=col['queryIndex0'];source=qi if result['orientation']=='forward' else len(s)-1-qi
                ensure(col['quality']==inp['qualities'][source],inp['name']+': quality reversed incorrectly')
                ensure(col['readPos']==trim+source+1,inp['name']+': source position incorrect')
                seen.append(qi)
            else:ensure(col['quality'] is None,inp['name']+': invented deletion quality')
            if col['refUnwrappedIndex0'] is not None:
                ri=col['refUnwrappedIndex0']
                ensure(col['fastaPos']==ri%len(r)+1,inp['name']+': FASTA coordinate incorrect')
                ensure(col['rcrsPos']==(origin-1+ri)%len(r)+1,inp['name']+': canonical coordinate incorrect')
                ensure(col['refBase']==r[ri%len(r)],inp['name']+': reference base incorrect')
                if last is not None:ensure(ri==last+1,inp['name']+': skipped reference without deletion')
                last=ri
            else:ensure(col['rcrsPos'] is None,inp['name']+': insertion given own reference position')
        ensure(seen==list(range(len(s))),inp['name']+': query duplicated/omitted')
        for v in result['variants']:
            if v['type']=='SNP':ensure(v['ref'] in 'ACGT' and v['alt'] in 'ACGT',inp['name']+': ambiguous SNP')
            if v['type']=='DEL':ensure(v['quality'] is None,inp['name']+': deletion given own Phred')
        checks.append({'name':inp['name'],'score':result['score'],'oracle_scores':expected_scores,'orientation':result['orientation'],'summary':result['summary']})
    # Expected specific synthetic behavior is independent of implementation output.
    ensure(by_name['synthetic_fasta_join']['summary']['wrapsFasta'],'Failed FASTA end-to-start join')
    ensure(by_name['synthetic_coordinate_join']['summary']['wrapsRcrs'],'Failed rotated canonical join')
    ensure(not by_name['synthetic_nonhomologous']['summary']['accepted'],'Nonhomologous read accepted')
    snp=next(v for v in by_name['synthetic_snp_iupac']['variants'] if v['type']=='SNP')
    ensure((snp['rcrsPos'],snp['quality'],snp['lowQuality'])==(66,12,True),'Synthetic SNP coordinate/quality failed')
    ensure(any(v['placementAmbiguous'] for v in by_name['synthetic_homopolymer']['variants']),'Homopolymer indel not flagged')
    if record:
        abi=data['abi'];raw=record.annotations['abif_raw'];generation=2
        ensure(abi['sequence']==str(record.seq),'AB1 sequence differs from Biopython')
        ensure(abi['qualities']==record.letter_annotations['phred_quality'],'AB1 Phred differs from Biopython')
        ensure(abi['peaks']==list(raw[f'PLOC{generation}']),'AB1 peaks differ from Biopython')
        for i,base in enumerate(raw['FWO_1'].decode()):ensure(abi['channels'][base]==list(raw[f'DATA{i+9}']),'AB1 trace differs '+base)
        for cutoff in (.01,.05):
            names=[f'real_p{cutoff:g}_rotated',f'real_p{cutoff:g}_reverse']+([f'real_p{cutoff:g}_canonical'] if canonical else [])
            normalized=[]
            for name in names:
                a=by_name[name]
                snps=sorted((v['rcrsPos'],v['ref'],v['alt'],v['quality']) for v in a['variants'] if v['type']=='SNP')
                normalized.append((snps,a['summary']['inserted'],a['summary']['deleted']))
            ensure(all(x==normalized[0] for x in normalized),'SNP/indel quantities changed across reference rotation/reverse complement')
    report={'status':'PASS','biopython_version':Bio.__version__,'alignment_cases':len(cases),'mott_cases':len(mott),'parser_compared':bool(record),'reference_rotation_compared':bool(canonical),'trimming':trim_report,'checks':checks,
       'limits':['Optimal indel placements may differ in repeats; score, read reconstruction, SNP coordinates and gap lengths are checked.','Synthetic browser-independent scientific validation is not a physical Android usability test.']}
    if record:report['real_variants']={name:by_name[name]['variants'] for name in by_name if name.startswith('real_p')}
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('checks','real_variants')},indent=2))
    return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except (AssertionError,ValueError,subprocess.CalledProcessError) as exc:
        print('VALIDATION FAILED: '+str(exc),file=sys.stderr)
        if isinstance(exc,subprocess.CalledProcessError) and exc.stderr:print(exc.stderr,file=sys.stderr)
        raise SystemExit(1)
