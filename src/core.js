/* SangerLab scientific core. MIT License; see LICENSE.
 * Coordinates: internal indices are zero-based; exported positions are one-based.
 * No network requests, external dependencies, read clipping, or reference imputation.
 */
(function (root) {
  'use strict';

  const VERSION = '0.1.3';
  const DNA = /^[ACGTRYSWKMBDHVN]+$/;
  const COMPLEMENT = { A:'T', C:'G', G:'C', T:'A', R:'Y', Y:'R', S:'S', W:'W', K:'M', M:'K', B:'V', V:'B', D:'H', H:'D', N:'N' };
  const NEG = -1e15;

  function cleanSequence(sequence, label) {
    if (typeof sequence !== 'string') throw new Error(label + ': sequenza non valida.');
    const value = sequence.toUpperCase();
    if (!DNA.test(value)) throw new Error(label + ': usare soltanto basi DNA IUPAC, senza gap.');
    return value;
  }

  function parseFasta(text) {
    if (typeof text !== 'string') throw new Error('FASTA: contenuto non testuale.');
    const lines = text.replace(/^\uFEFF/, '').split(/\r?\n/);
    let header = null;
    const parts = [];
    for (const raw of lines) {
      const line = raw.trim();
      if (!line || line.startsWith(';')) continue;
      if (line.startsWith('>')) {
        if (header !== null) throw new Error('Caricare una reference FASTA con una sola sequenza.');
        header = line.slice(1).trim();
        if (!header) throw new Error('FASTA: intestazione vuota.');
      } else {
        if (header === null) throw new Error('FASTA: manca la prima intestazione >nome.');
        parts.push(line.replace(/\s/g, ''));
      }
    }
    if (header === null || !parts.length) throw new Error('FASTA: nessuna sequenza presente.');
    const sequence = cleanSequence(parts.join(''), 'FASTA');
    const name = header.split(/\s+/)[0];
    const ambiguous = (sequence.match(/[^ACGT]/g) || []).length;
    return { name, description:header, sequence, length:sequence.length,
      warnings: ambiguous ? [`La reference contiene ${ambiguous} basi IUPAC ambigue; queste posizioni non producono SNP definiti.`] : [] };
  }

  function reverseComplement(sequence) {
    const value = cleanSequence(sequence, 'Reverse complement');
    let result = '';
    for (let i=value.length-1; i>=0; --i) result += COMPLEMENT[value[i]];
    return result;
  }

  function validateQualities(qualities, length) {
    if (!qualities || qualities.length !== length) throw new Error('Il numero di qualità Phred deve corrispondere alle basi.');
    for (let i=0; i<length; ++i) {
      if (!Number.isFinite(qualities[i]) || qualities[i] < 0 || qualities[i] > 255) {
        throw new Error('Qualità Phred non valida alla base ' + (i+1) + '.');
      }
    }
  }

  /** Modified Mott/Kadane: maximum positive contiguous sum(cutoff - 10^(-Q/10)).
   * Equal-score ties retain the earliest finished segment; exact zero terminal
   * increments are omitted. A nonpositive optimum returns an empty interval.
   * This is NOT a filter requiring every retained base to exceed a Q threshold.
   */
  function mottTrim(qualities, cutoff=0.01) {
    if (!Number.isFinite(cutoff) || cutoff <= 0 || cutoff >= 1) throw new Error('Il cutoff Mott deve essere compreso tra 0 e 1, estremi esclusi.');
    validateQualities(qualities, qualities ? qualities.length : -1);
    let running=0, begin=0, best=0, start=0, end=0;
    for (let i=0; i<qualities.length; ++i) {
      const increment = cutoff - Math.pow(10, -qualities[i]/10);
      if (running <= 1e-12) { running=0; begin=i; }
      running += increment;
      if (running > best + 1e-12) { best=running; start=begin; end=i+1; }
      if (running <= 1e-12) { running=0; begin=i+1; }
    }
    return { start, end, length:end-start, score:best, cutoff, empty:end===start };
  }

  function coordinate(refIndex0, origin, length) {
    const index = ((refIndex0 % length) + length) % length;
    return { fastaPos:index+1, rcrsPos:((origin-1+index)%length)+1 };
  }

  function substitution(a,b,match,mismatch) {
    // Ambiguous symbols carry no mismatch claim and no positive match evidence.
    if (a==='N' || b==='N' || !'ACGT'.includes(a) || !'ACGT'.includes(b)) return 0;
    return a===b ? match : mismatch;
  }

  /** Gotoh-style affine DP, query-global/reference-local. Gaps cost
   * gapOpen + (length-1)*gapExtend. Reference prefix/suffix are free; every
   * query base participates. Traceback uses one byte/cell and rolling scores.
   */
  function semiglobal(query, target, parameters) {
    const m=query.length, n=target.length, stride=n+1;
    const cells=(m+1)*stride;
    if (cells > parameters.maxCells) throw new Error('Allineamento troppo grande: ridurre la lettura/reference. Limite ' + parameters.maxCells.toLocaleString('it-IT') + ' celle.');
    const trace = new Uint8Array(cells);
    let pm=new Float64Array(stride), pi=new Float64Array(stride), pd=new Float64Array(stride);
    let cm=new Float64Array(stride), ci=new Float64Array(stride), cd=new Float64Array(stride);
    pi.fill(NEG); pd.fill(NEG);
    const go=parameters.gapOpen, ge=parameters.gapExtend;
    for (let i=1; i<=m; ++i) {
      cm[0]=NEG; cd[0]=NEG; ci[0]=-go-(i-1)*ge;
      trace[i*stride]=(i===1 ? 0 : 1)<<2;
      const a=query[i-1];
      const knownA = a==='A'||a==='C'||a==='G'||a==='T';
      const row=i*stride;
      for (let j=1; j<=n; ++j) {
        let predecessor=0, best=pm[j-1];
        if (pi[j-1]>best) {best=pi[j-1]; predecessor=1;}
        if (pd[j-1]>best) {best=pd[j-1]; predecessor=2;}
        const b=target[j-1];
        const knownB=b==='A'||b==='C'||b==='G'||b==='T';
        cm[j]=best + (knownA && knownB ? (a===b ? parameters.match : parameters.mismatch) : 0);
        let bits=predecessor;
        best=pi[j]-ge; predecessor=1;
        if (pm[j]-go>best) {best=pm[j]-go; predecessor=0;}
        if (pd[j]-go>best) {best=pd[j]-go; predecessor=2;}
        ci[j]=best; bits|=predecessor<<2;
        best=cd[j-1]-ge; predecessor=2;
        if (cm[j-1]-go>best) {best=cm[j-1]-go; predecessor=0;}
        if (ci[j-1]-go>best) {best=ci[j-1]-go; predecessor=1;}
        cd[j]=best; bits|=predecessor<<4;
        trace[row+j]=bits;
      }
      [pm,cm]=[cm,pm]; [pi,ci]=[ci,pi]; [pd,cd]=[cd,pd];
    }
    let score=NEG, end=0, state=0;
    let endpointCount=0;
    // A final reference deletion is never needed when its suffix is free.
    for (let j=0; j<=n; ++j) {
      const value=Math.max(pm[j],pi[j]);
      if (value>score) {score=value; end=j; state=pm[j]>=pi[j]?0:1; endpointCount=1;}
      else if (value===score) ++endpointCount;
    }
    const uniqueEndpoints=new Set();
    for(let j=0;j<=n;++j) if(Math.max(pm[j],pi[j])===score) {
      uniqueEndpoints.add(parameters.circular ? j%parameters.referenceLength : j);
    }
    let i=m, j=end;
    const reversed=[];
    while (i>0) {
      const bits=trace[i*stride+j];
      if (state===0) {
        if (j===0) throw new Error('Traceback non valido: stato match senza reference.');
        reversed.push({refIndex0:j-1,queryIndex0:i-1});
        state=bits&3; --i; --j;
      } else if (state===1) {
        reversed.push({refIndex0:null,queryIndex0:i-1});
        state=(bits>>2)&3; --i;
      } else if (state===2) {
        if (j===0) throw new Error('Traceback non valido: delezione senza reference.');
        reversed.push({refIndex0:j-1,queryIndex0:null});
        state=(bits>>4)&3; --j;
      } else throw new Error('Traceback non valido.');
    }
    reversed.reverse();
    return {score,start:j,end,rawColumns:reversed,endpointCount,uniqueEndpointCount:uniqueEndpoints.size};
  }

  function anchorFor(index, reference, origin, circular, observed) {
    if (index===null || index===undefined || (!circular && (index<0 || index>=reference.length))) return null;
    const coords=coordinate(index,origin,reference.length);
    return {...coords, base:reference[coords.fastaPos-1], observed};
  }

  function makeVariants(columns, reference, parameters) {
    const variants=[];
    for (let i=0; i<columns.length; ++i) {
      const c=columns[i];
      if (c.op==='=') continue;
      if (c.op==='X' || c.op==='?') {
        variants.push({type:c.op==='X'?'SNP':'AMBIGUOUS', rcrsPos:c.rcrsPos,fastaPos:c.fastaPos,
          ref:c.refBase, alt:c.readBase, readPos:c.readPos, quality:c.quality, flankQuality:null,
          lowQuality:c.lowQuality, positions:[c.rcrsPos], fastaPositions:[c.fastaPos], readPositions:[c.readPos],
          leftAnchor:null,rightAnchor:null,placementAmbiguous:false,
          note:c.op==='?'?'Base IUPAC ambigua nella read o nella reference; nessuna chiamata di SNP definita.':(c.lowQuality?'Differenza in una base di bassa qualità: controllare il cromatogramma.':'Differenza rispetto alla reference; controllare il cromatogramma.')});
        continue;
      }
      const start=i, op=c.op;
      while (i+1<columns.length && columns[i+1].op===op) ++i;
      const block=columns.slice(start,i+1);
      let left=null,right=null;
      for (let k=start-1;k>=0;--k) if(columns[k].refIndex0!==null) {left=columns[k].refIndex0;break;}
      for (let k=i+1;k<columns.length;++k) if(columns[k].refIndex0!==null) {right=columns[k].refIndex0;break;}
      // Flanks outside aligned coverage are reference context, explicitly labelled.
      const hasLeft=left!==null, hasRight=right!==null;
      if(op==='D') {if(left===null) left=block[0].refIndex0-1; if(right===null) right=block.at(-1).refIndex0+1;}
      else {if(left===null && right!==null) left=right-1; if(right===null && left!==null) right=left+1;}
      const leftAnchor=anchorFor(left,reference,parameters.origin,parameters.circular,hasLeft);
      const rightAnchor=anchorFor(right,reference,parameters.origin,parameters.circular,hasRight);
      const ref=op==='D'?block.map(x=>x.refBase).join(''):'';
      const alt=op==='I'?block.map(x=>x.readBase).join(''):'';
      const allele=ref||alt;
      const placementAmbiguous=!!((leftAnchor && leftAnchor.base===allele.at(-1)) || (rightAnchor && rightAnchor.base===allele[0]));
      const qs=block.map(x=>x.quality).filter(x=>x!==null);
      const flankQs=[];
      for(let k=start-1;k>=0;--k) if(columns[k].quality!==null) {flankQs.push(columns[k].quality);break;}
      for(let k=i+1;k<columns.length;++k) if(columns[k].quality!==null) {flankQs.push(columns[k].quality);break;}
      const quality=qs.length?Math.min(...qs):null;
      const flankQuality=flankQs.length?Math.min(...flankQs):null;
      const lowQuality=quality!==null?quality<parameters.qualityThreshold:flankQuality===null||flankQuality<parameters.qualityThreshold;
      const positions=block.filter(x=>x.rcrsPos!==null).map(x=>x.rcrsPos);
      const fastaPositions=block.filter(x=>x.fastaPos!==null).map(x=>x.fastaPos);
      variants.push({type:op==='I'?'INS':'DEL',rcrsPos:positions[0] ?? leftAnchor?.rcrsPos ?? rightAnchor?.rcrsPos ?? null,
        fastaPos:fastaPositions[0] ?? leftAnchor?.fastaPos ?? rightAnchor?.fastaPos ?? null,
        ref,alt,readPos:block.find(x=>x.readPos!==null)?.readPos??null,quality,flankQuality,lowQuality,positions,fastaPositions,
        readPositions:block.filter(x=>x.readPos!==null).map(x=>x.readPos), leftAnchor,rightAnchor,placementAmbiguous,
        note:'Indel rappresentato come nell’allineamento, non normalizzato.'+(placementAmbiguous?' Posizionamento equivalente possibile in una ripetizione.':' Nessuna equivalenza locale semplice rilevata; controllare comunque i picchi.')+(op==='D'?' La delezione non ha un proprio Phred; flankQuality riguarda le basi adiacenti.':'')});
    }
    return variants;
  }

  function alignRead(sequence, qualities, reference, options={}) {
    const query=cleanSequence(sequence,'Read');
    const target=cleanSequence(reference,'Reference');
    validateQualities(qualities,query.length);
    const parameters={origin:1,circular:true,orientation:'auto',match:2,mismatch:-3,gapOpen:6,gapExtend:1,
      trimStart:0,qualityThreshold:20,minIdentity:0.80,minReadLength:50,maxCells:100000000,...options};
    const L=target.length;
    parameters.referenceLength=L;
    if(!Number.isInteger(parameters.origin)||parameters.origin<1||parameters.origin>L) throw new Error('La coordinata rCRS della prima base deve essere un intero da 1 alla lunghezza della reference.');
    if(!['auto','forward','reverse'].includes(parameters.orientation)) throw new Error('Orientamento non valido.');
    if(!Number.isInteger(parameters.trimStart)||parameters.trimStart<0) throw new Error('trimStart deve essere un indice zero-based non negativo.');
    for(const key of ['match','mismatch','gapOpen','gapExtend','qualityThreshold','minIdentity','minReadLength','maxCells']) if(!Number.isFinite(parameters[key])) throw new Error('Parametro numerico non valido: '+key);
    if(parameters.match<=0||parameters.mismatch>0||parameters.gapOpen<=0||parameters.gapExtend<=0) throw new Error('Punteggi di allineamento non validi.');
    if(parameters.minIdentity<0||parameters.minIdentity>1||parameters.minReadLength<1||parameters.qualityThreshold<0) throw new Error('Soglie di controllo non valide.');
    if(parameters.circular && query.length>L) throw new Error('Per questa versione la read deve essere più corta o uguale alla reference circolare.');
    const doubled=parameters.circular?target+target:target;
    let best=null,otherScore=null,selectedOrientation=null,orientedQuery=null,orientedQualities=null;
    const orientations=parameters.orientation==='auto'?['forward','reverse']:[parameters.orientation];
    for(const orientation of orientations) {
      const seq=orientation==='forward'?query:reverseComplement(query);
      const result=semiglobal(seq,doubled,parameters);
      if(best===null||result.score>best.score) {
        otherScore=best?.score??null;
        best=result;selectedOrientation=orientation;orientedQuery=seq;
        orientedQualities=orientation==='forward'?Array.from(qualities):Array.from(qualities).reverse();
      } else otherScore=result.score;
    }
    const referenceSpan=best.end-best.start;
    if(parameters.circular && referenceSpan>L) throw new Error('Allineamento esteso oltre un giro della reference: verificare read, trimming e parametri.');
    const columns=best.rawColumns.map(({refIndex0,queryIndex0})=>{
      const coords=refIndex0===null?{fastaPos:null,rcrsPos:null}:coordinate(refIndex0,parameters.origin,L);
      const refBase=refIndex0===null?'-':doubled[refIndex0];
      const readBase=queryIndex0===null?'-':orientedQuery[queryIndex0];
      const readIndex0=queryIndex0===null?null:parameters.trimStart+(selectedOrientation==='forward'?queryIndex0:query.length-1-queryIndex0);
      const quality=queryIndex0===null?null:orientedQualities[queryIndex0];
      let op;
      if(refIndex0===null) op='I';
      else if(queryIndex0===null) op='D';
      else if(!'ACGT'.includes(refBase)||!'ACGT'.includes(readBase)) op='?';
      else op=refBase===readBase?'=':'X';
      return {op,refIndex0:refIndex0===null?null:refIndex0%L,refUnwrappedIndex0:refIndex0,queryIndex0,...coords,
        readIndex0,readPos:readIndex0===null?null:readIndex0+1,refBase,readBase,quality,lowQuality:quality!==null&&quality<parameters.qualityThreshold};
    });
    // Variant context needs unwrapped indices across the circular join.
    const contextColumns=columns.map(c=>({...c,refIndex0:c.refUnwrappedIndex0}));
    const variants=makeVariants(contextColumns,target,parameters);
    const counts={'=':0,'X':0,'I':0,'D':0,'?':0};
    for(const c of columns) ++counts[c.op];
    const identity=counts['=']/Math.max(1,columns.length);
    const lowQuality=orientedQualities.filter(q=>q<parameters.qualityThreshold).length;
    const first=columns.find(c=>c.rcrsPos!==null),last=columns.findLast(c=>c.rcrsPos!==null);
    const warnings=[];
    const accepted=identity>=parameters.minIdentity&&query.length>=parameters.minReadLength&&referenceSpan>0;
    if(query.length<parameters.minReadLength) warnings.push('Read più corta della lunghezza minima di controllo ('+parameters.minReadLength+' basi).');
    if(identity<parameters.minIdentity) warnings.push('Identità inferiore alla soglia di controllo ('+Math.round(parameters.minIdentity*100)+'%): allineamento non accettato.');
    if(lowQuality) warnings.push(lowQuality+' basi mantenute hanno qualità inferiore a Q'+parameters.qualityThreshold+'. Il trimming Mott può mantenere basi interne di bassa qualità.');
    if(counts['?']) warnings.push(counts['?']+' colonne contengono simboli IUPAC ambigui: nessuna chiamata di SNP definita per queste colonne.');
    if(otherScore!==null&&otherScore===best.score) warnings.push('Entrambi gli orientamenti hanno lo stesso punteggio; orientamento da verificare.');
    // A doubled circle has duplicate endpoints; count distinct circular loci.
    // This detects endpoint ambiguity, not every alternative internal gap path.
    if(best.uniqueEndpointCount>1) warnings.push('Più posizionamenti hanno lo stesso punteggio; verificare la specificità della read.');
    if(variants.some(v=>v.placementAmbiguous)) warnings.push('Almeno un indel cade in una ripetizione: la posizione mostrata dipende dall’allineamento e non è normalizzata.');
    if(!referenceSpan) warnings.push('Nessuna base della reference allineata.');
    let wrapsRcrs=false,wrapsFasta=false;
    for(let i=1;i<columns.length;++i) {
      const c=columns[i]; if(c.rcrsPos===null) continue;
      let k=i-1;while(k>=0&&columns[k].rcrsPos===null)--k;
      if(k>=0){if(c.rcrsPos<columns[k].rcrsPos)wrapsRcrs=true;if(c.fastaPos<columns[k].fastaPos)wrapsFasta=true;}
    }
    return {version:VERSION,orientation:selectedOrientation,score:best.score,otherOrientationScore:otherScore,
      sequence:orientedQuery,qualities:orientedQualities,columns,variants,
      alignment:{reference:columns.map(c=>c.refBase).join(''),query:columns.map(c=>c.readBase).join(''),
        matchLine:columns.map(c=>c.op==='='?'|':c.op==='?'?'?':' ').join('')},
      summary:{readLength:query.length,alignmentColumns:columns.length,referenceLength:L,referenceSpan,
        fastaStart:first?.fastaPos??null,fastaEnd:last?.fastaPos??null,rcrsStart:first?.rcrsPos??null,rcrsEnd:last?.rcrsPos??null,
        matches:counts['='],mismatches:counts.X,inserted:counts.I,deleted:counts.D,ambiguous:counts['?'],lowQuality,
        identity,readCoverage:1,referenceCoverage:referenceSpan/L,wrapsRcrs,wrapsFasta,accepted,
        equivalentEndpointCount:best.uniqueEndpointCount,
        meanQuality:orientedQualities.reduce((a,b)=>a+b,0)/orientedQualities.length},
      warnings,parameters};
  }

  const api={VERSION,parseFasta,reverseComplement,mottTrim,coordinate,alignRead};
  root.SangerCore=api;
  if(typeof module!=='undefined'&&module.exports) module.exports=api;
})(typeof globalThis!=='undefined'?globalThis:this);
