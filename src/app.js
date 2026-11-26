/* SangerLab UI, MIT. Inputs and results stay in memory in this browser. */
(function () {
  'use strict';
  const $ = id => document.getElementById(id);
  const Core = globalThis.SangerCore;
  const Parser = globalThis.SangerAbif;
  if(window.matchMedia('(max-width: 680px)').matches) $('view-count').value='30';
  const bundled = JSON.parse($('sanger-reference').textContent);
  const N = bundled.sequence.length;
  const canonical = bundled.sequence.slice(N - bundled.origin + 1) + bundled.sequence.slice(0, N - bundled.origin + 1);
  const colors = {A:'#19854b',C:'#326fbe',G:'#334155',T:'#c94652',N:'#8b7b89',R:'#8b7b89',Y:'#8b7b89'};
  const state = {read:null,reference:null,origin:bundled.origin,result:null,trimMethod:'Mott',trimCutoff:0.01,autoEmpty:false,
    readHash:null,referenceHash:null,viewStart:0,selected:null,alignmentPage:0,worker:null,run:0,load:0,busy:false,mapping:new Map(),elapsed:null};
  const number = n => n === null || n === undefined ? '—' : Number(n).toLocaleString('it-IT');
  const text = (id,value) => { $(id).textContent = value; };
  const element = (tag,content,className) => {const e=document.createElement(tag);if(content!==undefined)e.textContent=content;if(className)e.className=className;return e;};
  function status(message,type='') { text('status',message); $('status').className='status'+(type?' '+type:''); }
  async function sha256(bytes) {
    if(!globalThis.crypto?.subtle) return null;
    try {const a=await crypto.subtle.digest('SHA-256',bytes);return [...new Uint8Array(a)].map(n=>n.toString(16).padStart(2,'0')).join('');} catch {return null;}
  }
  function region() {
    if(!state.read) throw new Error('Carica prima un cromatogramma.');
    const start=Number($('trim-start').value),end=Number($('trim-end').value);
    if(state.autoEmpty||!Number.isInteger(start)||!Number.isInteger(end)||start<1||end<start||end>state.read.sequence.length) throw new Error('Trimming non valido: indica prima e ultima base, incluse, fra 1 e '+state.read.sequence.length+'.');
    return {start:start-1,end,length:end-start+1};
  }
  function toggleControls() {
    const ready=!!state.read;
    for(const id of ['trim-start','trim-end','reset-trim','jump-base','jump-button','prev-trace','next-trace','trace-slider']) $(id).disabled=!ready||state.busy;
    $('auto-trim').disabled=!ready||!state.read.qualities||state.busy;
    let valid=false;try {region();valid=true;} catch {}
    $('analyze').disabled=!ready||!state.read.qualities||!valid||state.busy;
    $('orientation').disabled=state.busy;$('cutoff').disabled=state.busy;
  }
  function invalidate() {
    state.run++;if(state.worker){state.worker.terminate();state.worker=null;}
    state.busy=false;state.result=null;state.mapping.clear();state.elapsed=null;
    $('result-content').hidden=true;$('result-empty').hidden=false;
    text('alignment-state','Da allineare');$('alignment-state').className='pill';
    $('analyze').replaceChildren(document.createTextNode('Allinea alla rCRS →'));
    renderStats();toggleControls();drawTrace();renderInspector();
  }
  function setReference(reference) {
    if(reference.sequence.length!==N) throw new Error('Questa versione richiede la rCRS completa di '+number(N)+' posizioni, originale o ruotata.');
    const offset=(canonical+canonical).indexOf(reference.sequence);
    if(offset<0||offset>=N) throw new Error('Il FASTA non è una rotazione esatta di NC_012920.1. Non posso assegnare automaticamente coordinate rCRS affidabili. La reference precedente resta selezionata.');
    invalidate();state.reference=reference;state.origin=offset+1;
    state.referenceHash=null;
    const thisReference=reference;
    sha256(new TextEncoder().encode(reference.sequence)).then(hash=>{if(state.reference===thisReference)state.referenceHash=hash;});
    text('reference-name',reference.name);text('reference-meta',number(N)+' posizioni · origine rCRS '+number(state.origin));
    const join=N-state.origin+1;
    text('coordinate-example',`Prima base FASTA → rCRS ${number(state.origin)}. `+(state.origin===1?'Il FASTA segue la numerazione canonica.':`Basi FASTA ${number(join)} e ${number(join+1)} → rCRS ${number(N)} e 1.`));
  }
  function trimAuto(announce=true) {
    if(!state.read?.qualities) return;
    const cutoff=Number($('cutoff').value);const trim=Core.mottTrim(state.read.qualities,cutoff);
    invalidate();state.trimMethod='Mott';state.trimCutoff=cutoff;state.autoEmpty=trim.empty;
    $('trim-start').value=trim.start+1;$('trim-end').value=trim.end;
    if(trim.empty)status('Nessun segmento ha un punteggio Mott positivo. Controlla le qualità; puoi impostare manualmente un intervallo per ispezionarlo.','error');
    else if(announce)status(`Trimming Mott p=${cutoff}: conservate le basi ${trim.start+1}–${trim.end}. Controlla i picchi e avvia l'allineamento.`);
    updateTrimInfo();
  }
  function updateTrimInfo() {
    let r=null;try {r=region();} catch {}
    if(r)text('trim-info',`${number(r.length)} basi conservate · ${r.start} escluse a sinistra, ${state.read.sequence.length-r.end} a destra · ${state.trimMethod}`);
    else text('trim-info','Intervallo non valido o nessun segmento conservato.');
    renderStats();toggleControls();drawTrace();renderInspector();
  }
  async function loadRead(file) {
    if(!file)return;const token=++state.load;
    if(file.size>64*1024*1024){status('File oltre 64 MiB: scegli un singolo cromatogramma AB1.','error');return;}
    status('Lettura del file AB1 e delle qualità…','busy');
    try {
      const bytes=await file.arrayBuffer();const parsed=Parser.parseAb1(bytes,file.name);const hash=await sha256(bytes);
      if(token!==state.load)return;
      invalidate();state.read=parsed;state.readHash=hash;state.selected=null;state.autoEmpty=false;
      $('trim-start').max=parsed.sequence.length;$('trim-end').max=parsed.sequence.length;$('jump-base').max=parsed.sequence.length;
      $('trim-start').value=1;$('trim-end').value=parsed.sequence.length;
      text('read-name',file.name);text('read-meta',`${number(parsed.sequence.length)} basi · ${parsed.tags.baseCalls} / ${parsed.tags.qualities||'Phred assenti'} · 4 canali`);
      $('trace-empty').hidden=true;
      if(parsed.qualities){trimAuto(false);state.viewStart=Math.max(0,Number($('trim-start').value)-6);}
      else {state.viewStart=0;state.trimMethod='manuale';updateTrimInfo();}
      if(parsed.warnings.length)status(parsed.warnings.join(' '),parsed.qualities?'':'error');
      else if(!state.autoEmpty)status(`File aperto: ${file.name}. Trimming proposto; puoi regolare le estremità e premere «Allinea alla rCRS».`);
      syncTraceControls();drawTrace();toggleControls();
    } catch(error){if(token===state.load)status(error.message,'error');}
  }
  function renderStats() {
    if(!state.read)return;
    text('stat-length',number(state.read.sequence.length));text('stat-length-note',`Phred ${state.read.qualities?'originali':'non disponibili'} · ${state.read.channelOrder}`);
    let r;try{r=region();}catch{text('stat-trim','—');text('stat-q','—');text('stat-map','—');return;}
    text('stat-trim',number(r.length));text('stat-trim-note',`basi AB1 ${r.start+1}–${r.end}`);
    if(state.read.qualities){const q=state.read.qualities.slice(r.start,r.end);text('stat-q',(100*q.filter(x=>x>=20).length/q.length).toLocaleString('it-IT',{maximumFractionDigits:1})+'%');text('stat-q-note',`Q ≥ 20 · Q media ${(q.reduce((a,b)=>a+b,0)/q.length).toLocaleString('it-IT',{maximumFractionDigits:1})}`);}
    else{text('stat-q','—');text('stat-q-note','Phred assenti nel file');}
    if(state.result){const s=state.result.summary;text('stat-map',`${number(s.rcrsStart)} → ${number(s.rcrsEnd)}`);text('stat-map-note',s.accepted?(s.wrapsRcrs?'attraversa 16.569 → 1':'nell’orientamento della reference'):'ipotesi non accettata: vedi avvisi');}
    else{text('stat-map','—');text('stat-map-note','disponibili dopo l’allineamento');}
  }
  function syncTraceControls() {
    if(!state.read)return;
    const count=Number($('view-count').value);const max=Math.max(0,state.read.sequence.length-count);
    state.viewStart=Math.max(0,Math.min(max,state.viewStart));
    $('trace-slider').max=max;$('trace-slider').value=state.viewStart;
    const end=Math.min(state.viewStart+count,state.read.sequence.length);
    text('trace-range',`Basi AB1 ${state.viewStart+1}–${end} di ${number(state.read.sequence.length)}`);
    text('slider-position',`${state.viewStart+1} / ${state.read.sequence.length}`);
  }
  let plotGeometry=null;
  function drawTrace() {
    const canvas=$('trace-canvas'),rect=canvas.getBoundingClientRect(),dpr=Math.min(window.devicePixelRatio||1,2);
    if(rect.width<1)return;
    canvas.width=Math.round(rect.width*dpr);canvas.height=Math.round(rect.height*dpr);
    const ctx=canvas.getContext('2d');ctx.scale(dpr,dpr);const w=rect.width,h=rect.height;ctx.clearRect(0,0,w,h);
    if(!state.read){plotGeometry=null;return;}
    syncTraceControls();const read=state.read,lo=state.viewStart,hi=Math.min(read.sequence.length,lo+Number($('view-count').value));
    const peaks=read.peakPositions,pad=hi-lo>1?(peaks[hi-1]-peaks[lo])/(hi-lo-1)/2:5;
    const sampleLo=Math.max(0,peaks[lo]-pad),sampleHi=Math.min(read.sampleCount-1,peaks[hi-1]+pad);
    const left=36,right=w-12,width=right-left,scale=width/Math.max(1,sampleHi-sampleLo),x=s=>left+(s-sampleLo)*scale;
    let minimum=0,maximum=1;
    for(const b of ['A','C','G','T'])for(let i=Math.floor(sampleLo);i<=Math.ceil(sampleHi);i++){const value=read.channels[b][i];if(value>maximum)maximum=value;if(value<minimum)minimum=value;}
    let trim={start:0,end:0};try{trim=region();}catch{}
    for(let i=lo;i<hi;i++) {
      const a=x(i? (peaks[i-1]+peaks[i])/2:sampleLo),b=x(i+1<peaks.length?(peaks[i]+peaks[i+1])/2:sampleHi);
      if(i<trim.start||i>=trim.end){ctx.fillStyle='#e5e4eab8';ctx.fillRect(Math.max(left,a),0,Math.min(right,b)-Math.max(left,a),h);}
      if(state.selected===i){ctx.fillStyle='#d8c5f396';ctx.fillRect(Math.max(left,a),0,Math.min(right,b)-Math.max(left,a),h);}
    }
    ctx.fillStyle='#656274';ctx.font='8px system-ui';ctx.textAlign='left';ctx.fillText('AB1',3,15);ctx.fillText('rCRS',3,29);
    ctx.strokeStyle='#e3e0e8';ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(left,211);ctx.lineTo(right,211);ctx.moveTo(left,284);ctx.lineTo(right,284);ctx.stroke();
    const y=value=>209-(value-minimum)/(maximum-minimum)*137;
    ctx.save();ctx.beginPath();ctx.rect(left,65,width,149);ctx.clip();
    for(const b of ['G','C','A','T']){ctx.strokeStyle=colors[b];ctx.lineWidth=1.25;ctx.beginPath();for(let i=Math.floor(sampleLo);i<=Math.ceil(sampleHi);i++){if(i===Math.floor(sampleLo))ctx.moveTo(x(i),y(read.channels[b][i]));else ctx.lineTo(x(i),y(read.channels[b][i]));}ctx.stroke();}
    ctx.restore();
    const qmax=Math.max(60,...(read.qualities||[]));
    ctx.setLineDash([3,3]);ctx.strokeStyle='#b9b2c5';ctx.beginPath();ctx.moveTo(left,284-20/qmax*50);ctx.lineTo(right,284-20/qmax*50);ctx.stroke();ctx.setLineDash([]);
    ctx.textAlign='left';ctx.font='8px system-ui';ctx.fillStyle='#656274';ctx.fillText('Q20',2,284-20/qmax*50+3);
    for(let i=lo;i<hi;i++){
      const px=x(peaks[i]),base=read.sequence[i],q=read.qualities?.[i],mapped=state.mapping.get(i);
      ctx.textAlign='center';ctx.font='8px system-ui';ctx.fillStyle='#656274';
      if(i%5===0||(i===lo&&(5-i%5)%5>=3)){ctx.fillText(String(i+1),px,15);if(mapped)ctx.fillText(mapped.rcrsPos===null?'+':String(mapped.rcrsPos),px,29);}
      ctx.font=hi-lo>80?'10px ui-monospace,monospace':'bold 12px ui-monospace,monospace';ctx.fillStyle=colors[base]||'#8b7b89';ctx.fillText(base,px,53);
      if(mapped&&mapped.op!=='='&&mapped.op!=='?'){ctx.fillStyle='#c1922c';ctx.beginPath();ctx.arc(px,62,2.4,0,2*Math.PI);ctx.fill();}
      if(q!==undefined&&q!==null){const bw=Math.max(2,Math.min(12,width/(hi-lo)*.65));ctx.fillStyle=q<20?'#d89767':'#8d70b8';ctx.fillRect(px-bw/2,284-q/qmax*50,bw,q/qmax*50);if(hi-lo<=50&&(width/(hi-lo)>12||i%5===0)){ctx.fillStyle='#656274';ctx.font='8px system-ui';ctx.fillText(String(q),px,231);}}
    }
    plotGeometry={left,right,sampleLo,sampleHi,lo,hi,width};
  }
  function renderInspector() {
    const box=$('base-inspector');
    if(!state.read||state.selected===null){box.replaceChildren(element('span','Seleziona una base nel cromatogramma. Le barre inferiori rappresentano i Phred score originali; la linea tratteggiata indica Q20.','muted'));return;}
    const i=state.selected,read=state.read,q=read.qualities?.[i],c=state.mapping.get(i);
    const fields=[['Base chiamata',read.sequence[i]],['Posizione AB1',String(i+1)],['Phred originale',q===undefined?'assente':'Q'+q],['P(errore) nominale',q===undefined?'—':(100*Math.pow(10,-q/10)).toLocaleString('it-IT',{maximumSignificantDigits:3})+'%'],['Posizione rCRS',c?(c.rcrsPos===null?'inserzione':String(c.rcrsPos)):'non allineata'],['Posizione FASTA',c?(c.fastaPos===null?'inserzione':String(c.fastaPos)):'—']];
    box.replaceChildren();fields.forEach(([label,value],index)=>{const item=element('div',undefined,'inspect-item');item.append(element('small',label));const strong=element('strong',value,index===0?'base-large':'');if(index===0)strong.style.color=colors[read.sequence[i]]||'#8b7b89';item.append(strong);box.append(item);});
    if(c&&state.result.orientation==='reverse')box.append(element('span','Nel confronto con la reference questa base è complementata: '+c.readBase+'.','small muted'));
  }
  function jumpToRead(index) {
    if(!state.read)return;state.selected=Math.max(0,Math.min(state.read.sequence.length-1,index));
    state.viewStart=Math.max(0,state.selected-Math.floor(Number($('view-count').value)/2));$('jump-base').value=state.selected+1;
    syncTraceControls();drawTrace();renderInspector();
  }
  function analyze() {
    let r;try{r=region();}catch(error){status(error.message,'error');return;}
    if(!state.read.qualities){status('Phred assenti: non è possibile eseguire questa analisi con qualità.','error');return;}
    invalidate();state.busy=true;const token=state.run;toggleControls();text('analyze','Allineamento in corso…');
    status('Allineamento circolare in corso. Confronto della lettura e del reverse complement…','busy');
    const payload={sequence:state.read.sequence.slice(r.start,r.end),qualities:state.read.qualities.slice(r.start,r.end),reference:state.reference.sequence,
      options:{origin:state.origin,circular:true,orientation:$('orientation').value,trimStart:r.start}};
    const source=$('sanger-core').textContent+'\nself.onmessage=function(e){try { self.postMessage({result:SangerCore.alignRead(e.data.sequence,e.data.qualities,e.data.reference,e.data.options)}); } catch(err){self.postMessage({error:err.message});}};';
    const url=URL.createObjectURL(new Blob([source],{type:'text/javascript'}));const started=performance.now();
    try {
      const worker=new Worker(url);URL.revokeObjectURL(url);state.worker=worker;
      const finish=()=>{worker.terminate();if(state.worker===worker)state.worker=null;state.busy=false;toggleControls();text('analyze','Allinea alla rCRS →');};
      worker.onmessage=event=>{if(token!==state.run)return;finish();if(event.data.error){status(event.data.error,'error');return;}state.elapsed=performance.now()-started;state.result=event.data.result;state.analysisTrim={...r,method:state.trimMethod,cutoff:state.trimMethod==='Mott'?state.trimCutoff:null};state.mapping=new Map(state.result.columns.filter(c=>c.readIndex0!==null).map(c=>[c.readIndex0,c]));state.alignmentPage=0;renderResults();renderStats();drawTrace();renderInspector();const count=state.result.variants.length;status(state.result.summary.accepted?`Allineamento completato: ${count} differenze candidate. Controlla qualità e picchi prima di interpretarle.`:'Allineamento non accettato dalle soglie di controllo. Esamina gli avvisi e il trimming.',state.result.summary.accepted?'':'error');};
      worker.onerror=()=>{if(token!==state.run)return;finish();status('Il browser non ha completato il calcolo. Prova a ricaricare il file o apri il tool da un sito HTTPS.','error');};
      worker.postMessage(payload);
    } catch(error){URL.revokeObjectURL(url);state.busy=false;toggleControls();text('analyze','Allinea alla rCRS →');status('Impossibile avviare il calcolo nel browser: '+error.message,'error');}
  }
  function renderResults() {
    const result=state.result,s=result.summary;$('result-empty').hidden=true;$('result-content').hidden=false;
    text('alignment-state',s.accepted?'Allineato · da revisionare':'Non accettato');$('alignment-state').className='pill '+(s.accepted?'ok':'warn');
    const summary=$('alignment-summary');summary.replaceChildren();
    for(const [label,value] of [['Orientamento',result.orientation==='forward'?'forward':'reverse complement'],['Identità',(100*s.identity).toLocaleString('it-IT',{maximumFractionDigits:2})+'%'],['Basi della lettura',number(s.readLength)],['Score',String(result.score)],['Passaggio origine',s.wrapsRcrs?'16.569 → 1':'no']]){const span=element('span',label+': ');span.append(element('strong',value));summary.append(span);}
    const warnings=$('result-warnings');warnings.replaceChildren();for(const warning of result.warnings)warnings.append(element('p',warning));
    if(result.variants.some(v=>v.quality!==null&&v.quality>=20&&v.quality<30))warnings.append(element('p','Alcune differenze hanno Q fra 20 e 29. Meritano una revisione particolarmente attenta dei picchi.'));
    text('variant-count',result.variants.length);renderVariants();renderAlignment();
  }
  function anchorText(v,coord) {
    if(v.type==='INS')return (v.leftAnchor?.[coord]??'?')+' | '+(v.rightAnchor?.[coord]??'?');
    const positions=coord==='rcrsPos'?v.positions:v.fastaPositions;
    return positions?.length>1?positions[0]+'–'+positions.at(-1):String(v[coord]??'—');
  }
  function variantIndex(v) {
    if(v.readPos!==null)return v.readPos-1;
    const c=state.result.columns.find(c=>c.readPos!==null&&(c.rcrsPos===v.leftAnchor?.rcrsPos||c.rcrsPos===v.rightAnchor?.rcrsPos));return c?c.readPos-1:null;
  }
  function renderVariants() {
    const body=$('variants-body');body.replaceChildren();const variants=state.result.variants;
    $('no-variants').hidden=variants.length!==0;$('variants-table').hidden=variants.length===0;
    for(const v of variants){const tr=element('tr');const type=element('td');type.append(element('span',v.type==='AMBIGUOUS'?'AMB':v.type,'variant-type'+(['INS','DEL'].includes(v.type)?' indel':'')));tr.append(type);
      tr.append(element('td',anchorText(v,'rcrsPos')),element('td',anchorText(v,'fastaPos')),element('td',(v.ref||'−')+' → '+(v.alt||'−')));
      const q=element('td',v.quality===null?(v.flankQuality===null?'—':'fianchi ≥ '+v.flankQuality):String(v.quality),v.lowQuality?'low-quality':'');tr.append(q);
      let note=v.type==='SNP'?(v.quality<20?'Qualità bassa · rivedere':v.quality<30?'Q20–29 · rivedere':'Candidata · verificare il picco'):v.type==='AMBIGUOUS'?'Base ambigua · non risolta':v.placementAmbiguous?'Ripetizione · posizione non univoca':'Indel · rappresentazione non normalizzata';
      if(v.type==='INS')note+='; coordinate fra le due ancore';if(v.type==='DEL')note+='; Q dei fianchi';
      tr.append(element('td',note,'note-cell'));const action=element('td');const index=variantIndex(v);if(index!==null){const button=element('button','Vedi picco','button small-button');button.addEventListener('click',()=>{jumpToRead(index);$('trace-title').scrollIntoView({behavior:'smooth',block:'start'});});action.append(button);}else action.textContent='—';tr.append(action);body.append(tr);
    }
  }
  function renderAlignment() {
    if(!state.result)return;const columns=state.result.columns,size=40;const max=Math.max(0,Math.floor((columns.length-1)/size));state.alignmentPage=Math.max(0,Math.min(max,state.alignmentPage));
    const start=state.alignmentPage*size,end=Math.min(start+size,columns.length),slice=columns.slice(start,end);text('align-range',`Colonne ${start+1}–${end} / ${columns.length}`);$('align-prev').disabled=state.alignmentPage===0;$('align-next').disabled=state.alignmentPage===max;
    const table=element('table',undefined,'alignment-grid');const body=element('tbody');
    for(const row of [{name:'rCRS',key:'rcrsPos',coord:true},{name:'FASTA',key:'fastaPos',coord:true},{name:'Reference',key:'refBase'},{name:'Lettura',key:'readBase'},{name:'Q',key:'quality'},{name:'AB1',key:'readPos',coord:true}]){const tr=element('tr');tr.append(element('th',row.name,'row-name'));for(const c of slice){let val=c[row.key];if(val===null)val=row.key.includes('Pos')&&c.op==='I'?'+':'−';const td=element('td',String(val),(row.coord?'coord ':'')+(['X','I','D','?'].includes(c.op)?'mismatch ':'')+(c.lowQuality?'lowq ':'')+(c.rcrsPos===1?'origin-cell ':'')+(c.readIndex0!==null?'base-click':''));if(row.key==='readBase'||row.key==='refBase')td.style.color=colors[c[row.key]]||'#787471';if(c.readIndex0!==null)td.addEventListener('click',()=>jumpToRead(c.readIndex0));tr.append(td);}body.append(tr);}table.append(body);$('alignment-grid').replaceChildren(table);
  }
  function download(name,content,type='text/plain;charset=utf-8') {
    const url=URL.createObjectURL(new Blob([content],{type}));const a=element('a');a.href=url;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);
  }
  const wrap=seq=>(seq.match(/.{1,80}/g)||[]).join('\n');
  const cleanName=name=>name.replace(/\.[^.]+$/,'').replace(/[^A-Za-z0-9_-]/g,'_').slice(0,100)||'lettura';
  const cell=value=>String(value??'').replace(/[\t\r\n]/g,' ');
  const tsv=rows=>rows.map(row=>row.map(cell).join('\t')).join('\n')+'\n';
  function report() {
    return {tool:'SangerLab',version:Core.VERSION,generatedAt:new Date().toISOString(),input:{fileName:state.read.fileName,sha256:state.readHash,baseCount:state.read.sequence.length,tags:state.read.tags},reference:{name:state.reference.name,accession:'NC_012920.1',sequenceSha256:state.referenceHash,length:N,originRcrs1:state.origin,exactRotationVerified:true},coordinateConvention:'All public positions 1-based; trim.start zero-based inclusive, trim.end zero-based exclusive. Insertions between leftAnchor/rightAnchor; not mitochondrial-normalized.',trim:state.analysisTrim,analysis:state.result,elapsedMs:state.elapsed,scope:'Single original basecaller read; no consensus, no heteroplasmy estimate, no haplogroup assignment.'};
  }
  function exportResults(kind) {
    if(!state.result)return;const name=cleanName(state.read.fileName),r=state.result,s=r.summary,trim=state.analysisTrim;
    const label=`${name} trim_AB1=${trim.start+1}-${trim.end} orientation=${r.orientation} rCRS=${s.rcrsStart}-${s.rcrsEnd} wraps_origin=${s.wrapsRcrs} accepted=${s.accepted}`;
    if(kind==='fasta')download(name+'_trimmed.fasta','>'+label+'\n'+wrap(r.sequence)+'\n');
    if(kind==='alignment')download(name+'_alignment.fasta','>rCRS_aligned rCRS='+s.rcrsStart+'-'+s.rcrsEnd+'\n'+wrap(r.alignment.reference)+'\n>'+label+'\n'+wrap(r.alignment.query)+'\n');
    if(kind==='variants'){
      const rows=[['type','rcrs_position_or_between_anchors','fasta_position_or_between_anchors','ref','alt','read_pos_AB1','basecaller_phred','flank_phred_min','low_quality_Q20','placement_ambiguous','left_anchor_rcrs','right_anchor_rcrs','left_anchor_observed','right_anchor_observed','alignment_accepted','note']];
      for(const v of r.variants)rows.push([v.type,anchorText(v,'rcrsPos'),anchorText(v,'fastaPos'),v.ref,v.alt,v.readPos,v.quality,v.flankQuality,v.lowQuality,v.placementAmbiguous,v.leftAnchor?.rcrsPos,v.rightAnchor?.rcrsPos,v.leftAnchor?.observed,v.rightAnchor?.observed,s.accepted,v.note]);
      download(name+'_differences.tsv',tsv(rows),'text/tab-separated-values;charset=utf-8');
    }
    if(kind==='quality'){
      const rows=[['read_position_AB1','original_base','original_phred','retained','fasta_position','rcrs_position','aligned_base','alignment_op']];
      for(let i=0;i<state.read.sequence.length;i++){const c=state.mapping.get(i);rows.push([i+1,state.read.sequence[i],state.read.qualities?.[i],i>=trim.start&&i<trim.end,c?.fastaPos,c?.rcrsPos,c?.readBase,c?.op]);}
      download(name+'_qualities.tsv',tsv(rows),'text/tab-separated-values;charset=utf-8');
    }
    if(kind==='report')download(name+'_report.json',JSON.stringify(report(),null,2)+'\n','application/json');
  }
  $('ab1-file').addEventListener('change',event=>loadRead(event.target.files[0]));
  $('fasta-file').addEventListener('change',async event=>{const file=event.target.files[0];if(!file)return;try{if(file.size>1000000)throw new Error('FASTA troppo grande per questa versione mtDNA.');setReference(Core.parseFasta(await file.text()));status('Reference verificata base per base. Prima posizione FASTA = rCRS '+state.origin+'. Ricalcola l’allineamento.');}catch(error){status(error.message,'error');}});
  $('restore-reference').addEventListener('click',()=>{setReference({name:bundled.name,sequence:bundled.sequence});status('Reference rCRS_15990 ripristinata.');});
  $('auto-trim').addEventListener('click',()=>trimAuto());$('cutoff').addEventListener('change',()=>trimAuto());
  for(const id of ['trim-start','trim-end'])$(id).addEventListener('input',()=>{invalidate();state.autoEmpty=false;state.trimMethod='manuale';updateTrimInfo();status('Intervallo modificato. Premi «Allinea alla rCRS» per aggiornare i risultati.');});
  $('orientation').addEventListener('change',()=>{invalidate();status('Orientamento modificato. Ricalcola l’allineamento.');});
  $('reset-trim').addEventListener('click',()=>{invalidate();state.autoEmpty=false;state.trimMethod='nessuno';$('trim-start').value=1;$('trim-end').value=state.read.sequence.length;updateTrimInfo();status('Tutte le basi sono conservate. Controlla le estremità prima di allineare.');});
  $('analyze').addEventListener('click',analyze);
  $('view-count').addEventListener('change',()=>{syncTraceControls();drawTrace();});
  $('trace-slider').addEventListener('input',event=>{state.viewStart=Number(event.target.value);drawTrace();});
  $('prev-trace').addEventListener('click',()=>{state.viewStart-=Math.ceil(Number($('view-count').value)*.7);drawTrace();});
  $('next-trace').addEventListener('click',()=>{state.viewStart+=Math.ceil(Number($('view-count').value)*.7);drawTrace();});
  const jump=()=>{const n=Number($('jump-base').value);if(Number.isInteger(n)&&n>=1&&n<=state.read.sequence.length)jumpToRead(n-1);};
  $('jump-button').addEventListener('click',jump);$('jump-base').addEventListener('keydown',e=>{if(e.key==='Enter')jump();});
  $('trace-canvas').addEventListener('pointerdown',event=>{if(!plotGeometry)return;const rect=$('trace-canvas').getBoundingClientRect(),x=event.clientX-rect.left,g=plotGeometry;const sample=g.sampleLo+(x-g.left)/g.width*(g.sampleHi-g.sampleLo);let best=g.lo;for(let i=g.lo+1;i<g.hi;i++)if(Math.abs(state.read.peakPositions[i]-sample)<Math.abs(state.read.peakPositions[best]-sample))best=i;state.selected=best;$('jump-base').value=best+1;drawTrace();renderInspector();});
  for(const [tab,other,panel,otherPanel] of [['tab-variants','tab-alignment','variants-panel','alignment-panel'],['tab-alignment','tab-variants','alignment-panel','variants-panel']])$(tab).addEventListener('click',()=>{$(tab).className='tab active';$(tab).setAttribute('aria-selected','true');$(other).className='tab';$(other).setAttribute('aria-selected','false');$(panel).hidden=false;$(otherPanel).hidden=true;});
  $('align-prev').addEventListener('click',()=>{state.alignmentPage--;renderAlignment();});$('align-next').addEventListener('click',()=>{state.alignmentPage++;renderAlignment();});
  for(const kind of ['fasta','alignment','variants','quality','report'])$('export-'+kind).addEventListener('click',()=>exportResults(kind));
  const drop=$('read-drop');for(const name of ['dragenter','dragover'])drop.addEventListener(name,e=>{e.preventDefault();drop.classList.add('dragging');});
  for(const name of ['dragleave','drop'])drop.addEventListener(name,e=>{e.preventDefault();drop.classList.remove('dragging');if(name==='drop')loadRead(e.dataTransfer.files[0]);});
  new ResizeObserver(()=>drawTrace()).observe($('trace-canvas'));
  setReference({name:bundled.name,sequence:bundled.sequence});
  // Read-only snapshot for reproducible browser QA; no internal ABI metadata exposed.
  globalThis.SangerLab={getReport:()=>state.result?report():null};
})();
