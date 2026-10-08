'use strict';
// Experimental retrospective validation. Input is loaded before inference;
// the answer file is not opened until the prediction has been frozen.
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
require('../dist/sim-core.js');require('../dist/mechanism-engine.js');
const C=globalThis.PrimingCore,E=globalThis.MechanismEngine,dir=path.join(__dirname,'mechanism_validation');
const files={extca:'jp286282_2mm_ca',iono:'jp286282_ionomycin',pdbu:'jp286282_pdbu',hk:'sciadv_hk'};
const id=process.argv[2]||'extca',stem=files[id];assert.ok(stem,'Unknown validation case');
const input=JSON.parse(fs.readFileSync(path.join(dir,stem+'_input.json'),'utf8'));
const reference=C.resolvePreset(input.referencePreset,'publication');
const keys=[...E.FAMILIES.priming,...E.FAMILIES.calcium.filter(k=>k!=='caRestReference')];
const mode=process.env.SEARCH_MODE||'fast',maxChanged=Number(process.env.MAX_CHANGED||2);
const request={referenceName:input.referencePreset,reference,locks:input.locks,keys,bounds:Object.fromEntries(keys.map(k=>[k,E.defaultBounds(reference,k)])),targets:input.targets,protocol:input.protocol,settings:{mode,maxChanged,lambdaChange:.005,lambdaCount:.01,nearTolerance:.1,minimalCriterion:'near',seed:2026}};
assert.ok(!('published' in request));assert.ok(!('expectedDirections' in request));
const started=Date.now(),result=E.infer(request);
const frozen={case:id,frozenAt:new Date().toISOString(),appVersion:E.VERSION,request,result};
const predictionPath=process.env.REPORT_PATH||path.join('/tmp','priminglab-'+id+'-prediction.json');
fs.writeFileSync(predictionPath,JSON.stringify(frozen));
// The prediction now exists on disk. Only now may the published answer be read.
const answer=JSON.parse(fs.readFileSync(path.join(dir,stem+'_answer.json'),'utf8'));
const assessed=Object.entries(answer.expectedDirections).map(([key,direction])=>{const fold=result.minimal.params[key]/reference[key],inferred=fold>1.01?'up':fold<.99?'down':'same';return{key,direction,inferred,fold,directionCorrect:inferred===direction}});
const ranked=[...result.familyModels].filter(x=>x.bestError!==null).sort((a,b)=>b.relativeSupport-a.relativeSupport);
const summary={kind:input.kind,case:id,mode,maxChanged,elapsedSeconds:Math.round((Date.now()-started)/1000),predictionPath,changed:result.minimal.changed,coverage:result.minimal.coverage,targetsMet:result.minimal.rows.filter(x=>x.met).length+'/'+result.minimal.rows.length,expectedFamily:answer.expectedFamily,topFamily:ranked[0]?.name||null,assessed};
console.log(JSON.stringify(summary));
if(process.env.ASSERT_SCIENCE==='1'){
 assert.equal(summary.topFamily,answer.expectedFamily,'Expected mechanism family was not top-ranked');
 assert.ok(assessed.every(x=>x.directionCorrect),'Expected key mechanism directions not recovered');
}
