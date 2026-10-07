'use strict';
// Blind synthetic validation: only this fixture knows the published HK rates.
// Source: Science Advances aea0449 Supplementary Table S2, PDF p.19.
const fs=require('node:fs');require('../dist/sim-core.js');require('../dist/mechanism-engine.js');
const C=PrimingCore,E=MechanismEngine,wt=C.resolvePreset('sciadv','publication');
const publication={...wt,nSites:3150,pRel0:.379,kf1:.419,kb1:.194,slope1:221000,kf2:.2916,kb2:.3514,slope2:1810000};
const protocol={pulses:40,pprFrequency:200,recoveryFrequency:200,frequencies:[5,20,100,200],delays:[.02,.05,.1,.2,.3,.5,1,2,5,10]};
const ids=['p1','ppr','ss:5','ss:20','ss:100','ss:200','r:0.1','r:0.3','r:1'];
const target=E.extractPhenotype(publication,protocol,ids);
const allowed=[...E.DEFAULT_KEYS,'nSites'];
const request={referenceName:'Sci Adv Table S2 WT/-',reference:wt,locks:{pRel0:.379},keys:allowed,bounds:Object.fromEntries(allowed.map(k=>[k,E.defaultBounds(wt,k)])),protocol,targets:[...ids.map(id=>({id,kind:'quantitative',mean:target.values[id],unit:'absolute',uncertainty:'none',weight:1})),{id:'frp',kind:'quantitative',mean:2506,unit:'absolute',uncertainty:'scale',error:172,weight:3}],settings:{mode:process.env.SEARCH_MODE||'standard',maxChanged:3,lambdaChange:.005,lambdaCount:.01,nearTolerance:.1,seed:2026}};
let stage;const result=E.infer(request,p=>{if(p.stage!==stage){stage=p.stage;console.log('Stage',stage,'evaluations',p.evaluations)}});
const summary={kind:'Hybrid validation: simulated Table S2 phenotype plus experimental Table S1 FRP; not an independent experimental fit',source:'https://www.science.org/doi/10.1126/sciadv.aea0449',table:'Supplementary Tables S1 and S2',publication,coverage:result.minimal.coverage,numericalCoverage:result.numerical.coverage,changed:result.minimal.changed,comparison:[...E.DEFAULT_KEYS,'nSites'].map(key=>({key,wt:wt[key],publication:publication[key],inferred:result.minimal.params[key]})),targets:result.minimal.rows,identifiability:result.identifiability,evaluations:result.evaluations};
console.log(JSON.stringify(summary,null,2));
if(process.env.REPORT_PATH){fs.writeFileSync(process.env.REPORT_PATH,JSON.stringify({summary,result},null,2))}
