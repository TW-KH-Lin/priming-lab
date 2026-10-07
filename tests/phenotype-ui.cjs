/* Run with jsdom available via NODE_PATH; no browser or network access. */
"use strict";
const assert=require("node:assert/strict"), fs=require("node:fs"), path=require("node:path"), {JSDOM}=require("jsdom");
const dist=path.join(__dirname,"../dist"),dom=new JSDOM(fs.readFileSync(path.join(dist,"index.html"),"utf8"),{runScripts:"outside-only",url:"http://localhost/"}),w=dom.window;
Object.defineProperty(w.SVGElement.prototype,"viewBox",{get(){const [x=0,y=0,width=640,height=330]=(this.getAttribute("viewBox")||"").split(/\s+/).map(Number);return {baseVal:{x,y,width,height}}}});
let workerCount=0;
w.Worker=class {
  constructor(){workerCount++;this.cancelled=false}
  terminate(){this.cancelled=true}
  postMessage({base,request}){setTimeout(()=>{if(this.cancelled)return;try{const result=w.PhenotypeEngine.compare(base,request);this.onmessage({data:{type:"result",result}})}catch(error){this.onmessage({data:{type:"error",message:error.message}})}});}
};
for(const file of ["sim-core.js","phenotype-engine.js","phenotype-ui.js","app.js"])w.eval(fs.readFileSync(path.join(dist,file),"utf8")+(file==="app.js"?";globalThis.__readParams=()=>preset();globalThis.__readState=()=>state;":""));
const $=s=>w.document.querySelector(s),param=()=>JSON.parse(JSON.stringify(w.__readParams()));
function input(selector,value){$(selector).value=String(value);$(selector).dispatchEvent(new w.Event("input",{bubbles:true}));}
const settled=()=>new Promise(resolve=>setTimeout(resolve,30));
(async()=>{
  assert.match($("#status").textContent,/Complete/);
  $("#preset").value="sciadv";$("#preset").dispatchEvent(new w.Event("change"));
  const baseline=param();
  $("#togglePhenotype").click();
  assert.equal($(".phenotype-panel").hidden,false);
  assert.equal($("#guideResults").textContent,"");
  assert.equal(w.document.querySelectorAll("#guideRun").length,1);
  input("#guidePFusion",1.1);input("#guidePool",1);
  input('[data-observation="first"]',"ratio");input('[data-ratio="first"]',1.5);
  input('[data-observation="steady"]',"down");input('[data-observation="recovery"]',"down");
  $("#guideRun").click();await settled();
  assert.match($("#guideStatus").textContent,/observations met/);
  const first=param();assert.equal(first.pRel0,baseline.pRel0*1.1);assert.equal(first.nSites,baseline.nSites);
  assert.ok(w.__readState().phenotypeComparison.tested > 0);
  $("#runButton").click();assert.ok(w.__readState().phenotypeComparison.tested > 0,"Exports retain the comparison after rerunning plots");
  $("#guideRun").click();await settled();assert.deepEqual(param(),first,"Repeat must not accumulate changes");
  input('[data-ratio="first"]',0);$("#guideRun").click();await settled();
  assert.deepEqual(param(),first,"Invalid input must not alter simulation");
  assert.match($("#guideStatus").textContent,/Check/);
  $("#guideUndo").click();assert.deepEqual(param(),baseline,"Restore must recover exact baseline");
  input('[data-ratio="first"]',1.5);$("#guideRun").click();$("#guideRun").click();await settled();
  assert.deepEqual(param(),baseline,"Cancel must leave baseline unchanged");
  $("#preset").value="jp-control";$("#preset").dispatchEvent(new w.Event("change"));
  assert.match($("#guideBaseline").textContent,/J Physiol/);
  assert.equal($("#guideResults").textContent,"");
  $("#togglePhenotype").click();assert.equal($(".phenotype-panel").hidden,true);
  assert.match($("#status").textContent,/Complete/);
  console.log(`Phenotype DOM checks passed: ${workerCount} worker requests; repeat, fixed inputs, cancellation, restoration and preset switching.`);
  w.close();
})().catch(error=>{console.error(error);process.exitCode=1;w.close()});
