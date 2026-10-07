'use strict';
importScripts('sim-core.js','mechanism-engine.js');
self.onmessage=({data})=>{try{const result=MechanismEngine.infer(data,p=>self.postMessage({type:'progress',...p}));self.postMessage({type:'result',result})}catch(e){self.postMessage({type:'error',message:e.message})}};
