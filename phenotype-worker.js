"use strict";
importScripts("sim-core.js", "phenotype-engine.js");
self.onmessage = ({data}) => {
  try {
    const result = PhenotypeEngine.compare(data.base, data.request, (done,total) => self.postMessage({type:"progress",done,total}));
    self.postMessage({type:"result",result});
  } catch (error) { self.postMessage({type:"error",message:error.message}); }
};
