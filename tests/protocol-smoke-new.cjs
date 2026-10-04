"use strict";
require("../dist/sim-core.js");const C=globalThis.PrimingCore,p=C.resolvePreset("jp-control","publication");
const train=C.simulateTrain(200,40,p),recovery=C.simulateRecovery(train,[.125,.25,.5,1,2,4,8,16],p,"last5"),conditional=C.conditioningEstimator(p,{preFreq:10,prePulses:4,testFreq:200,testPulses:40,window:10});
if(recovery.points.length!==8||!recovery.points.every(x=>Number.isFinite(x.release)&&Number.isFinite(x.pFusion)&&Number.isFinite(x.mass)))process.exit(1);
if(conditional.test.length!==40||!Number.isFinite(conditional.pFusionEstimate)||conditional.timing[4]!==.5)process.exit(1);
console.log(JSON.stringify({recoveryPoints:recovery.points.length,conditionedPoints:conditional.test.length,pFusionEstimate:conditional.pFusionEstimate}));
