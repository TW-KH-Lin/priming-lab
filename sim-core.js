"use strict";

/* Pure simulation engine. No DOM dependencies. */
(function(global){
  const clamp=(v,a,b)=>Math.min(b,Math.max(a,v));
  const clone=o=>JSON.parse(JSON.stringify(o));

  const BASE={
    frequencies:[5,20,50,100,200],pulses:40,nSites:2639,pRel0:.393,
    kf1:.402517,kb1:.18472,slope1:1.23939e7,kf2:.207254,kb2:.248014,slope2:1.27727e7,
    tauTsl:.09,fractionTsl:.16,kHalf1:2.8e-7,kRefract:5000,useMM:true,model:2,
    caRestActual:5e-8,caRestReference:5e-8,useMultiCa:false,caTau:.06,caAmpl:1.1e-7,
    caTauFast:.06,caTauSlow:.23,caFracSlow:.15,caAmplGlobal:4.54e-7,caTauLocal:.00025,caAmplLocal:3e-5,
    caApScaleWithY:true,yInc:.39,zDec:.4,yMax:1.32,zMin:.75,tauY:.014,tauZ:3,yPower:4.5,
    pFusionMode:"dynamic",odeStep:.001,source:"10.1073/pnas.2207987119"
  };

  const PRESETS={
    pnas2022:{label:"PNAS 2022 · Figure 4",...BASE},
    pnas2024:{label:"PNAS 2024 · PC-FSIN mean",...BASE,nSites:25.5,pRel0:.6,kf1:.61,kb1:.38,slope1:1.06430155210643e7,kf2:.24,kb2:.3,slope2:2.2172949002217297e7,tauTsl:.0725,fractionTsl:.18,kHalf1:3.3e-7,caTau:.1025,source:"10.1073/pnas.2322550121"},
    "jp-control":{label:"JP286282 · ionomycin control",...BASE,nSites:2622,pRel0:.22,kf1:.37,kb1:.221,slope1:2.245e6,kf2:.199,kb2:.253,slope2:1.014e6,fractionTsl:.1,kRefract:2.6,useMM:false,model:3,useMultiCa:true,caApScaleWithY:false,yMax:1.31,zMin:.87,tauY:.017,yPower:4,odeStep:.0001,source:"10.1113/JP286282"},
    "jp-iono":{label:"JP286282 · +2.5 µM ionomycin",...BASE,nSites:2622,pRel0:.3,kf1:.37,kb1:.221,slope1:2.245e6,kf2:.199,kb2:.253,slope2:1.014e6,fractionTsl:.15,kRefract:4.2,useMM:false,model:3,useMultiCa:true,caRestActual:1.3e-7,caRestReference:5e-8,caApScaleWithY:false,yMax:1.16,zMin:1,tauY:.017,yPower:4,odeStep:.0001,source:"10.1113/JP286282"},
    "jp-extca":{label:"JP286282 · 2.0 mM external Ca",...BASE,nSites:2910,pRel0:.42,kf1:.374015,kb1:.215,slope1:2.181e6,kf2:.191,kb2:.264,slope2:1.042e6,fractionTsl:.1,kRefract:2.9,useMM:false,model:3,useMultiCa:true,caRestActual:9e-8,caRestReference:5e-8,caAmplGlobal:5.45e-7,caAmplLocal:3.6e-5,caApScaleWithY:false,yMax:1.28,zMin:.87,tauY:.017,yPower:4,odeStep:.0001,source:"10.1113/JP286282"},
    "jp-pdbu":{label:"J Physiol 2025 · Rat calyx · 1 µM PDBu",...BASE,nSites:3297,pRel0:.29,kf1:.373,kb1:.210,slope1:2.051e6,kf2:.378,kb2:.271,slope2:2.436e6,fractionTsl:.11,kRefract:2.8,useMM:false,model:3,useMultiCa:true,caRestActual:5e-8,caRestReference:5e-8,caApScaleWithY:false,yMax:1.27,zMin:.87,tauY:.017,yPower:4,odeStep:.0001,source:"10.1113/JP286282"},
    sciadv:{label:"Science Advances 2026 · WT/−",...BASE,nSites:3300,pRel0:.341,kf1:.5,kb1:.231,slope1:8.59e5,kf2:.1357,kb2:.3514,slope2:1.22e6,tauTsl:.07,fractionTsl:.08,kRefract:3.3,useMM:false,model:3,useMultiCa:true,caTauFast:.05,caTauSlow:.24,caFracSlow:.12,caAmplGlobal:3.9e-7,caTauLocal:.00036,caAmplLocal:3.3e-5,caApScaleWithY:false,yInc:.36,yMax:1.22,zMin:.77,tauY:.014,yPower:4.5,odeStep:.0001,source:"10.1126/sciadv.aea0449"}
  };

  function resolvePreset(key,target="publication"){
    const p=clone(PRESETS[key]||PRESETS.pnas2022);p.key=key;p.target=target;
    if(target==="python"){
      p.caRestReference=p.caRestActual;
      if(key==="pnas2024"){p.pFusionMode="dynamic";p.caApScaleWithY=true}
    }else{
      if(key==="pnas2024"){p.pFusionMode="constant";p.caApScaleWithY=false}
      if(key.startsWith("jp-")||key==="sciadv")p.caRestReference=5e-8;
    }
    return p;
  }

  function effectiveCa(s,p){return p.caRestActual+s.caFast+(p.useMultiCa?s.caSlow+s.caLocal:0)}
  function rates(ca,p){
    const above=Math.max(ca-p.caRestReference,0);let k1=p.kf1+p.slope1*above;
    if(p.useMM)k1/=1+above/Math.max(p.kHalf1,1e-18);
    return [Math.max(k1,0),Math.max(p.kf2+p.slope2*above,0)];
  }
  function steadyState(p){
    const [k1,k2]=rates(p.caRestActual,p),den=p.kb1*p.kb2+k1*(p.kb2+k2);
    const ts=p.nSites*k1*k2/den,ls=p.nSites*p.kb2*k1/den;
    return {es:p.nSites-ts-ls,ls,ts,tsl:0,ers:0,caFast:0,caSlow:0,caLocal:0};
  }
  function deriv(s,p){
    const [k1,k2]=rates(effectiveCa(s,p),p),b3=1/Math.max(p.tauTsl,1e-12);let d;
    if(p.model===1)d={ts:k2*s.ls-p.kb2*s.ts,es:-k1*s.es+p.kb1*s.ls,ls:k1*s.es+p.kb2*s.ts-(p.kb1+k2)*s.ls,tsl:0,ers:0};
    else if(p.model===2)d={ts:k2*s.ls-p.kb2*s.ts,es:-k1*s.es+p.kb1*s.ls,ls:k1*s.es+p.kb2*s.ts-(p.kb1+k2)*s.ls+b3*s.tsl,tsl:-b3*s.tsl,ers:0};
    else d={ts:k2*s.ls-p.kb2*s.ts,es:-k1*s.es+p.kb1*s.ls+p.kRefract*s.ers,ls:k1*s.es+p.kb2*s.ts-(p.kb1+k2)*s.ls+b3*s.tsl,tsl:-b3*s.tsl,ers:-p.kRefract*s.ers};
    d.caFast=-s.caFast/Math.max(p.useMultiCa?p.caTauFast:p.caTau,1e-12);
    d.caSlow=p.useMultiCa?-s.caSlow/Math.max(p.caTauSlow,1e-12):0;
    d.caLocal=p.useMultiCa?-s.caLocal/Math.max(p.caTauLocal,1e-12):0;return d;
  }
  const STATE_KEYS=["es","ls","ts","tsl","ers","caFast","caSlow","caLocal"];
  function add(a,b,h){const o={};for(const k of STATE_KEYS)o[k]=a[k]+h*b[k];return o}
  function normalize(s,p){
    for(const k of ["es","ls","ts","tsl","ers"])s[k]=Math.max(0,s[k]);
    const total=s.es+s.ls+s.ts+s.tsl+s.ers,scale=p.nSites/Math.max(total,1e-18);
    for(const k of ["es","ls","ts","tsl","ers"])s[k]*=scale;return s;
  }
  function rk4(s,h,p){
    const k1=deriv(s,p),k2=deriv(add(s,k1,h/2),p),k3=deriv(add(s,k2,h/2),p),k4=deriv(add(s,k3,h),p),o={};
    for(const k of STATE_KEYS)o[k]=s[k]+h*(k1[k]+2*k2[k]+2*k3[k]+k4[k])/6;return normalize(o,p);
  }
  function integrate(s,duration,p){
    if(!Number.isFinite(duration)||duration<0)throw new RangeError("Integration duration must be finite and non-negative");
    if(duration===0)return clone(s);let hTarget=Math.max(p.odeStep||.001,1e-9);
    if(p.model===3)hTarget=Math.min(hTarget,.2/Math.max(p.kRefract,1e-12));
    if(p.useMultiCa){if(!Number.isFinite(p.caTauLocal)||p.caTauLocal<=0)throw new RangeError("Local calcium time constant must be positive");hTarget=Math.min(hTarget,.25*p.caTauLocal)}
    const n=Math.max(1,Math.ceil(duration/hTarget)),h=duration/n;let out=clone(s);for(let i=0;i<n;i++)out=rk4(out,h,p);return out;
  }
  function fusion(y,z,p){return clamp(p.pRel0*(p.pFusionMode==="constant"?1:Math.pow(y,p.yPower)*z),0,1)}
  function applyAP(s,y,z,p){
    const pf=fusion(y,z,p),ts=s.ts,tsl=s.tsl,ls=s.ls;let release;
    if(p.model===1){release=ts*pf;s.ts=ts*(1-pf);s.es+=release}
    else{release=(ts+tsl)*pf;s.ts=ts*(1-pf);s.ls=ls*(1-p.fractionTsl);s.tsl=tsl*(1-pf)+ls*p.fractionTsl;if(p.model===2)s.es+=release;else s.ers+=release}
    const scale=p.caApScaleWithY?y:1;
    if(p.useMultiCa){s.caFast+=scale*p.caAmplGlobal*(1-p.caFracSlow);s.caSlow+=scale*p.caAmplGlobal*p.caFracSlow;s.caLocal+=scale*p.caAmplLocal}else s.caFast+=scale*p.caAmpl;
    normalize(s,p);return {state:s,release,pFusion:pf,yPost:y+p.yInc*(p.yMax-y),zPost:z-p.zDec*(z-p.zMin)};
  }
  function buildTrain(freq,n){return Array.from({length:n},(_,i)=>(i+1)/freq)}
  function buildConditioned(preFreq,preN,testFreq,testN){
    let t=0,out=[];for(let i=0;i<preN;i++){t+=1/preFreq;out.push(t)}t+=1/preFreq;out.push(t);for(let i=1;i<testN;i++){t+=1/testFreq;out.push(t)}return out;
  }
  function simulateTiming(times,p){
    let s=steadyState(p),y=1,z=1,prev=0;const rows=[];
    for(let i=0;i<times.length;i++){
      const dt=times[i]-prev;s=integrate(s,dt,p);if(i>0){y=1+(y-1)*Math.exp(-dt/p.tauY);z=1+(z-1)*Math.exp(-dt/p.tauZ)}
      const pre=clone(s),ca=effectiveCa(s,p),ap=applyAP(s,y,z,p);s=ap.state;
      rows.push({pulse:i+1,time:times[i],release:ap.release,pFusion:ap.pFusion,y,z,ca,caFast:pre.caFast,caSlow:pre.caSlow,caLocal:pre.caLocal,ES:pre.es,LS:pre.ls,TS:pre.ts,TSL:pre.tsl,ERS:pre.ers,mass:pre.es+pre.ls+pre.ts+pre.tsl+pre.ers});
      y=ap.yPost;z=ap.zPost;prev=times[i];
    }
    const m1=Math.max(rows[0].release,1e-18);rows.forEach(r=>r.releaseNorm=r.release/m1);
    return {rows,finalState:clone(s),yPost:y,zPost:z,endTime:prev,params:p};
  }
  function simulateTrain(freq,n,p){if(!Number.isFinite(freq)||freq<=0)throw new RangeError("Frequency must be positive");if(!Number.isInteger(n)||n<1)throw new RangeError("Pulse count must be a positive integer");const out=simulateTiming(buildTrain(freq,n),p);out.frequency=freq;return out}

  function simulateRecovery(train,intervals,p,baselineMode="last5"){
    const sorted=[...new Set(intervals.filter(v=>Number.isFinite(v)&&v>0))].sort((a,b)=>a-b),last=train.rows.at(-1).release;
    const last5=train.rows.slice(-Math.min(5,train.rows.length)).reduce((a,b)=>a+b.release,0)/Math.min(5,train.rows.length),baseline=baselineMode==="last"?last:last5,m1=train.rows[0].release;
    const points=sorted.map(interval=>{
      const s=integrate(train.finalState,interval,p),y=1+(train.yPost-1)*Math.exp(-interval/p.tauY),z=1+(train.zPost-1)*Math.exp(-interval/p.tauZ),pf=fusion(y,z,p),release=(s.ts+s.tsl)*pf;
      const recoveryDen=m1-baseline;
      return {interval,release,overM1:Math.abs(m1)>1e-18?release/m1:NaN,overLast:Math.abs(last)>1e-18?release/last:NaN,recovered:Math.abs(recoveryDen)>1e-18?(release-baseline)/recoveryDen:NaN,pFusion:pf,ES:s.es,LS:s.ls,TS:s.ts,TSL:s.tsl,ERS:s.ers,ca:effectiveCa(s,p),caFast:s.caFast,caSlow:s.caSlow,caLocal:s.caLocal,y,z,mass:s.es+s.ls+s.ts+s.tsl+s.ers};
    });
    return {points,m1,last,last5,baselineMode,train};
  }
  function linearFit(x,y,n=x.length){
    const count=Math.min(n,x.length,y.length),xs=x.slice(0,count),ys=y.slice(0,count);if(count<2)return null;
    const mx=xs.reduce((a,b)=>a+b,0)/count,my=ys.reduce((a,b)=>a+b,0)/count;let num=0,den=0;
    for(let i=0;i<count;i++){num+=(xs[i]-mx)*(ys[i]-my);den+=(xs[i]-mx)**2}if(den<1e-18)return null;
    const slope=num/den;return {slope,intercept:my-slope*mx,count};
  }
  function conditioningEstimator(p,{preFreq=10,prePulses=4,testFreq=200,testPulses=40,window=10}={}){
    const cond=simulateTiming(buildConditioned(preFreq,prePulses,testFreq,testPulses),p),ctrl=simulateTrain(testFreq,testPulses,p),test=cond.rows.slice(-testPulses);
    const difference=test.map((r,i)=>r.releaseNorm-ctrl.rows[i].releaseNorm),cumulative=difference.map((_,i)=>i?difference.slice(0,i).reduce((a,b)=>a+b,0):0),fit=linearFit(cumulative,difference,window);
    return {conditioned:cond,control:ctrl,test,difference,cumulative,fit,pFusionEstimate:fit?-fit.slope:NaN,actualPFusion:test.map(r=>r.pFusion),timing:cond.rows.map(r=>r.time),prePulses,testPulses};
  }
  function eq31(p,{frequency=10,pulses=40}={}){
    const train=simulateTrain(frequency,pulses,p),ppr=train.rows[1].release/train.rows[0].release,dm=train.rows.at(-1).release/train.rows[0].release,den=1-dm,estimate=Math.abs(den)>1e-18?(1-ppr)/den:NaN;
    return {train,ppr,dm,estimate,modelPFusion:train.rows[0].pFusion};
  }
  function invariants(result,p){
    const rows=result.rows||[],massError=Math.max(0,...rows.map(r=>Math.abs(r.mass-p.nSites))),minPool=Math.min(Infinity,...rows.flatMap(r=>[r.ES,r.LS,r.TS,r.TSL,r.ERS])),pfOk=rows.every(r=>r.pFusion>=0&&r.pFusion<=1);
    return {massError,minPool,pFusionBounded:pfOk,model1Branches:p.model!==1||rows.every(r=>Math.abs(r.TSL)<1e-10&&Math.abs(r.ERS)<1e-10),model3ERS:p.model!==3||rows.every(r=>Math.abs(r.ERS)<1e-10)};
  }

  global.PrimingCore={PRESETS,resolvePreset,steadyState,simulateTrain,simulateTiming,simulateRecovery,conditioningEstimator,eq31,invariants,buildConditioned,effectiveCa,integrate};
})(typeof window!=="undefined"?window:globalThis);
