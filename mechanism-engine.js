(function(g){
  'use strict';
  const C=g.PrimingCore, VERSION='49.0.0';
  const DEFAULT_KEYS=['kf1','kb1','slope1','kf2','kb2','slope2'];
  const FAMILIES={priming:DEFAULT_KEYS,calcium:['caRestActual','caRestReference','caAmplGlobal','caAmplLocal','caTauFast','caTauSlow','caFracSlow'],stp:['yInc','yMax','zDec','zMin','tauY','tauZ','yPower'],refractory:['fractionTsl','tauTsl','kRefract']};
  const NAMES={nSites:'Total release-site capacity N',pRel0:'Initial pFusion',kf1:'Basal k₁',kb1:'b₁',slope1:'Ca slope σ₁',kf2:'Basal k₂',kb2:'b₂',slope2:'Ca slope σ₂',caRestActual:'Actual resting Ca²⁺',caRestReference:'Rate-reference Ca²⁺',caAmplGlobal:'AP global Ca amplitude',caAmplLocal:'AP local Ca amplitude',caTauFast:'Fast Ca decay',caTauSlow:'Slow Ca decay',caFracSlow:'Slow Ca fraction',yInc:'Facilitation increment',yMax:'Maximum facilitation',zDec:'Depression decrement',zMin:'Minimum depression',tauY:'Facilitation decay',tauZ:'Depression decay',yPower:'Facilitation exponent',fractionTsl:'TSL fraction',tauTsl:'TSL lifetime',kRefract:'ERS recovery b₄'};
  const QUAL={strongDown:[0,.4],down:[.4,.7],mildDown:[.7,.9],same:[.9,1.1],mildUp:[1.1,1.3],up:[1.3,1.6],strongUp:[1.6,Infinity]};
  const SIM_NUMERIC=['nSites','pRel0','kf1','kb1','slope1','kf2','kb2','slope2','tauTsl','fractionTsl','kHalf1','kRefract','caRestActual','caRestReference','caTau','caAmpl','caTauFast','caTauSlow','caFracSlow','caAmplGlobal','caTauLocal','caAmplLocal','yInc','zDec','yMax','zMin','tauY','tauZ','yPower','odeStep'];
  const clone=x=>JSON.parse(JSON.stringify(x)), clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
  function subsets(keys,size){const out=[];function walk(i,a){if(a.length===size){out.push(a);return}for(let j=i;j<keys.length;j++)walk(j+1,[...a,keys[j]])}walk(0,[]);return out}
  function defaultBounds(p,key){if(key==='nSites')return[p.nSites*.5,p.nSites*2];if(key==='pRel0')return[Math.max(1e-6,p.pRel0*.2),Math.min(1,p.pRel0*2)];if(key==='caRestActual')return[1e-8,1e-6];if(key==='caRestReference')return[1e-8,1e-6];if(key==='caAmplGlobal'||key==='caAmplLocal')return[p[key]*.2,p[key]*5];if(key==='fractionTsl'||key==='caFracSlow'||key==='yInc'||key==='zDec'||key==='zMin')return[0,1];if(key==='yMax')return[1,Math.max(2,p.yMax*2)];if(key==='yPower')return[.2,p.yPower*5];if(key.startsWith('slope'))return[0,p[key]*10];if(key.startsWith('tau')||key.startsWith('caTau'))return[p[key]*.2,p[key]*5];return[p[key]*.1,p[key]*10]}
  function safeParameters(p){
    const positive=['nSites','kf1','kb1','kf2','kb2','tauTsl','kRefract','caTau','caTauFast','caTauSlow','caTauLocal','tauY','tauZ','odeStep','kHalf1','caRestActual','caRestReference'];
    if(positive.some(k=>!(p[k]>0))||['slope1','slope2','caAmpl','caAmplGlobal','caAmplLocal','caRestActual','caRestReference'].some(k=>!(p[k]>=0)))throw Error('Rates, calcium and time constants must be physically valid.');
    if(['pRel0','fractionTsl','caFracSlow','yInc','zDec','zMin'].some(k=>!(p[k]>=0&&p[k]<=1))||p.yMax<1||p.yPower<0||Math.min(p.tauTsl,p.caTau,p.caTauFast,p.caTauSlow,p.caTauLocal,p.tauY,p.tauZ)<1e-6)throw Error('Invalid probability or time constant.');
    // Conservative stability guard for the unchanged simulator's RK4 step.
    const accumulation=(a,t)=>a/(-Math.expm1(-1/333/t));
    const ca=(p.useMultiCa?accumulation(p.caAmplGlobal*(1-p.caFracSlow),p.caTauFast)+accumulation(p.caAmplGlobal*p.caFracSlow,p.caTauSlow)+accumulation(p.caAmplLocal,p.caTauLocal):accumulation(p.caAmpl,p.caTau))*(p.caApScaleWithY?p.yMax:1)+Math.max(0,p.caRestActual-p.caRestReference);
    const h=Math.min(p.odeStep,p.model===3?.2/p.kRefract:Infinity,p.useMultiCa?.25*p.caTauLocal:Infinity);
    if(h*(p.kf1+p.kb1+p.kf2+p.kb2+(p.slope1+p.slope2)*ca+(p.model>1?1/p.tauTsl:0)+(p.model===3?p.kRefract:0))>1.5)throw Error('Parameter bounds are too stiff for the current simulator step; narrow the ranges.');
  }
  function computeSteadyStateDepth(train){
    const rows=train?.rows||[],last=rows.slice(-Math.min(5,rows.length));
    return last.length?last.reduce((sum,row)=>sum+row.releaseNorm,0)/last.length:NaN;
  }
  function trainPeak(rows){
    let index=0;for(let i=1;i<rows.length;i++)if(rows[i].releaseNorm>rows[index].releaseNorm)index=i;
    return {index,value:rows[index]?.releaseNorm??NaN};
  }
  function computeDepressionT50(train,frequencyHz=train?.frequency){
    const rows=train?.rows||[],ssOverP1=computeSteadyStateDepth(train),peak=trainPeak(rows),amplitude=peak.value-ssOverP1;
    if(rows.length<2||!(frequencyHz>0)||!(amplitude>1e-9))return {t50Dep:null,n50Dep:null,reached:false,peakPulse:rows[peak.index]?.pulse??null,peakResponse:peak.value,ssOverP1};
    const threshold=ssOverP1+.5*amplitude,start=rows[0].time;
    for(let i=peak.index;i<rows.length;i++)if(rows[i].releaseNorm<=threshold)return {t50Dep:rows[i].time-start,n50Dep:rows[i].pulse,reached:true,peakPulse:rows[peak.index].pulse,peakResponse:peak.value,ssOverP1};
    return {t50Dep:null,n50Dep:null,reached:false,peakPulse:rows[peak.index].pulse,peakResponse:peak.value,ssOverP1};
  }
  function computeTimeToSteadyState(train,frequencyHz=train?.frequency,tolerance=.1,minConsecutive=3){
    const rows=train?.rows||[],ssOverP1=computeSteadyStateDepth(train),peak=trainPeak(rows),amplitude=peak.value-ssOverP1;
    if(rows.length<minConsecutive||!(frequencyHz>0)||!(amplitude>1e-9))return {tSS:null,nSS:null,reached:false,tolerance,minConsecutive};
    const band=Math.max(Math.abs(amplitude)*tolerance,1e-9),start=rows[0].time;
    for(let i=peak.index;i<=rows.length-minConsecutive;i++){
      const inBand=rows.slice(i).every(row=>Math.abs(row.releaseNorm-ssOverP1)<=band);
      if(inBand)return {tSS:rows[i].time-start,nSS:rows[i].pulse,reached:true,tolerance,minConsecutive};
    }
    return {tSS:null,nSS:null,reached:false,tolerance,minConsecutive};
  }
  function computeTrainPhenotype(train,frequencyHz=train?.frequency,tolerance=.1){
    return {ssOverP1:computeSteadyStateDepth(train),...computeDepressionT50(train,frequencyHz),...computeTimeToSteadyState(train,frequencyHz,tolerance,3)};
  }
  function recoverySeries(curve,maxTime=4){
    const points=(curve?.points||curve||[]).filter(p=>Number.isFinite(p.interval)&&Number.isFinite(p.recovered)&&p.interval>=0&&p.interval<=maxTime).map(p=>({time:p.interval,value:clamp(p.recovered,0,1)})).sort((a,b)=>a.time-b.time);
    if(!points.length||points[0].time>0)points.unshift({time:0,value:0});
    const all=(curve?.points||curve||[]).filter(p=>Number.isFinite(p.interval)&&Number.isFinite(p.recovered)).sort((a,b)=>a.interval-b.interval),exact=all.find(p=>p.interval===maxTime);
    if(exact&&!points.some(p=>p.time===maxTime))points.push({time:maxTime,value:clamp(exact.recovered,0,1)});
    if(!points.some(p=>p.time===maxTime)){
      const before=[...all].reverse().find(p=>p.interval<maxTime),after=all.find(p=>p.interval>maxTime);
      if(before&&after){const f=(maxTime-before.interval)/(after.interval-before.interval);points.push({time:maxTime,value:clamp(before.recovered+f*(after.recovered-before.recovered),0,1)})}
    }
    return points.sort((a,b)=>a.time-b.time);
  }
  function computeRecoveryAUC(curve,maxTime=4){
    const points=recoverySeries(curve,maxTime);if(!points.some(p=>p.time===maxTime))return null;
    let area=0;for(let i=1;i<points.length;i++)area+=(points[i].time-points[i-1].time)*(points[i].value+points[i-1].value)/2;return area;
  }
  function computeRecoveryThresholdTime(curve,threshold,maxTime=4){
    const points=recoverySeries(curve,maxTime);for(let i=1;i<points.length;i++)if(points[i].value>=threshold){const a=points[i-1],b=points[i],span=b.value-a.value,f=span>1e-12?(threshold-a.value)/span:1;return {time:a.time+clamp(f,0,1)*(b.time-a.time),censored:false,reached:true}}
    return {time:maxTime,censored:true,reached:false};
  }
  function computeRecoveryPhenotype(curve,maxTime=4){
    const valueAt=time=>{const points=curve?.points||curve||[],exact=points.find(p=>Math.abs(p.interval-time)<1e-9);return exact?.recovered??null},t50=computeRecoveryThresholdTime(curve,.5,maxTime),t80=computeRecoveryThresholdTime(curve,.8,maxTime);
    return {r025:valueAt(.25),r05:valueAt(.5),r1:valueAt(1),r2:valueAt(2),r4:valueAt(4),auc0to4:computeRecoveryAUC(curve,maxTime),t50:t50.time,t80:t80.time,t50Censored:t50.censored,t80Censored:t80.censored,maxTime};
  }
  function extractPhenotype(p,protocol,ids,full=false){
    const values={},trains={},trainMetrics={},frequencies=new Set();
    const train=f=>{if(!trains[f])trains[f]=C.simulateTrain(f,protocol.pulses,p);return trains[f]};
    for(const id of ids){if(id.startsWith('ss:')||id.startsWith('cum:')||id.startsWith('train:'))frequencies.add(+id.split(':')[1]);if(id.startsWith('dep:'))frequencies.add(+id.split(':')[2]);if(['p1','ppr'].includes(id))frequencies.add(protocol.pprFrequency)}
    if(full)protocol.frequencies.forEach(f=>frequencies.add(f));
    for(const f of frequencies){const t=train(f),r=t.rows,m1=r[0].release,metric=computeTrainPhenotype(t,f,protocol.steadyTolerance??.1),timeBoundary=r.at(-1).time-r[0].time,pulseBoundary=r.at(-1).pulse;trainMetrics[f]=metric;values['ss:'+f]=metric.ssOverP1;values['dep:t50:'+f]=metric.t50Dep??timeBoundary;values['dep:tss:'+f]=metric.tSS??timeBoundary;values['dep:n50:'+f]=metric.n50Dep??pulseBoundary;values['dep:nss:'+f]=metric.nSS??pulseBoundary;values['cum:'+f]=r.reduce((s,x)=>s+x.release,0);if(f===protocol.pprFrequency){values.p1=m1;values.ppr=r[1].release/m1}}
    for(const id of ids)if(id.startsWith('train:')){const [,f,j,unit]=id.split(':'),r=train(+f).rows[+j-1];values[id]=r?(unit==='raw'?r.release:r.releaseNorm):NaN}
    const needsRecoverySummary=ids.some(x=>x.startsWith('rec:')),standardRecovery=needsRecoverySummary?[.25,.5,1,2,4]:[],delays=[...new Set(ids.filter(x=>x.startsWith('r:')).map(x=>+x.slice(2)).concat(standardRecovery,full?protocol.delays:[]))];
    let recovery=null,recoveryMetrics=null;if(delays.length){recovery=C.simulateRecovery(train(protocol.recoveryFrequency),delays,p,'last5');for(const row of recovery.points)values['r:'+row.interval]=row.recovered;recoveryMetrics=computeRecoveryPhenotype(recovery,4);values['rec:auc4']=recoveryMetrics.auc0to4;values['rec:t50']=recoveryMetrics.t50;values['rec:t80']=recoveryMetrics.t80}
    const resting=C.steadyState(p);values.restTS=resting.ts/p.nSites;values.frp=resting.ls+resting.ts;values.ts=resting.ts;values.ls=resting.ls;values.es=resting.es;values['ca:rest']=p.caRestActual;values['ca:ap_global']=p.caAmplGlobal;values['ca:ap_local']=p.caAmplLocal;
    return {values,resting,trainMetrics,recoveryMetrics, ...(full?{trains,recovery}: {})};
  }
  function compileTargets(targets,wt,qualRanges=QUAL){
    return targets.filter(t=>t.weight>0&&t.kind!=='unknown').map(t=>{
      const ref=wt[t.id];if(!Number.isFinite(ref))throw Error('Undefined reference readout: '+t.id);
      let low,high,scale;
      if(t.kind==='qualitative'){
        if(!(ref>0))throw Error('Use an absolute quantitative target for '+t.id+': WT is zero or negative.');
        const range=qualRanges[t.level];if(!range||!Number.isFinite(range[0])||range[0]<0||!(range[1]===null||range[1]>range[0]))throw Error('Invalid qualitative range');low=ref*range[0];high=range[1]===null?Infinity:ref*range[1];scale=Math.max(Math.abs(ref)*.1,1e-8);
      }else{
        const factor=t.unit==='ratio'?ref:1;if(t.unit==='ratio'&&ref<=0)throw Error('Ratio target needs positive WT: '+t.id);
        if(!Number.isFinite(t.mean))throw Error('Invalid mean for '+t.id);
        const mean=t.mean*factor;
        if(t.uncertainty==='ci'){
          if(!(Number.isFinite(t.low)&&Number.isFinite(t.high)&&t.low<t.high&&t.mean>=t.low&&t.mean<=t.high))throw Error('Invalid 95% interval for '+t.id);
          scale=Math.abs((t.high-t.low)*factor)/3.92;
        }else if(['sd','sem','scale'].includes(t.uncertainty)){
          if(!(t.error>0&&Number.isFinite(t.error)))throw Error('Uncertainty must be positive for '+t.id);scale=t.error*Math.abs(factor);
        }else scale=Math.max(Math.abs(mean),Math.abs(ref),1e-12)*.1;
        low=high=mean;
      }
      return {...t,low,high,scale,block:t.block||targetBlock(t.id)};
    });
  }
  function targetBlock(id){if(id==='ppr')return'ppr';if(id.startsWith('ss:'))return'trainDepth';if(id.startsWith('dep:'))return'trainKinetics';if(id.startsWith('train:'))return'trainCurve';if(id.startsWith('r:')||id.startsWith('rec:'))return'recoveryKinetics';if(id.startsWith('ca:'))return'calcium';return'release'}
  function error(values,targets,blockWeights={}){
    const groups=new Map();const rows=targets.map(t=>{
      const s=values[t.id],d=!Number.isFinite(s)?Infinity:s<t.low?s-t.low:s>t.high?s-t.high:0,residual=d/t.scale;
      const contribution=t.weight*residual*residual,block=t.block||targetBlock(t.id),group=groups.get(block)||{loss:0,weight:0};group.loss+=contribution;group.weight+=t.weight;groups.set(block,group);
      return {id:t.id,simulated:s,targetLow:t.low,targetHigh:Number.isFinite(t.high)?t.high:null,scale:t.scale,residual,contribution,met:Number.isFinite(residual)&&(t.kind==='qualitative'?d===0:Math.abs(residual)<=1)};
    });let loss=0,weights=0;for(const [block,group]of groups){const w=blockWeights[block]??(['trainCurve','recoveryCurve'].includes(block)?.25:1);if(w>0){loss+=w*group.loss/Math.max(group.weight,1e-12);weights+=w}}
    return {value:weights?loss/weights:Infinity,rows,blocks:Object.fromEntries([...groups].map(([block,group])=>[block,{error:group.loss/Math.max(group.weight,1e-12),weight:blockWeights[block]??(['trainCurve','recoveryCurve'].includes(block)?.25:1)}]))};
  }
  const countChanged=(p,wt,keys)=>keys.filter(k=>Math.abs(p[k]-wt[k])>Math.max(Math.abs(wt[k])*.01,1e-12));
  function complexity(p,wt,keys){const changed=countChanged(p,wt,keys);return {changed,change:keys.reduce((s,k)=>s+Math.abs(Math.log(Math.max(p[k],1e-12)/Math.max(wt[k],1e-12))),0),count:changed.length}}
  function validate(r){
    if(!r.reference||!r.targets?.length)throw Error('Add at least one measured phenotype.');
    const p=r.protocol;if(!p||!(p.pulses>=2&&p.pulses<=200&&Number.isInteger(p.pulses)))throw Error('Stimuli must be an integer from 2 to 200.');
    if(p.steadyTolerance!==undefined&&!(Number.isFinite(p.steadyTolerance)&&p.steadyTolerance>=.01&&p.steadyTolerance<=.5))throw Error('Steady-state tolerance must be 1–50%.');
    if(!Array.isArray(p.frequencies)||!p.frequencies.length||!Array.isArray(p.delays))throw Error('Choose at least one trace frequency and a valid delay list.');
    if(![...p.frequencies,p.pprFrequency,p.recoveryFrequency].every(f=>Number.isFinite(f)&&f>=.5&&f<=333))throw Error('Frequencies must be between 0.5 and 333 Hz.');
    if(!p.delays.every(t=>Number.isFinite(t)&&t>0&&t<=16))throw Error('Recovery delays must be positive and at most 16 s.');
    if(!(r.reference.pRel0>0&&r.reference.pRel0<=1&&r.reference.nSites>0))throw Error('The reference needs positive N and a valid initial fusion probability.');
    for(const k of SIM_NUMERIC)if(!Number.isFinite(r.reference[k]))throw Error('Missing or non-finite reference parameter: '+k);
    if(![1,2,3].includes(r.reference.model)||r.reference.odeStep<1e-6)throw Error('Invalid reference model or integration step.');
    safeParameters(r.reference);
    for(const [k,v] of Object.entries(r.locks||{}))if(!Object.keys(NAMES).includes(k)||!Number.isFinite(v)||v<0||k==='pRel0'&&(v<=0||v>1)||k==='nSites'&&v<=0)throw Error('Invalid locked value: '+k);
    if(!r.keys.length||r.keys.length>24||new Set(r.keys).size!==r.keys.length)throw Error('Choose 1–24 distinct inferable parameters.');
    const fractions=['fractionTsl','caFracSlow','yInc','zDec','zMin'];
    for(const k of r.keys){if(!NAMES[k]||k in (r.locks||{}))throw Error('Locked or unsupported inferable parameter: '+k);const b=r.bounds[k];if(!b||b.length!==2||!b.every(Number.isFinite)||b[0]<0||b[1]<=b[0]||(!k.startsWith('slope')&&!fractions.includes(k)&&b[0]<=0)||fractions.includes(k)&&b[1]>1||k==='pRel0'&&b[1]>1||k==='yMax'&&b[0]<1)throw Error('Invalid bounds: '+k)}
    for(const k of r.keys){const extreme={...r.reference,...r.locks,[k]:r.bounds[k][k.startsWith('tau')||k.startsWith('caTau')?0:1]};safeParameters(extreme)}
    if(!['fast','standard','exhaustive'].includes(r.settings.mode)||!(r.settings.maxChanged>=1&&r.settings.maxChanged<=r.keys.length)||!Number.isInteger(r.settings.maxChanged))throw Error('Invalid search settings.');
    if(![r.settings.lambdaChange,r.settings.lambdaCount].every(v=>Number.isFinite(v)&&v>=0)||!(r.settings.nearTolerance>=0&&r.settings.nearTolerance<=1))throw Error('Invalid complexity/tolerance setting.');
    if(r.blockWeights&&Object.values(r.blockWeights).some(v=>!Number.isFinite(v)||v<0))throw Error('Phenotype block weights must be non-negative.');
    if(r.settings.iterations!==undefined&&!(Number.isInteger(r.settings.iterations)&&r.settings.iterations>=1&&r.settings.iterations<=200))throw Error('Invalid iteration budget.');
    if(r.settings.seed!==undefined&&!(Number.isInteger(r.settings.seed)&&r.settings.seed>=0&&r.settings.seed<=4294967295))throw Error('Seed must be an integer between 0 and 4294967295.');
    for(const t of r.targets){if(!Number.isFinite(t.weight)||t.weight<0)throw Error('Weights must be non-negative');if(!['p1','ppr','restTS','frp','ts','ls','es','ca:rest','ca:ap_global','ca:ap_local','rec:auc4','rec:t50','rec:t80'].includes(t.id)&&! /^(ss|cum|r):[\d.]+$/.test(t.id)&&!/^dep:(t50|tss|n50|nss):[\d.]+$/.test(t.id)&&!/^train:[\d.]+:\d+:(norm|raw)$/.test(t.id))throw Error('Invalid readout: '+t.id);if(t.id==='frp'&&t.kind==='quantitative'&&!(t.mean>0))throw Error('Measured FRP must be positive');if(['ts','ls','es'].includes(t.id)&&t.kind==='quantitative'&&!(t.mean>=0))throw Error('Pool target must be non-negative');if(t.id.startsWith('r:')&&!(+t.id.slice(2)>0&&+t.id.slice(2)<=16))throw Error('Invalid recovery delay');if(/^(ss|cum|train):/.test(t.id)&&!(+t.id.split(':')[1]>=.5&&+t.id.split(':')[1]<=333))throw Error('Invalid readout frequency');if(t.id.startsWith('dep:')&&!(+t.id.split(':')[2]>=.5&&+t.id.split(':')[2]<=333))throw Error('Invalid depression frequency');if(t.id.startsWith('train:')&&!(+t.id.split(':')[2]>=1&&+t.id.split(':')[2]<=p.pulses))throw Error('Train exceeds stimulus count')}
  }
  function coverage(e,wtError){return wtError<1e-12?null:clamp(1-e/wtError,0,1)}
  function optimizeSubset(subset,seedParams,context,budget,startCount){
    const {r,evaluate,rng}=context;
    const encode=(k,v)=>Math.log1p(v/Math.max(r.reference[k],1e-9)),decode=(k,z)=>Math.max(r.reference[k],1e-9)*Math.expm1(z);
    const limits=subset.map(k=>r.bounds[k].map(v=>encode(k,v)));
    const bound=x=>x.map((v,i)=>clamp(v,...limits[i]));
    const evalX=x=>{const p={...context.fixed};subset.forEach((k,i)=>p[k]=decode(k,x[i]));return {x,solution:evaluate(p)}};
    const starts=[subset.map(k=>encode(k,seedParams[k])),subset.map(k=>encode(k,r.reference[k]))];
    while(starts.length<startCount)starts.push(limits.map(([a,b])=>a+rng()*(b-a)));
    const found=[];
    for(const start of starts.slice(0,startCount)){
      const initial=bound(start),simplex=[evalX(initial)];
      subset.forEach((_,i)=>{const x=initial.slice();const step=.15*(limits[i][1]-limits[i][0]);x[i]+=x[i]+step>limits[i][1]?-step:step;simplex.push(evalX(bound(x)))});
      for(let iter=0;iter<budget;iter++){
        simplex.sort((a,b)=>a.solution.objective-b.solution.objective);const n=subset.length,best=simplex[0],worst=simplex[n];
        if(n===0)break;
        const center=subset.map((_,i)=>simplex.slice(0,n).reduce((s,v)=>s+v.x[i],0)/n);
        const reflected=evalX(bound(center.map((v,i)=>2*v-worst.x[i])));
        if(reflected.solution.objective<best.solution.objective){const expanded=evalX(bound(center.map((v,i)=>v+2*(reflected.x[i]-v))));simplex[n]=expanded.solution.objective<reflected.solution.objective?expanded:reflected}
        else if(reflected.solution.objective<simplex[n-1].solution.objective)simplex[n]=reflected;
        else{const outside=reflected.solution.objective<worst.solution.objective,contracted=evalX(bound(center.map((v,i)=>v+.5*((outside?reflected.x[i]:worst.x[i])-v))));if(contracted.solution.objective<(outside?reflected.solution.objective:worst.solution.objective))simplex[n]=contracted;else for(let i=1;i<=n;i++)simplex[i]=evalX(bound(simplex[i].x.map((v,j)=>best.x[j]+.5*(v-best.x[j]))))}
      }
      simplex.sort((a,b)=>a.solution.objective-b.solution.objective);found.push(simplex[0].solution);
    }
    return found;
  }
  function infer(r,progress=()=>{}){
    validate(r);r=clone(r);const wt={...r.reference},fixed={...wt,...r.locks},ids=[...new Set(r.targets.filter(t=>t.weight>0&&t.kind!=='unknown').map(t=>t.id))];
    const wtSim=extractPhenotype(wt,r.protocol,ids),targets=compileTargets(r.targets,wtSim.values,r.qualitativeRanges||QUAL);
    if(!targets.length)throw Error('All measurements are missing or zero-weight.');
    const wtError=error(wtSim.values,targets,r.blockWeights).value,cache=new Map(),pool=[],subsetResults=[];let evaluations=0,seed=(r.settings.seed??2026)>>>0;
    const rng=()=>{seed=(1664525*seed+1013904223)>>>0;return seed/4294967296};
    const evaluate=params=>{
      const key=r.keys.map(k=>params[k].toPrecision(12)).join('|');if(cache.has(key))return cache.get(key);
      const sim=extractPhenotype(params,r.protocol,ids),fit=error(sim.values,targets,r.blockWeights),c=complexity(params,wt,r.keys),s={params:{...params},values:sim.values,error:fit.value,rows:fit.rows,blockErrors:fit.blocks,...c};
      s.objective=s.error+r.settings.lambdaChange*c.change+r.settings.lambdaCount*c.count;s.coverage=coverage(s.error,wtError);
      if(!Number.isFinite(s.error))s.objective=Infinity;
      cache.set(key,s);pool.push(s);evaluations++;return s;
    };
    const baseline=evaluate(fixed),context={r,evaluate,fixed,rng};let warm=baseline;
    const mode=r.settings.mode,budget=r.settings.iterations??({fast:12,standard:28,exhaustive:55}[mode]),starts=mode==='exhaustive'?3:mode==='standard'?2:1;
    const sensitivity=[];
    // Stage 1 screens independent perturbations before subset optimization.
    for(const k of r.keys){const lo={...fixed,[k]:clamp(wt[k]*.7,...r.bounds[k])},hi={...fixed,[k]:clamp(wt[k]*1.3,...r.bounds[k])};const a=evaluate(lo),b=evaluate(hi);sensitivity.push({key:k,effect:Object.fromEntries(ids.map(id=>[id,(b.values[id]-a.values[id])/Math.max(Math.abs(wtSim.values[id]),.01)]))})}
    let stageKeys=r.keys;
    for(let size=1;size<=r.settings.maxChanged;size++){
      if(mode==='fast'&&size>2)break;
      if(size===2&&r.keys.length>8){
        const singles=subsetResults.filter(s=>s.subset.length===1).sort((a,b)=>a.best.objective-b.best.objective);
        const ceiling=mode==='fast'?6:mode==='standard'?10:Math.min(r.keys.length,14);
        const pool=new Set(singles.slice(0,ceiling).map(s=>s.subset[0]));
        for(const family of Object.values(FAMILIES)){const ranked=singles.filter(s=>family.includes(s.subset[0]));if(ranked.length)pool.add(ranked[0].subset[0])}
        stageKeys=r.keys.filter(k=>pool.has(k));
      }else if(mode==='fast'&&size===2)stageKeys=subsetResults.filter(s=>s.subset.length===1).sort((a,b)=>a.best.error-b.best.error).slice(0,4).map(s=>s.subset[0]);
      if(size>=3&&stageKeys.length>8){
        const ranked=subsetResults.filter(s=>s.subset.length===2).sort((a,b)=>a.best.objective-b.best.objective).slice(0,8).flatMap(s=>s.subset);
        const top=new Set(ranked.slice(0,8));stageKeys=r.keys.filter(k=>top.has(k));
      }
      const groups=subsets(stageKeys,size);
      for(let index=0;index<groups.length;index++){
        const subset=groups[index],start={...fixed};subset.forEach(k=>start[k]=warm.params[k]);
        const found=optimizeSubset(subset,start,context,budget,starts),best=found.sort((a,b)=>a.error-b.error)[0];subsetResults.push({subset,best});if(best.error<warm.error)warm=best;
        progress({stage:size,subset:index+1,total:groups.length,evaluations,bestError:warm.error});
      }
      if(mode!=='exhaustive'&&size>=2&&warm.rows.every(x=>x.met)&&warm.coverage>=.995)break;
    }
    const finite=pool.filter(s=>Number.isFinite(s.error));if(!finite.length)throw Error('No finite simulation satisfies the supplied setup.');
    finite.sort((a,b)=>a.error-b.error);const numerical=finite[0],threshold=numerical.error+Math.max(1e-6,numerical.error*r.settings.nearTolerance);
    const qualifying=finite.filter(s=>r.settings.minimalCriterion==='coverage'?s.coverage!==null&&s.coverage>=.9:s.error<=threshold);
    const minimal=qualifying.sort((a,b)=>a.count-b.count||a.objective-b.objective||a.error-b.error)[0]||numerical;
    const near=finite.filter(s=>s.error<=threshold),distinct=[];
    // Remove near-duplicate optimizer samples before computing descriptive ranges.
    for(const s of near){if(!distinct.some(t=>r.keys.every(k=>Math.abs(Math.log((s.params[k]+1e-12)/(t.params[k]+1e-12)))<.08)))distinct.push(s);if(distinct.length>=200)break}
    const alternativePool=[...subsetResults.map(s=>s.best).filter(s=>s.error<=threshold||s.rows.every(t=>t.met)).sort((a,b)=>a.error-b.error),...near];
    const alternates=[];for(const s of [minimal,...alternativePool]){if(!alternates.some(t=>r.keys.every(k=>Math.abs(Math.log((s.params[k]+1e-12)/(t.params[k]+1e-12)))<.15)))alternates.push({...s,withinNearOptimal:s.error<=threshold});if(alternates.length===5)break}
    progress({stage:'attribution',evaluations});
    const distance=(a,b)=>targets.reduce((s,t)=>s+t.weight*((a[t.id]-b[t.id])/t.scale)**2,0);
    const denom=distance(baseline.values,minimal.values);
    const attribution=r.keys.map(k=>{
      const alone=evaluate({...fixed,[k]:minimal.params[k]}),ablated=evaluate({...minimal.params,[k]:wt[k]});
      const suff=coverage(alone.error,baseline.error),necessary=denom>1e-12?distance(ablated.values,minimal.values)/denom:null;
      return {key:k,sufficiency:suff,necessity:necessary,aloneError:alone.error,ablatedError:ablated.error,phenotype:Object.fromEntries(targets.map(t=>[t.id,{sufficiencyDelta:(alone.values[t.id]-baseline.values[t.id])/t.scale,ablationDelta:(minimal.values[t.id]-ablated.values[t.id])/t.scale}]))};
    });
    const sufficientCombinations=[];for(let n=0;n<=minimal.changed.length;n++)for(const keys of subsets(minimal.changed,n)){const p={...fixed};keys.forEach(k=>p[k]=minimal.params[k]);const s=evaluate(p);sufficientCombinations.push({keys,error:s.error,coverage:s.coverage})}
    const percentile=(a,f)=>{const i=(a.length-1)*f,j=Math.floor(i);return a[j]+(a[Math.ceil(i)]-a[j])*(i-j)};
    const minimalNear=distinct.filter(s=>s.count===Math.min(...distinct.map(x=>x.count)));
    const identifiability=r.keys.map(key=>{
      const v=distinct.map(s=>s.params[key]).sort((a,b)=>a-b),ratios=v.map(x=>x/Math.max(wt[key],1e-12));
      const up=ratios.filter(x=>x>1.01).length,down=ratios.filter(x=>x<.99).length,directionConsistency=Math.max(up,down)/Math.max(1,v.length),low=percentile(v,.05),high=percentile(v,.95),a=attribution.find(x=>x.key===key);
      const confidence=v.length<5?'Insufficient sampling':directionConsistency>.9&&(high-low)/Math.max(percentile(v,.5),1e-12)<.5&&(a.necessity||0)>.2?'High':directionConsistency>.7?'Medium':'Low';
      return {key,median:percentile(v,.5),low,high,direction:up>=down?'increase':'decrease',directionConsistency,inclusion:minimalNear.filter(s=>s.changed.includes(key)).length/Math.max(1,minimalNear.length),confidence,samples:v.length};
    });
    const frontier=[];for(let n=0;n<=r.keys.length;n++){const ss=finite.filter(s=>s.count<=n);if(ss.length){const s=ss.reduce((a,b)=>a.error<b.error?a:b);frontier.push({count:n,error:s.error,coverage:s.coverage})}}
    const familyOptions=[['Calcium handling',['calcium']],['Priming kinetics',['priming']],['STP / refractory',['stp','refractory']],['Calcium + priming',['calcium','priming']],['Priming + STP',['priming','stp','refractory']],['Calcium + STP',['calcium','stp','refractory']],['Combined',['calcium','priming','stp','refractory']]];
    const parameterFamily=k=>Object.keys(FAMILIES).find(f=>FAMILIES[f].includes(k));
    const familyModels=familyOptions.map(([name,families])=>{const models=subsetResults.filter(s=>{const active=[...new Set(s.best.changed.map(parameterFamily).filter(Boolean))];return active.length===families.length&&active.every(f=>families.includes(f))});const best=models.length?models.reduce((a,b)=>a.best.objective<b.best.objective?a:b).best:null;return{name,families,bestError:best?.error??null,objective:best?.objective??null,changed:best?.changed??[],searchedSubsets:models.length}});
    const validFamilies=familyModels.filter(x=>x.objective!==null),baseObjective=Math.min(...validFamilies.map(x=>x.objective)),temperature=Math.max(baseObjective*.1,.01),totalSupport=validFamilies.reduce((s,x)=>s+Math.exp(-(x.objective-baseObjective)/temperature),0);
    for(const item of familyModels)item.relativeSupport=item.objective===null?0:Math.exp(-(item.objective-baseObjective)/temperature)/totalSupport;
    const traces={wt:extractPhenotype(wt,r.protocol,ids,true),minimal:extractPhenotype(minimal.params,r.protocol,ids,true)};
    let nextExperiment=null;if(alternates.length>1){
      const possible=['ss:10','ss:50','ss:333','r:0.05','r:0.5','r:2'].filter(id=>!ids.includes(id));
      for(const id of possible){const predictions=alternates.map(s=>extractPhenotype(s.params,r.protocol,[id]).values[id]),mean=predictions.reduce((s,v)=>s+v,0)/predictions.length,spread=(Math.max(...predictions)-Math.min(...predictions))/Math.max(Math.abs(mean),.05);if(!nextExperiment||spread>nextExperiment.spread)nextExperiment={id,predictions,spread,recoveryFrequency:r.protocol.recoveryFrequency,pulses:r.protocol.pulses}}
    }
    const interpretations=minimal.changed.map(k=>`${NAMES[k]} ${minimal.params[k]>wt[k]?'increases':'decreases'}: ${Math.round(100*(minimal.params[k]/Math.max(wt[k],1e-12)-1))}% relative to WT.`);
    const restWT=traces.wt.resting.ts/wt.nSites,restM=traces.minimal.resting.ts/minimal.params.nSites;
    interpretations.push(`Resting TS fraction changes from ${(restWT*100).toFixed(1)}% to ${(restM*100).toFixed(1)}%.`);
    if(minimal.changed.includes('slope1'))interpretations.push(`The selected Ca-dependent ES→LS slope is ${minimal.params.slope1<wt.slope1?'lower':'higher'}; its contributions to replenishment and recovery are quantified below.`);
    if(minimal.changed.includes('caRestActual'))interpretations.push(`A change in actual resting Ca is consistent with altered basal occupancy while the rate-reference Ca remains ${minimal.changed.includes('caRestReference')?'also variable':'at its control value'}.`);
    if(minimal.changed.some(k=>['caAmplGlobal','caAmplLocal'].includes(k)))interpretations.push('A change in AP-associated Ca amplitude is sufficient within this model for part of the activity-dependent phenotype; it does not by itself prove altered Ca influx.');
    if(minimal.changed.includes('kf2'))interpretations.push('The candidate basal LS→TS forward rate differs from control, affecting the pre-stimulation TS pool.');
    if(minimal.changed.includes('slope2'))interpretations.push('The candidate Ca-dependent LS→TS slope differs from control, affecting activity-dependent maturation.');
    return {version:VERSION,createdAt:new Date().toISOString(),request:r,targets,wtError,constrainedBaseline:baseline,numerical,minimal,alternatives:alternates,identifiability,attribution,sensitivity,sufficientCombinations,frontier,familyModels,traces,nextExperiment,interpretations,nearOptimalCount:distinct.length,threshold,evaluations,subsets:subsetResults.map(s=>({keys:s.subset,error:s.best.error,objective:s.best.objective,params:s.best.params})),searchLimit:'Bounded finite-iteration search; large parameter spaces use hierarchical screening, not every subset. Relative model support is a search heuristic, not biological causality. Sampled ranges are not statistical confidence intervals.'};
  }
  g.MechanismEngine={VERSION,DEFAULT_KEYS,FAMILIES,NAMES,QUAL,defaultBounds,computeSteadyStateDepth,computeDepressionT50,computeTimeToSteadyState,computeTrainPhenotype,computeRecoveryAUC,computeRecoveryThresholdTime,computeRecoveryPhenotype,extractPhenotype,compileTargets,error,coverage,subsets,validate,infer};
})(globalThis);
