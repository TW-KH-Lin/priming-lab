"use strict";

const PRESETS = {
  pnas2022: { nSites:2639,pRel0:.393,kf1:.402517,kb1:.18472,slope1:1.23939e7,kf2:.207254,kb2:.248014,slope2:1.27727e7,fractionTsl:.16,tauTsl:.09,kRefract:5000,caRest:50e-9,caAmpl:110e-9,caTau:.06,kHalf1:280e-9,useMM:true,yInc:.39,zDec:.4,yMax:1.32,zMin:.75,tauY:.014,tauZ:3,yPower:4.5,model:3 },
  "jp-control": { nSites:2622,pRel0:.22,kf1:.370,kb1:.221,slope1:2.245e6,kf2:.199,kb2:.253,slope2:1.014e6,fractionTsl:.10,tauTsl:.09,kRefract:2.6,caRest:50e-9,caAmpl:454e-9,caTau:.06,kHalf1:280e-9,useMM:false,yInc:.39,zDec:.4,yMax:1.31,zMin:.87,tauY:.017,tauZ:3,yPower:4,model:2 },
  "jp-iono": { nSites:2622,pRel0:.30,kf1:.370,kb1:.221,slope1:2.245e6,kf2:.199,kb2:.253,slope2:1.014e6,fractionTsl:.15,tauTsl:.09,kRefract:4.2,caRest:130e-9,caAmpl:454e-9,caTau:.06,kHalf1:280e-9,useMM:false,yInc:.39,zDec:.4,yMax:1.16,zMin:1,tauY:.017,tauZ:3,yPower:4,model:2 }
};
const COLORS=["#15a89f","#ef6b5b","#5675d6","#ca9b25","#8a5cc8","#2f8bc9","#c25b98","#65727e"];
const POOL_COLORS={ES:"#94a3b8",LS:"#18b7ae",TS:"#315f96",TSL:"#f27160",ERS:"#d6a33a"};
let currentResults=[];

const $=s=>document.querySelector(s);
const $$=s=>Array.from(document.querySelectorAll(s));
const clamp=(v,a,b)=>Math.min(b,Math.max(a,v));

function applyPreset(key){
  const p=PRESETS[key]; if(!p)return;
  $$('[data-param]').forEach(el=>{const v=p[el.dataset.param];if(v!==undefined)el.value=String(v)});
  $('#modelVariant').value=String(p.model);
}

function readParams(){
  const base={...PRESETS.pnas2022};
  $$('[data-param]').forEach(el=>{const v=Number(el.value);if(Number.isFinite(v))base[el.dataset.param]=v});
  base.model=Number($('#modelVariant').value);
  return base;
}

function steadyState(p){
  const den=p.kb1*p.kb2+p.kf1*(p.kb2+p.kf2);
  const ts=p.nSites*p.kf1*p.kf2/den,ls=p.nSites*p.kb2*p.kf1/den;
  return {es:p.nSites-ts-ls,ls,ts,tsl:0,ers:0,ca:0};
}

function rates(ca,p){
  const above=Math.max(ca-p.caRest,0);
  let k1=p.kf1+p.slope1*above;
  if(p.useMM)k1/=1+above/Math.max(p.kHalf1,1e-18);
  return [Math.max(k1,0),Math.max(p.kf2+p.slope2*above,0)];
}

function deriv(s,p){
  const [k1,k2]=rates(p.caRest+s.ca,p), b3=1/Math.max(p.tauTsl,1e-12);
  let d;
  if(p.model===1)d={ts:k2*s.ls-p.kb2*s.ts,es:-k1*s.es+p.kb1*s.ls,ls:k1*s.es+p.kb2*s.ts-(p.kb1+k2)*s.ls,tsl:0,ers:0};
  else if(p.model===3)d={ts:k2*s.ls-p.kb2*s.ts,es:-k1*s.es+p.kb1*s.ls,ls:k1*s.es+p.kb2*s.ts-(p.kb1+k2)*s.ls+b3*s.tsl,tsl:-b3*s.tsl,ers:0};
  else d={ts:k2*s.ls-p.kb2*s.ts,es:-k1*s.es+p.kb1*s.ls+p.kRefract*s.ers,ls:k1*s.es+p.kb2*s.ts-(p.kb1+k2)*s.ls+b3*s.tsl,tsl:-b3*s.tsl,ers:-p.kRefract*s.ers};
  d.ca=-s.ca/Math.max(p.caTau,1e-12);return d;
}

function addState(a,b,h){const o={};for(const k of ['es','ls','ts','tsl','ers','ca'])o[k]=a[k]+h*b[k];return o}
function rk4Step(s,h,p){
  const k1=deriv(s,p),k2=deriv(addState(s,k1,h/2),p),k3=deriv(addState(s,k2,h/2),p),k4=deriv(addState(s,k3,h),p),o={};
  for(const k of ['es','ls','ts','tsl','ers','ca'])o[k]=s[k]+h*(k1[k]+2*k2[k]+2*k3[k]+k4[k])/6;
  for(const k of ['es','ls','ts','tsl','ers'])o[k]=Math.max(0,o[k]);
  const total=o.es+o.ls+o.ts+o.tsl+o.ers,scale=p.nSites/Math.max(total,1e-18);
  for(const k of ['es','ls','ts','tsl','ers'])o[k]*=scale;
  return o;
}

function integrate(s,duration,p){
  const maxStep=.001,n=Math.max(1,Math.ceil(duration/maxStep)),h=duration/n;
  for(let i=0;i<n;i++)s=rk4Step(s,h,p);return s;
}

function simulate(freq,pulses,p){
  let s=steadyState(p),y=1,z=1,prevT=0;
  const rows=[];
  for(let i=0;i<pulses;i++){
    const t=(i+1)/freq,isi=t-prevT;
    s=integrate(s,isi,p);
    if(i>0){y=1+(y-1)*Math.exp(-isi/p.tauY);z=1+(z-1)*Math.exp(-isi/p.tauZ)}
    const pRel=clamp(p.pRel0*Math.pow(y,p.yPower)*z,0,1);
    const ts0=s.ts,tsl0=s.tsl,ls0=s.ls;
    let release;
    if(p.model===1){release=ts0*pRel;s.ts=ts0*(1-pRel);s.es+=release}
    else{release=(ts0+tsl0)*pRel;s.ts=ts0*(1-pRel);s.ls=ls0*(1-p.fractionTsl);s.tsl=tsl0*(1-pRel)+ls0*p.fractionTsl;if(p.model===3)s.es+=release;else s.ers+=release}
    s.ca+=p.caAmpl*y;
    const total=s.es+s.ls+s.ts+s.tsl+s.ers,scale=p.nSites/Math.max(total,1e-18);
    for(const k of ['es','ls','ts','tsl','ers'])s[k]*=scale;
    rows.push({pulse:i+1,time:t,release,pRel,ES:s.es,LS:s.ls,TS:s.ts,TSL:s.tsl,ERS:s.ers,Ca:p.caRest+s.ca});
    y+=p.yInc*(p.yMax-y);z-=p.zDec*(z-p.zMin);prevT=t;
  }
  const first=Math.max(rows[0].release,1e-18);rows.forEach(r=>r.releaseNorm=r.release/first);
  return {frequency:freq,rows};
}

function parseFrequencies(){
  const vals=$('#frequencies').value.split(/[;,\s]+/).map(Number).filter(v=>Number.isFinite(v)&&v>0).slice(0,8);
  return [...new Set(vals)];
}

function svgEl(tag,attrs={},text=''){const e=document.createElementNS('http://www.w3.org/2000/svg',tag);Object.entries(attrs).forEach(([k,v])=>e.setAttribute(k,v));if(text)e.textContent=text;return e}
function drawChart(svg,series,{yMin=0,yMax=null,yLabel='',xLabel='Pulse'}={}){
  svg.replaceChildren();const W=svg.viewBox.baseVal.width||900,H=svg.viewBox.baseVal.height||390,m={l:50,r:18,t:18,b:42};
  const all=series.flatMap(s=>s.values.filter(Number.isFinite));if(!all.length){svg.append(svgEl('text',{x:W/2,y:H/2,class:'empty-chart'},'Run a simulation to view results'));return}
  const maxLen=Math.max(...series.map(s=>s.values.length));yMax=yMax??Math.max(...all)*1.08;if(yMax<=yMin)yMax=yMin+1;
  const x=i=>m.l+(maxLen===1?0:i/(maxLen-1))*(W-m.l-m.r),y=v=>H-m.b-(v-yMin)/(yMax-yMin)*(H-m.t-m.b);
  for(let i=0;i<=4;i++){const yy=m.t+i*(H-m.t-m.b)/4,val=yMax-i*(yMax-yMin)/4;svg.append(svgEl('line',{x1:m.l,y1:yy,x2:W-m.r,y2:yy,class:'grid-line'}));svg.append(svgEl('text',{x:m.l-8,y:yy+4,'text-anchor':'end',class:'axis-label'},formatTick(val)))}
  for(let i=0;i<5;i++){const idx=Math.round(i*(maxLen-1)/4),xx=x(idx);svg.append(svgEl('text',{x:xx,y:H-17,'text-anchor':'middle',class:'axis-label'},String(idx+1)))}
  svg.append(svgEl('line',{x1:m.l,y1:m.t,x2:m.l,y2:H-m.b,class:'axis-line'}));svg.append(svgEl('line',{x1:m.l,y1:H-m.b,x2:W-m.r,y2:H-m.b,class:'axis-line'}));
  svg.append(svgEl('text',{x:(m.l+W-m.r)/2,y:H-2,'text-anchor':'middle',class:'axis-label'},xLabel));
  if(yLabel)svg.append(svgEl('text',{x:13,y:H/2,transform:`rotate(-90 13 ${H/2})`,'text-anchor':'middle',class:'axis-label'},yLabel));
  series.forEach(s=>{const d=s.values.map((v,i)=>(i?'L':'M')+x(i).toFixed(2)+' '+y(v).toFixed(2)).join(' ');svg.append(svgEl('path',{d,class:'chart-path',stroke:s.color}))});
}
function formatTick(v){if(Math.abs(v)>=1000)return Math.round(v).toLocaleString();if(Math.abs(v)<.01&&v!==0)return v.toExponential(1);return Number(v.toFixed(2)).toString()}

function render(){
  const legend=$('#releaseLegend');legend.replaceChildren();currentResults.forEach((r,i)=>{const s=document.createElement('span'),dot=document.createElement('i');dot.style.background=COLORS[i];s.append(dot,document.createTextNode(`${r.frequency} Hz`));legend.append(s)});
  drawChart($('#releaseChart'),currentResults.map((r,i)=>({values:r.rows.map(x=>x.releaseNorm),color:COLORS[i]})),{yLabel:'Normalized release'});
  const sel=$('#detailFrequency'),old=sel.value;sel.replaceChildren();currentResults.forEach((r,i)=>{const o=document.createElement('option');o.value=String(i);o.textContent=`${r.frequency} Hz`;sel.append(o)});if(old&&Number(old)<currentResults.length)sel.value=old;renderDetail();
  const r=currentResults[0],rows=r.rows,first=rows[0].release,ppr=rows.length>1?rows[1].release/first:NaN,last=rows.slice(-Math.min(5,rows.length)).reduce((a,b)=>a+b.release,0)/Math.min(5,rows.length)/first;
  const cards=$$('#metrics article');cards[0].querySelector('strong').textContent=formatTick(first);cards[1].querySelector('strong').textContent=Number(ppr.toFixed(3));cards[2].querySelector('strong').textContent=Number(last.toFixed(3));
}
function renderDetail(){
  if(!currentResults.length)return;const r=currentResults[Number($('#detailFrequency').value)||0],rows=r.rows;
  drawChart($('#poolChart'),Object.entries(POOL_COLORS).map(([k,c])=>({values:rows.map(x=>x[k]),color:c})),{yLabel:'Sites'});
  drawChart($('#probChart'),[{values:rows.map(x=>x.pRel),color:'#f27160'}],{yMin:0,yMax:1,yLabel:'Probability'});
}
function run(){
  const freqs=parseFrequencies(),pulses=clamp(Math.round(Number($('#pulses').value)||40),2,200);if(!freqs.length){$('#status').innerHTML='<span style="background:#ef6b5b"></span> Enter a valid frequency';return}
  $('#status').innerHTML='<span></span> Computing';
  requestAnimationFrame(()=>{try{const p=readParams();currentResults=freqs.map(f=>simulate(f,pulses,p));render();$('#exportButton').disabled=false;$('#status').innerHTML='<span></span> Simulation complete'}catch(e){console.error(e);$('#status').innerHTML='<span style="background:#ef6b5b"></span> Check parameters'}});
}
function exportCsv(){
  if(!currentResults.length)return;const head=['frequency_hz','pulse','time_s','release','release_normalized','p_release','ES','LS','TS','TSL','ERS','Ca_M'];const lines=[head.join(',')];
  currentResults.forEach(r=>r.rows.forEach(x=>lines.push([r.frequency,x.pulse,x.time,x.release,x.releaseNorm,x.pRel,x.ES,x.LS,x.TS,x.TSL,x.ERS,x.Ca].join(','))));
  const blob=new Blob([lines.join('\n')],{type:'text/csv'}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download='priminglab_simulation.csv';a.click();setTimeout(()=>URL.revokeObjectURL(a.href),1000);
}

$('#preset').addEventListener('change',e=>{if(e.target.value!=='custom')applyPreset(e.target.value);run()});
$$('[data-param],#modelVariant,#pulses,#frequencies').forEach(el=>el.addEventListener('change',()=>{$('#preset').value='custom'}));
$('#runButton').addEventListener('click',run);$('#exportButton').addEventListener('click',exportCsv);$('#detailFrequency').addEventListener('change',renderDetail);
applyPreset('pnas2022');run();
