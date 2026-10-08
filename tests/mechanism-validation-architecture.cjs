'use strict';
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path');
const root=path.join(__dirname,'..'),dist=path.join(root,'dist');
const production=['mechanism.html','mechanism-ui.js','mechanism-engine.js','mechanism-worker.js'].map(name=>[name,fs.readFileSync(path.join(dist,name),'utf8')]);
for(const [name,text] of production){
 assert.ok(!/mechanism[_-]validation|publication-solutions|publishedSolution|Reveal publication|PDBu benchmark|ionomycin benchmark|HK benchmark/i.test(text),name+' references validation answers or benchmark UI');
}
assert.ok(!fs.existsSync(path.join(dist,'mechanism-publication-solutions.json')),'Answer file must not be in production build');
assert.ok(!fs.existsSync(path.join(dist,'mechanism-benchmarks.js')),'Benchmark loader must not be in production build');
assert.deepEqual(fs.readdirSync(dist).filter(name=>/benchmark|validation/i.test(name)),[],'Validation artifacts must not be shipped in dist');
const validation=path.join(__dirname,'mechanism_validation'),files=fs.readdirSync(validation);
for(const stem of ['jp286282_2mm_ca','jp286282_ionomycin','jp286282_pdbu','sciadv_hk']){
 assert.ok(files.includes(stem+'_input.json'));assert.ok(files.includes(stem+'_answer.json'));
 const input=JSON.parse(fs.readFileSync(path.join(validation,stem+'_input.json'),'utf8'));
 assert.ok(!('published' in input));assert.ok(!('expectedDirections' in input));
}
console.log('Mechanism architecture checks passed: production is publication-agnostic; validation inputs and answers are isolated.');
