(function (global) {
  "use strict";
  const names = {ppr:"Paired-pulse ratio",depression:"Depression at train end",steady:"Steady-state release",recovery:"Recovery at probe delay"};
  const rateNames = {pRel0:"pFusion",nSites:"Total sites N",kf1:"k₁",kb1:"b₁",kf2:"k₂",kb2:"b₂",slope1:"Ca slope σ₁",slope2:"Ca slope σ₂"};
  const fmt = v => Number.isFinite(v) ? (Math.abs(v) >= 1000 ? v.toLocaleString("en",{maximumFractionDigits:1}) : v.toFixed(3)) : "Not defined";
  let hooks, panel, base, worker, result, applying = false, lastApplied = null;
  const $ = selector => panel.querySelector(selector);
  function stop() { if (worker) worker.terminate(); worker = null; $("#guideRun").textContent = "Find a matching simulation"; }
  function capture(params) {
    stop(); base = {...params}; result = null; lastApplied = null;
    $("#guideResults").replaceChildren(); $("#guideUndo").hidden = true;
    $("#guideBaseline").textContent = `Baseline: ${hooks.baselineName()}`;
    $("#guideStatus").textContent = "Choose your observations, then compare simulated outcomes.";
  }
  function observation(key, label) {
    return `<label class="guide-observation"><span>${label}</span><select data-observation="${key}" aria-label="${label}"><option value="none">No data</option><option value="same">Unchanged</option><option value="up">${key === "recovery" ? "More recovered" : key === "depression" ? "More depression" : "Higher"}</option><option value="down">${key === "recovery" ? "Less recovered" : key === "depression" ? "Less depression" : "Lower"}</option><option value="ratio">Measured ratio…</option></select><input data-ratio="${key}" type="number" min="0.01" max="20" step="any" placeholder="e.g. 0.50" aria-label="${label}, condition / control ratio" hidden></label>`;
  }
  function init(callbacks) {
    hooks = callbacks; panel = document.querySelector(".phenotype-panel");
    panel.className = "phenotype-panel phenotype-guide";
    panel.innerHTML = `<summary>Match observed changes</summary>
      <p id="guideBaseline" class="guide-baseline"></p>
      <fieldset><legend>1 · What did you measure?</legend>
        <p class="guide-help">Condition ÷ control: 1 = unchanged, 2 = doubled. Leave unknown measurements blank.</p>
        <div class="guide-measurements"><label>pFusion ratio<input id="guidePFusion" type="number" min="0.01" max="20" step="any" placeholder="Baseline" aria-label="Measured pFusion ratio"></label><label>Total pool ratio<input id="guidePool" type="number" min="0.01" max="20" step="any" placeholder="Baseline" aria-label="Measured total pool ratio"></label></div>
        <p class="guide-help">Unknown pFusion and total sites N stay at baseline. Total pool here means N, not the readily releasable pool.</p>
        ${observation("first", "First EPSC")}
      </fieldset>
      <fieldset><legend>2 · What changed?</legend>${Object.entries(names).map(([key,label]) => observation(key,label)).join("")}</fieldset>
      <details class="guide-protocol"><summary id="guideProtocolSummary">Measurement conditions</summary><div class="guide-measurements"><label>Train frequency (Hz)<input id="guideFrequency" type="number" min="0.5" max="333" value="200" step="any"></label><label>Stimuli<input id="guidePulses" type="number" min="2" max="200" value="40" step="1"></label><label>Probe delay (s)<input id="guideDelay" type="number" min="0.01" max="16" value="0.5" step="any"></label></div><p class="guide-help">PPR = pulse 2 / pulse 1. Steady state = mean last 5 releases (absolute). Depression = 1 − last-5 / first. Recovery = (probe − last-5) / (first − last-5). Recovery is the amount at this delay, not a fitted speed.</p></details>
      <button id="guideRun" class="button primary" type="button">Find a matching simulation</button>
      <p id="guideStatus" role="status" aria-live="polite"></p><section id="guideResults" aria-label="Comparison with observations"></section>
      <button id="guideUndo" class="guide-undo" type="button" hidden>Restore baseline</button>`;
    const toggle = document.querySelector("#togglePhenotype");
    toggle.textContent = "Describe my phenotype";
    toggle.addEventListener("click", () => {
      if (panel.hidden) {panel.hidden = false; panel.open = true; capture(hooks.params());}
      else {restore(); panel.hidden = true;}
      toggle.setAttribute("aria-expanded", String(!panel.hidden));
      toggle.textContent = panel.hidden ? "Describe my phenotype" : "Close phenotype guide";
    });
    panel.addEventListener("input", event => {
      if (event.target.matches("input,select")) {
        stop();
        if (event.target.dataset.observation) $( `[data-ratio="${event.target.dataset.observation}"]`).hidden = event.target.value !== "ratio";
        $("#guideStatus").textContent = result ? "Inputs changed. Run again to compare with these observations." : "Ready to compare with your observations.";
        if (result) $("#guideResults").classList.add("guide-stale");
        protocolSummary();
      }
    });
    $("#guideRun").addEventListener("click", () => worker ? (stop(), $("#guideStatus").textContent = "Comparison cancelled.") : start());
    $("#guideUndo").addEventListener("click", restore);
    protocolSummary();
  }
  function protocolSummary() {
    $("#guideProtocolSummary").textContent = `${$("#guideFrequency").value} Hz · ${$("#guidePulses").value} stimuli · probe at ${$("#guideDelay").value} s`;
  }
  function request() {
    const read = id => { const input = $(id); if (!input.checkValidity()) throw Error("Check the measurement values and protocol limits."); return input.value.trim() === "" ? null : Number(input.value); };
    const targets = {};
    panel.querySelectorAll("[data-observation]").forEach(select => {
      const key = select.dataset.observation;
      targets[key] = {mode:select.value, value:select.value === "ratio" ? read(`[data-ratio="${key}"]`) : null};
    });
    return {protocol:{frequency:read("#guideFrequency"),pulses:read("#guidePulses"),delay:read("#guideDelay")},measured:{pFusion:read("#guidePFusion"),pool:read("#guidePool")},targets};
  }
  function start() {
    try {
      const observations = request();
      global.PhenotypeEngine.validate(base, observations);
      worker = new Worker("phenotype-worker.js");
      $("#guideRun").textContent = "Cancel comparison";
      $("#guideStatus").textContent = "Simulating candidate mechanisms…";
      const activeWorker = worker;
      worker.onmessage = ({data}) => {
        if (worker !== activeWorker) return;
        if (data.type === "progress") $("#guideStatus").textContent = `Comparing ${data.done} of ${data.total} candidates…`;
        else if (data.type === "error") {stop(); $("#guideStatus").textContent = data.message;}
        else if (data.type === "result") {
          stop(); result = data.result;
          applying = true;
          try {hooks.apply(result.params, result); lastApplied = JSON.stringify(hooks.params());} finally {applying = false;}
          render();
        }
      };
      worker.onerror = () => { if (worker !== activeWorker) return; stop(); $("#guideStatus").textContent = "The comparison could not finish. Try again or use a shorter protocol."; };
      worker.postMessage({base,request:observations});
    } catch (error) {stop(); $("#guideStatus").textContent = error.message;}
  }
  function render() {
    const {assessment, params, tested, request} = result;
    const box = $("#guideResults"); box.classList.remove("guide-stale");
    $("#guideStatus").textContent = `${assessment.matched} of ${assessment.rows.length} observations met by the best tested candidate. Simulation plots updated.`;
    const wanted = row => row.mode === "ratio" ? `${row.target.toFixed(2)}×` : {same:"Unchanged",up:"Higher",down:"Lower"}[row.mode];
    box.innerHTML = `<h3>Your observations and simulation</h3><p class="guide-help">Ratios versus baseline. Numbers: ±5% tolerance. Arrows: at least 5% change; magnitude is unconstrained.</p><div class="guide-outcomes">${assessment.rows.map(row => `<div class="guide-outcome"><strong>${row.label}</strong><span>Wanted: ${wanted(row)}</span><span>Simulated: ${row.ratio === null ? "Not defined" : `${row.ratio.toFixed(2)}×`}</span><span class="${row.matched ? "guide-met" : "guide-miss"}">${row.error === null ? "Cannot assess" : row.matched ? row.mode === "up" || row.mode === "down" ? "Direction met" : "Within tolerance" : "Still different"}</span></div>`).join("")}</div>
      <p class="guide-help">pFusion: ${(params.pRel0 / base.pRel0).toFixed(3)}× · total sites: ${(params.nSites / base.nSites).toFixed(3)}× (fixed to your inputs or baseline).</p>
      <details class="guide-explanation"><summary>How the model explains this</summary><p>${tested} finite candidates compared. Calcium sensitivity and priming transitions can change; b₄, calcium waveforms and modulation parameters stay at baseline. This is one possible explanation, not a unique estimate of the biological rates.</p><div class="guide-param-table">${Object.entries(rateNames).map(([key,label]) => `<div><span>${label}</span><span>${fmt(base[key])} → ${fmt(params[key])}</span></div>`).join("")}</div><p>Protocol: ${request.protocol.frequency} Hz, ${request.protocol.pulses} stimuli, ${request.protocol.delay} s recovery delay. Unknown observations do not influence selection. No mutant publication parameter set is used to choose candidates.</p></details>`;
    $("#guideUndo").hidden = false;
  }
  function restore() {
    stop(); if (!base) return;
    applying = true;
    try { hooks.restore(base); } finally { applying = false; }
    result = null; lastApplied = null;
    $("#guideResults").replaceChildren(); $("#guideUndo").hidden = true;
    $("#guideStatus").textContent = "Baseline restored. Your observations are kept.";
  }
  function onSimulation(params) {
    if (!panel || panel.hidden || applying) return;
    const fingerprint = JSON.stringify(params);
    if (fingerprint !== (lastApplied || JSON.stringify(base))) capture(params);
    return result;
  }
  global.PhenotypeGuide = {init, onSimulation};
})(globalThis);
