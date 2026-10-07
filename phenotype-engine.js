/* Finite mechanism comparison. No publication mutant parameters or global fitting. */
(function (global) {
  "use strict";
  const C = global.PrimingCore;
  const adjustable = ["kf1", "kb1", "kf2", "kb2", "slope1", "slope2"];
  const labels = {first:"First EPSC",ppr:"Paired-pulse ratio",depression:"Depression at train end",steady:"Steady-state release",recovery:"Recovery at probe delay"};
  function validate(base, request) {
    if (!(base.pRel0 > 0 && base.pRel0 <= 1 && base.nSites > 0)) throw Error("Choose a baseline with positive pFusion and total sites before comparing ratios.");
    const {frequency, pulses, delay} = request.protocol;
    if (!(frequency >= .5 && frequency <= 333 && Number.isInteger(pulses) && pulses >= 2 && pulses <= 200 && delay >= .01 && delay <= 16)) throw Error("Use 0.5–333 Hz, 2–200 stimuli and a 0.01–16 s probe delay.");
    for (const [key, value] of Object.entries(request.measured)) {
      if (value !== null && !(Number.isFinite(value) && value > 0 && value <= 20)) throw Error("Measured ratios must be positive numbers, up to 20×.");
    }
    if (base.pRel0 * (request.measured.pFusion ?? 1) > 1) throw Error("This pFusion ratio exceeds probability 1. Check the baseline and measurement.");
    for (const target of Object.values(request.targets)) {
      if (!["none", "same", "up", "down", "ratio"].includes(target.mode)) throw Error("Invalid observation.");
      if (target.mode === "ratio" && !(Number.isFinite(target.value) && target.value > 0 && target.value <= 20)) throw Error("Enter a positive target ratio, up to 20×.");
    }
    if (!Object.values(request.targets).some(t => t.mode !== "none")) throw Error("Choose at least one observed change, or enter a first EPSC ratio.");
  }
  function metrics(p, protocol, needRecovery = true) {
    const train = C.simulateTrain(protocol.frequency, protocol.pulses, p);
    const rows = train.rows, tail = rows.slice(-Math.min(5, rows.length));
    const first = rows[0].release, steady = tail.reduce((sum, row) => sum + row.release, 0) / tail.length;
    const recovery = needRecovery ? C.simulateRecovery(train, [protocol.delay], p, "last5").points[0].recovered : null;
    return {first, ppr: rows[1].release / first, steady, depression: 1 - steady / first, recovery};
  }
  function firstAtRest(p) { return C.steadyState(p).ts * p.pRel0; }
  function constrainFirst(p, release, via = "kf2") {
    const above = Math.max(p.caRestActual - p.caRestReference, 0);
    let k1 = p.kf1 + p.slope1 * above;
    if (p.useMM) k1 /= 1 + above / Math.max(p.kHalf1, 1e-18);
    const fraction = release / (p.nSites * p.pRel0);
    if (!(fraction > 0 && fraction < 1 && k1 > 0)) return false;
    if (via === "kb2") p.kb2 = (p.kf2 + p.slope2 * above) * k1 * (1 - fraction) / (fraction * (p.kb1 + k1));
    else p.kf2 = fraction * p.kb2 * (p.kb1 + k1) / (k1 * (1 - fraction)) - p.slope2 * above;
    return adjustable.every(key => Number.isFinite(p[key]) && p[key] >= 0) && p.kb2 > 0;
  }
  function candidateList(base, request) {
    const fixed = {...base, pRel0: base.pRel0 * (request.measured.pFusion ?? 1), nSites: base.nSites * (request.measured.pool ?? 1)};
    const templates = [{}];
    const factors = [.25, .5, 1, 2, 4];
    for (const key of adjustable) for (const f of factors) if (f !== 1) templates.push({[key]:f});
    for (const a of [.125,.25,.5,1,2,4]) for (const b of factors) templates.push({slope1:a,slope2:b});
    for (const a of factors) for (const b of factors) templates.push({kf1:a,kb1:a,kf2:b,kb2:b});
    for (const a of [.125,.25,.5,1,2,4]) for (const b of factors) templates.push({slope1:a,kf2:b,kb2:b});
    const firstTarget = request.targets.first;
    const ratio = firstTarget?.mode === "ratio" ? firstTarget.value : firstTarget?.mode === "same" ? 1 : null;
    const targetRelease = ratio === null ? null : firstAtRest(base) * ratio;
    const candidates = [], seen = new Set();
    for (const factors of templates) {
      for (const via of targetRelease === null ? ["kf2"] : ["kf2", "kb2"]) {
        const p = {...fixed};
        for (const [key, factor] of Object.entries(factors)) p[key] *= factor;
        if (targetRelease !== null && !constrainFirst(p, targetRelease, via)) continue;
        if (!adjustable.every(key => base[key] === 0 ? p[key] === 0 : p[key] / base[key] >= 1 / 32 && p[key] / base[key] <= 32)) continue;
        const id = adjustable.map(key => p[key].toPrecision(12)).join(":");
        if (!seen.has(id)) { seen.add(id); candidates.push(p); }
      }
    }
    // Slow protocols remain cancellable; test a reproducible subset to limit work.
    const duration = request.protocol.pulses / request.protocol.frequency + request.protocol.delay;
    if (duration > 20 && candidates.length > 64) return Array.from({length:64}, (_,i) => candidates[Math.round(i * (candidates.length - 1) / 63)]);
    return candidates;
  }
  function assess(before, after, targets) {
    const rows = Object.entries(targets).filter(([,t]) => t.mode !== "none").map(([key,t]) => {
      const a = before[key], b = after[key];
      const valid = Number.isFinite(a) && Number.isFinite(b) && a > 1e-10;
      const ratio = valid ? b / a : null;
      const target = t.mode === "ratio" ? t.value : 1;
      let error = Infinity, matched = false;
      if (valid) {
        error = t.mode === "up" ? Math.max(0, 1.05 - ratio) : t.mode === "down" ? Math.max(0, ratio - .95) : Math.abs(ratio - target) / target;
        matched = t.mode === "up" || t.mode === "down" ? error < 1e-9 : error <= .05;
      }
      return {key, label:labels[key], mode:t.mode, target:t.mode === "ratio" ? t.value : null, before:a, after:b, ratio, matched, error:Number.isFinite(error) ? error : null};
    });
    return {rows, loss:rows.reduce((s,r) => s + (r.error === null ? 1e6 : r.error ** 2), 0), matched:rows.filter(r => r.matched).length};
  }
  function compare(base, request, progress = () => {}) {
    validate(base, request);
    const needRecovery = request.targets.recovery?.mode !== "none";
    const baseline = metrics(base, request.protocol, needRecovery), candidates = candidateList(base, request);
    if (!candidates.length) throw Error("No candidate satisfies the first EPSC and pool constraints within the tested range. Check those observations.");
    let best = null;
    candidates.forEach((params, index) => {
      const values = metrics(params, request.protocol, needRecovery), assessment = assess(baseline, values, request.targets);
      const distance = adjustable.reduce((sum,key) => sum + Math.abs(Math.log((params[key] + 1e-20) / (base[key] + 1e-20))), 0);
      const rank = assessment.loss + 1e-7 * distance;
      if (!best || rank < best.rank) best = {params, values, assessment, rank};
      if (index % 5 === 0 || index === candidates.length - 1) progress(index + 1, candidates.length);
    });
    return {baselineParams:base, baseline, ...best, request, tested:candidates.length, method:"finite-mechanism-candidates-v1"};
  }
  global.PhenotypeEngine = {compare, metrics, assess, candidateList, constrainFirst, validate, adjustable};
})(globalThis);
