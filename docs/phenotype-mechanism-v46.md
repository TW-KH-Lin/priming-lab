# Phenotype → Mechanism — v46

## Scope

Separate `/mechanism.html` module. The v45 forward simulator and `dist/sim-core.js` are unchanged. No AI service, account, upload or paid API is used. Computation is in a cancellable browser worker.

Core SHA-256: `485f6314ea0acca8a48545e82dd8f3a4890d041333054173d9e4164e8284939e`.

## Search and interpretation

- Six default kinetic parameters; optional TSL, ERS, calcium-decay and y/z time constants. Other parameters stay at the reference values. Measured N and pFusion are immutable.
- N is total release sites, **not** the measured FRP/RRP. Do not enter FRP as N without an explicit mapping assumption.
- Qualitative ranges produce zero loss inside the range. Missing observations have zero weight. Quantitative errors use supplied SD/SEM, a normal-approximation 95% CI scale, or an explicitly labelled 10% scale.
- Objective: mean weighted squared standardized residuals plus log-deviation and changed-count penalties. A 1% relative threshold defines the displayed changed count; exact values are exported.
- Single-parameter screening precedes subset optimization. Bounded Nelder–Mead search uses deterministic starts/seed. Fast checks singles and selected pairs; Standard stages larger subsets when needed; Exhaustive searches every subset up to the user-selected maximum, not all continuous parameter values.
- Default minimal model is the fewest changed parameters within 10% of the best tested error. An optional 90%-coverage rule is less strict and may miss individual targets. The UI explicitly reports unmet observations.
- Coverage = clamped relative error reduction from WT. It is not a fraction of biological mechanism explained. Undefined when WT error is zero.
- Parameter sufficiency measures improvement beyond WT **plus measured locks**. Necessity compares the candidate to its single-parameter ablation, normalized to its difference from WT plus locks. These scores can be conditional, non-additive and non-unique.
- Alternative candidates are labelled either within the near-optimal error tolerance, or outside it but within every observation's target range/uncertainty scale.
- Sampled 5–95% parameter ranges and confidence labels are search-conditional descriptive measures, not statistical confidence intervals. Fewer than five distinct near-optimal samples yield “Insufficient sampling”.
- Recovery uses the unchanged core's true post-conditioning state, including pools, calcium, y and z. Probe branches are non-destructive.
- Outputs include traces, resting occupancies, parameter attribution, necessity/sufficiency, Pareto envelope, alternatives and a simulated discriminating experiment. JSON and CSV retain settings, seed, bounds, targets, parameters and traces.

## HK validation

`tests/mechanism-hk.cjs` constructs target readouts from the independently transcribed Science Advances supplementary Table S2 HK parameters, then gives inference only WT, readouts and fixed N/pFusion. The optimizer does not receive HK kinetic rates. This is a **synthetic same-engine parameter-recovery test**, not an independent fit to experimental observations.

The legacy forward-app HK preset differs from Table S2. It is not used as validation truth and is omitted from the new module's reference menu. It remains unchanged in the original simulator to preserve that version. The verified comparison uses HK N = 3150, pFusion = 0.379, k1 = 0.419, b1 = 0.194, sigma1 = 221000, k2 = 0.2916, b2 = 0.3514, sigma2 = 1810000.

## Checks and limitations

Run `node tests/mechanism.cjs` and, with jsdom available, `node tests/mechanism-ui.cjs`. Original forward golden/protocol and v45 phenotype tests are retained. `REPORT_PATH=/tmp/hk.json node tests/mechanism-hk.cjs` regenerates the full benchmark.

Numerical and in-memory DOM tests do not substitute for real device visual QA. Browser-control verification was unavailable during this task; live clicking and phone rendering have not been verified. Long low-frequency protocols and exhaustive searches may be slow on phones; cancellation terminates the worker.

The optional recovery-tau fitting and general conditioning-target input are not included in this first module. Additional experiment suggestions currently compare unmeasured train frequencies and recovery intervals, not arbitrary protocols. The requested independent J Physiol phenotype-first comparison is paused at the user's instruction; no validation claim is made for those conditions.
