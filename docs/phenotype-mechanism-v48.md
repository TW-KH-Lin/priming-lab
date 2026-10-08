# Phenotype → Mechanism — v48 architecture cleanup

The production `/mechanism.html` page is publication-agnostic. It accepts a control model, hard-locked measured pFusion/N_total values, phenotype targets and uncertainty, then repeatedly calls the unchanged `dist/sim-core.js` forward simulator through a cancellable worker. FRP = LS + TS and TS are phenotype targets; neither is silently assigned to N_total.

The generic engine searches priming, actual-versus-reference calcium, AP-calcium amplitude/decay, and optional STP/refractory families. It reports minimal and best-numerical candidates separately, alternatives, identifiability, necessity, sufficiency and phenotype attribution. It does not classify a phenotype by a publication name. Treatment presets remain available only in the forward simulator; the mechanism reference menu excludes them.

Publication retrospective validation is isolated under `tests/mechanism_validation/`. Every study has an `*_input.json` containing control/constraints/experimental phenotype and a separate `*_answer.json` containing expected directions/published values. `tests/mechanism-blind.cjs` runs inference, writes a frozen prediction, and only then opens the answer. Neither file type is imported by the production HTML, UI, engine or worker, and answer files are absent from `dist/`.

Validation types must remain distinct:

- `tests/mechanism-blind.cjs`: experimental retrospective validation.
- `tests/mechanism_validation/synthetic_hk.cjs`: synthetic inverse-model regression (published parameters generate same-engine targets); it is not independent biological validation.

Run:

```text
node tests/mechanism.cjs
node tests/mechanism-validation-architecture.cjs
node tests/mechanism-blind.cjs extca|iono|pdbu|hk
node tests/mechanism-ui.cjs  # with jsdom available
```

Set `ASSERT_SCIENCE=1` for the retrospective mechanism-direction gate. A high curve-fit score is not sufficient for that gate. Current experimental summary targets do not uniquely recover every published mechanism; a failed scientific gate should remain visible rather than being converted into a production-page PASS.
