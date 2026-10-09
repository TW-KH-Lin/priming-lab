# Depression and recovery phenotype patch — v49

This focused patch leaves `dist/sim-core.js` and the forward-simulation equations unchanged. It extends the generic mechanism inference layer with independent train-depth, train-speed and early-recovery readouts.

## Train metrics

For each simulated frequency the mechanism result exports:

- `ssOverP1`: mean of the final five normalized responses;
- `t50Dep` and `n50Dep`: first threshold crossing halfway from the train peak to the final plateau;
- `tSS` and `nSS`: first point after which every remaining response stays inside the configured plateau band, with at least three responses available.

Time is measured from pulse 1. The train peak is used when facilitation precedes depression. If a stable plateau is not reached, `tSS`/`nSS` remain null and the UI reports “Not reached.” The default plateau tolerance is 10% of total depression amplitude and is editable from 1–50%.

## Recovery metrics

Recovery is summarized over 0–4 seconds using the points 0.25, 0.5, 1, 2 and 4 seconds:

- trapezoidal `auc0to4`, anchored at zero recovery at time zero;
- linearly interpolated `t50` and `t80` threshold times;
- a censoring flag when a threshold is not reached by 4 seconds.

Censored times are displayed as `>4 s`; they are never extrapolated. Points at 8 and 16 seconds remain available for visualization but carry no fitting weight unless a user explicitly enters them as targets.

## Block-normalized scoring

Targets are grouped into release magnitude, PPR, train depth, train kinetics, recovery kinetics, full-train shape and calcium blocks. Residuals are averaged within each block before the configurable total block weight is applied. The full normalized train defaults to weight 0.25. Consequently, a 40-pulse curve refines shape without counting forty times more than one PPR measurement.

`tests/mechanism-phenotype-metrics.cjs` covers independent depth/speed, non-stable plateaus, early recovery despite equal 16-second endpoints, recovery threshold censoring and 40-point block normalization.
