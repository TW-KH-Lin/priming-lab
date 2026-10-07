# Phenotype guide (v45)

The optional guide accepts observations rather than kinetic-rate inputs. Activate
“Describe my phenotype”, enter condition/control ratios if measured, choose
changes, check the frequency, train length and probe delay, then run one comparison.

Unknown pFusion and total release-site count N are held at the current baseline.
N is not the readily releasable pool. Unknown outcomes are omitted, whereas
“unchanged” is an explicit constraint. “Measured ratio” constrains an outcome to
within 5%; directional selections require at least 5% change and do not specify
its magnitude. The comparison reports each outcome separately, without a global
match percentage. Recovery is recovered fraction at one specified probe delay,
not a fitted recovery time constant. Steady state is the mean of the last five
pulses (or all pulses if fewer), not proof of asymptotic convergence.

Implementation: `dist/phenotype-engine.js` forward-simulates a finite catalog of
changes to k1, b1, k2, b2 and the two calcium slopes. Other mechanisms, including
b4 and calcium waveforms, remain at baseline. If first EPSC is constrained,
resting equilibrium determines k2 or b2 for each candidate; impossible candidates
are discarded. No global optimizer or mutant publication parameters are used.
Long protocols use a deterministic subset of 64 candidates. A worker keeps the
interface responsive and can be cancelled. The returned candidate minimizes
squared outcome residuals within this catalog, with a tiny preference for smaller
parameter changes when outcomes are equivalent. It is not a unique identification
of biological rate constants. Finite candidate coverage can leave targets unmet.

The captured baseline is immutable between comparisons. Main-model changes start
a fresh baseline. Restore or deactivate returns to that baseline. The simulation
JSON export contains observations, protocol, baseline, parameters, candidate count
and per-outcome residuals for reproducibility.

Checks: `node tests/phenotype-guide.cjs`, existing golden/protocol tests, and
`node tests/phenotype-ui.cjs` with jsdom available on NODE_PATH. The latter exercises
the page DOM and worker interactions in memory, not real-browser layout.

Synthetic-data recovery tests demonstrate internal consistency only. These tests
do not establish reproduction of any published experiment. Existing publication
presets and the simulation core are unchanged in this version.
