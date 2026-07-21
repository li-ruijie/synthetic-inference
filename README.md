# Replication code

Code for "Preserving Generalized Linear Model Inference in Synthetic
Data".  This is the same code used in the
development and testing of the paper's results.

The simulation results are reproducible in full from seed 6489.  Some
scripts also carry an application on a restricted dataset, which the
current paper does not include.  That dataset is not publicly available,
and the scripts that read it skip that section when the file is absent,
so everything else runs unchanged.

## Requirements

R 4.6.1.  Three packages are required: MASS, Matrix, and qrng.  Two
further packages are optional, and each skips gracefully when absent:
haven (for the restricted dataset, which is not distributed here) and
ggplot2 (for the diamonds demonstration in
`analysis/25-public-datasets.r`).  The NMES1988 application in
`analysis/73-nmes-twopart.r` reads its data from the AER package, which
that script requires.

Python 3.13.  The Python environment is pinned exactly in
`requirements.txt` (SDV 1.35.1, torch 2.11.0, numpy 2.4.4):

```bash
uv venv --python 3.13 .venv
uv pip install --python .venv/bin/python -r requirements.txt
```

The torch pin is required. CTGAN, TVAE and CopulaGAN train under the
torch random number generator, so the seeding in all four Python generators
(`generate/13-gen-ml-default.py`, `generate/23-gen-ml-tuned.py`, `generate/26-mc-ml-generate.py` and
`generate/72-mc-ml-largen.py`) depends on it.

## Layout

Every path below is relative to the repository root, and every script is run
from there.

```
lib/          sourced by other scripts, not run directly
simulate/     scenario generation and the public DP design constants
generate/     the Python ML generators (SDV, SmartNoise)
analysis/     single-run reports and tables
montecarlo/   the replication drivers
verify/       assertion suites that certify the mechanisms
results/      every CSV, .rds and .RData the scripts read or write
```

## Quick start

Run everything from the repository root.

```bash
# 1. Generate simulation scenarios (saves results/d.sim.scenarios.RData)
Rscript simulate/00-simulate.r

# 1b. Public DP constants from the DGP (needed by analysis/20-dp-precheck.r)
Rscript simulate/64-dp-design-constants.r

# 2. Export scenarios to CSV for Python ML generators
Rscript simulate/12-export-scenarios.r

# 3. Generate ML synthetic data (default hyperparameters)
python generate/13-gen-ml-default.py

# 4. Run single-seed comparison (6 methods, scenarios 1/3/4/6/8/9)
Rscript analysis/17-report-ml-comparison.r

# 5. Monte Carlo all-method comparison (100 reps)
#    Phase 1: ML generation (parallel Python, 2400 tasks).  Measured at
#    28213s (7.8h) on 14 workers, so allow proportionally longer at the
#    default MC_WORKERS=8.
python generate/26-mc-ml-generate.py
#    Phase 2: R diagnostics (parallel R)
Rscript montecarlo/25-mc-allmethod.r

# 6. Differential privacy simulation (mu-GDP releases behind the DP tables).
#    Needs results/dp-design-constants.rds from step 1b.
Rscript montecarlo/62-mc-dp-gdp.r             # Route 1 + Route 2 utility across links and n
Rscript montecarlo/32-mc-dp-projection-zi.r  # zero-inflated release with a private zero-model
Rscript montecarlo/63-dp-audit-mugdp.r       # empirical membership-inference audit of all routes
Rscript montecarlo/58e-mc-mixed-copula.r     # mixed-margin copula categorical fidelity (V4)

# 7. Verify the mechanisms are sound (each must print OK for every task)
Rscript verify/verify-gdp-mechanisms.R   # 12 tasks: the mu-GDP arithmetic
Rscript verify/verify-dp-soundness.R     #  8 tasks: the sensitivity invariants
Rscript verify/verify-projection.R       #  6 tasks: the projection identities
```

## File inventory

### Core libraries

| File | Description |
|------|-------------|
| `lib/00-shared.r` | Shared utilities, environment setup, seed 6489 |
| `lib/gen-syn-family.r` | Gaussian copula generator with Cornish-Fisher quantile approximation |
| `lib/parallel_setup.R` | Portable parallel backend (fork on Linux, PSOCK on Windows) |

### Differential privacy bounds

Every bound used to calibrate DP noise comes from a single named source and carries a
`source` string that is printed into the run log.

| File | Description |
|------|-------------|
| `lib/dp-public-bounds.r` | The only permitted source of a DP bound. Every bound is a functional of the published design, not a statistic of the analysis sample |
| `lib/dp-gaussian-mechanism.r` | The one Gaussian-mechanism calibration. Exact analytic mechanism (Balle & Wang 2018) for all eps, since the classical `sqrt(2 log(1.25/delta))/eps` bound is valid only for eps < 1 |
| `simulate/64-dp-design-constants.r` | Population constants computed from the published DGP (N = 1e6, public seed 20260713), independent of the n = 2000 analysis sample. Run before `analysis/20-dp-precheck.r` |

**Why this matters, and it is the paper's own argument.** A DP noise scale must
be calibrated to a sensitivity bound fixed *before* the data are seen. A bound
read off the sample (`max|y|`, or a min-max rescaling by `range(y)`) gives two
neighbouring datasets two different noise scales, so no sensitivity bound holds
and no DP theorem can be instantiated. That is not a weakened guarantee, it is
the absence of one. `lib/dp-public-bounds.r` exists to make that line impossible to
write by accident, and it *refuses* rather than improvising: where the response
is unbounded above, no finite `B_y` exists and it says so.

### Differential privacy generators

The mu-Gaussian DP releases behind the supplement's DP tables. Every noise scale
is calibrated to a public sensitivity from `lib/dp-public-bounds.r`, not to a
statistic of the analysis sample.

#### Routes and V-numbers are different labels

Both appear throughout the code and they mean different properties. A **Route** names
an assumption regime, meaning what you must be able to state publicly before the
data are seen. A **V-number** labels one generator implementation. The two
cross-cut, so neither can replace the other. Route 2 has no V-number, and V1, V4,
and V5 have no route.

| Label | Function | Route |
|-------|----------|-------|
| V1 | `gen.syn.dp.projected` | none, deprecated, carries no formal DP claim |
| V2 | `gen.syn.dp.ols.public` | Route 1, OLS |
| V3 | `gen.syn.dp.ridge.public` | Route 1, ridge |
| none | `gen.syn.dp.projected.subagg` | Route 2 |
| V4 | `gen.syn.dp.full` | none |
| V5 | `gen.syn.dp.logistic` | none |
| V6 | `gen.syn.dp.suffstats` | Route 3 |

The three sound routes differ in what they require:

- **Route 1** needs a public design and a public bound on `|y|`. It is exact and
  the most efficient of the three.
- **Route 2** needs a public design and a public coefficient box, but no bound on
  the response at all, so it survives an unbounded `y` where Route 1 does not
  exist.
- **Route 3** keeps the design private and needs a public clip on the response.

V4 and V5 release sensitive covariates through a DP copula and carry no route
number. V1 is the analysis-level baseline the paper contrasts the sound routes
against, and it makes no formal DP claim at any epsilon.

| File | Description | Paper reference |
|------|-------------|-----------------|
| `lib/28-dp-projection.r` | Subsample-and-aggregate (Route 2, needs no bound on the response) and the zero-inflated generator with a private zero-model | Supplement S6 |
| `lib/33-dp-projection-public-x.r` | Public-covariate exact release (Route 1), exact refit under a bounded response | Supplement S6 |
| `lib/34-dp-ridge-projection.r` | Ridge coefficient and variance channels under a public design | Supplement S6 |
| `lib/35-dp-full.r` | Sensitive-covariate release (V4), three channels composed in one quadrature. The covariate copula is mixed-margin, so categoricals join it through indicator scores and a calibrated probit/logit back-map with the public biserial de-attenuation (assumption A6) | Supplement S6 |
| `lib/53-dp-logistic.r` | Sensitive-covariate logistic release (V5). Its opt-in copula shares the mixed-margin construction | Supplement S6 |
| `lib/66-dp-suffstat.r` | Whitened sufficient-statistic release (Route 3), vech/X'y/y'y channels in one quadrature | Supplement S6 |

### Generators

| File | Description | Paper reference |
|------|-------------|-----------------|
| `lib/07-calibrated-generator.r` | Empirical Bayes shrinkage calibration | Section 5.1, Theorem S1, Supplement S5 |
| `lib/09-symmetric-generator.r` | Direction-agnostic covariance calibration | Section 5.2, Theorem S2, Supplement S5 |
| `lib/11-projection-generator.r` | Sufficient statistics projection | Section 5.3, Theorem 5, Supplement S1 |
| `generate/13-gen-ml-default.py` | SDV ML generators (CTGAN, TVAE, GaussianCopula, CopulaGAN) | Section 6.2 |

### Simulation and analysis

| File | Description | Paper reference |
|------|-------------|-----------------|
| `simulate/00-simulate.r` | 11 simulation scenarios with known DGP | Section 6.1 |
| `simulate/12-export-scenarios.r` | Export simulation data to CSV for Python | Section 6 |
| `analysis/14-report-inferential.r` | Single-run inferential fidelity on the five criteria | Section 6.3, Supplement S1 |
| `analysis/16-report-privacy.r` | Distance-to-closest-record and membership-inference AUC | Section 7 |
| `analysis/17-report-ml-comparison.r` | Six-method inferential fidelity comparison | Section 6.3, Supplement S2 |
| `analysis/17b-ml-coeftable.r` | Per-coefficient table for all seven methods (race3, age, six-method) | Supplement S2 |
| `montecarlo/18-mc-projection.r` | Projection Monte Carlo, 100 replications | None. The paper reports the 100-replication projection results of `montecarlo/25-mc-allmethod.r` (Supplement S11) |
| `generate/23-gen-ml-tuned.py` | CTGAN with tuned hyperparameters | Supplement S2 |
| `analysis/24-report-tuned.r` | Tuned against default CTGAN and the projection. Run `generate/23-gen-ml-tuned.py` first | Section 6.2, Supplement S2 |
| `analysis/25-public-datasets.r` | Diamonds (ggplot2) and Adult (UCI, downloaded on first run) | Section 6.4, Supplement S3 |
| `montecarlo/59-mc-cross-family.r` | Cross-family refit under full synthesis and public covariates | Supplement S2 |
| `montecarlo/60-mc-binary.r` | Scenario 11 binomial-logit Monte Carlo | Section 6.3 |
| `analysis/60b-binary-coeftable.r` | Scenario 11 per-coefficient table | Supplement S1 |
| `montecarlo/18b-fixed-design-se.r` | Standard errors under full synthesis against the fixed design | Supplement S13 |
| `montecarlo/18c-twopart-se.r` | Two-part projection standard errors, single run | None. Supplement S13 reports the 100-replication run of `montecarlo/18d-twopart-mc.r` |
| `montecarlo/18d-twopart-mc.r` | Two-part projection standard errors, 100 replications | Supplement S13 |
| `montecarlo/19-high-dim-simulation.r` | High-dimensional strata simulation (S proportional to n) | Supplement S10 |
| `analysis/20-dp-precheck.r` | Differential privacy analytical pre-check. Run `simulate/64-dp-design-constants.r` first | Supplement S6 |
| `analysis/21-dp-quantify-existing.r` | Per-person influence (DFBETA). Descriptive only, carries no DP guarantee | Supplement S6 |
| `analysis/22-premasking-demo.r` | Pre-masking privacy improvement | Supplement S6 |
| `montecarlo/68-kwise-fidelity.r` | k-wise distributional fidelity of the projection, coupling bound asserted per draw | Supplement S11 |
| `montecarlo/69-design-gram.r` | Design-Gram projection under full synthesis, Gram identity asserted per draw | Supplement S11 |
| `analysis/71-phi-utilisation.r` | Healthcare-utilisation two-part model application (restricted data, not publicly available, and the script skips if absent) | None. The current paper does not include this application |
| `72-mc-ml-largen-{gen.r, .py, report.r}` | Larger-n ML robustness check (n = 20000, structural failure not a small-sample artefact) | Supplement S11 |
| `analysis/73-nmes-twopart.r` | Two-part model of outpatient visits on NMES1988 (`AER` package), original fit against the copula and two-part projection releases | Section 6.4, Supplement S4 |
| `analysis/75-selection-commutation.r` | All-subset model selection by AIC and BIC on the original data against the release | Supplement S11 |
| `analysis/76-misspec-identity.r` | Misspecified analyst models fitted on the original data and on the release | Supplement S11 |

### Monte Carlo pipeline

| File | Description |
|------|-------------|
| `generate/26-mc-ml-generate.py` | Phase 1: parallel ML synthetic data generation (100 reps) |
| `montecarlo/25-mc-allmethod.r` | Phase 2: R diagnostics for all 7 methods (100 reps) |
| `generate/56-mc-dpml-generate.py` | Phase 1 of the DP machine learning arm: DP-SGD CTGAN and PATE-CTGAN under SmartNoise. Needs a GPU and its own environment of smartnoise-synth, opacus, and Python 3.13, with the recipe given in the script header |
| `analysis/57-mc-dpml-report.r` | Phase 2 of the DP machine learning arm: R diagnostics over the Phase 1 output |

### Differential privacy simulation and verification

The drivers behind the supplement's empirical DP tables, and the three suites
that certify the mechanisms are sound. Run `simulate/64-dp-design-constants.r` first.

| File | Description | Paper reference |
|------|-------------|-----------------|
| `montecarlo/36-mc-all-variants.r` | The single DP-variant driver (V1 to V4, kappa sweep), bounds from `lib/dp-public-bounds.r` | Supplement S6 |
| `montecarlo/42-baseline-comparison.r` | Gaussian mechanism on X'y as the analytic baseline, against the projection routes | Supplement S6 |
| `montecarlo/62-mc-dp-gdp.r` | Route 1 and Route 2 utility across links, at a baseline width, and over n | Supplement S6 |
| `montecarlo/32-mc-dp-projection-zi.r` | Fully private zero-inflated release with a private zero-model | Supplement S6 |
| `montecarlo/63-dp-audit-mugdp.r` | Empirical membership-inference AUC against the mu-GDP bound | Supplement S6 |
| `montecarlo/58e-mc-mixed-copula.r` | Categorical-to-continuous association recovered by the mixed-margin copula, against independent margins | Supplement S6 |
| `montecarlo/66b-mc-dp-suffstat.r` | Route 3 utility across the bounded-response scenarios and the budget grid | Supplement S6 |
| `montecarlo/65-dp-suffstat-survey.r` | Route 3 scale survey (a decision instrument, from which the paper takes only the two S6 condition numbers) | Supplement S6 |
| `verify/verify-gdp-mechanisms.R` | 12 tasks checking the mu-GDP arithmetic (composition, exact dual, calibration) | None |
| `verify/verify-dp-soundness.R` | 8 tasks checking the sensitivity invariants (each fails on an earlier, unsound release) | None |
| `verify/verify-projection.R` | 6 tasks checking the projection identities (design-Gram identity, exact triple, two-part blocks, loud failures, byte-identity) | None |

## Data

The simulation scenarios are generated programmatically
(`simulate/00-simulate.r`) and are reproducible from seed 6489.

Some scripts also carry a healthcare-utilisation application on a
restricted dataset, which the current paper does not include.  That
dataset is **not publicly available**.  It is not distributed in this
repository and cannot be redistributed, and access requests go to the
corresponding author.

Nine scripts read it (`analysis/14-report-inferential.r`,
`analysis/16-report-privacy.r`, `montecarlo/18-mc-projection.r`,
`analysis/20-dp-precheck.r`, `analysis/21-dp-quantify-existing.r`,
`analysis/22-premasking-demo.r`, `montecarlo/68-kwise-fidelity.r`,
`montecarlo/69-design-gram.r`, `analysis/71-phi-utilisation.r`).  Each
tests for the file first and skips that section when it is absent, so
every script runs here without it, and no number in the current paper
depends on it.

The diamonds dataset is publicly available via the `ggplot2` R
package, and the NMES1988 dataset via the `AER` R package.

The Adult (Census Income) dataset is committed at
`results/adult.data.csv` so the repository is self-contained and the run does
not depend on a third-party URL. It is redistributed under its UCI licence,
Creative Commons Attribution 4.0 International (CC by 4.0), which permits
sharing and adaptation for any purpose provided credit is given.

> Becker, B. and Kohavi, R. (1996). *Adult*. UCI Machine Learning Repository.
> <https://doi.org/10.24432/C5XW20>. Derived from 1994 US Census Bureau data.
> Licensed CC by 4.0.

## Reproducibility

All results use master seed 6489.  The Monte Carlo pipeline
(`montecarlo/25-mc-allmethod.r` + `generate/26-mc-ml-generate.py`) generates 100
replications per (scenario, method) combination with deterministic
per-task seeds derived from the master seed.

Seeding the SDV generators takes three calls, and all three are needed.
`np.random.seed` alone reaches neither the torch generator that CTGAN,
TVAE and CopulaGAN train under, nor the fixed internal generator that
`copulas` gives GaussianCopula.  The scripts therefore call
`np.random.seed(seed)` and `torch.manual_seed(seed)` before `fit`, and
`synth._set_random_state(seed)` after it.  Without the last call
GaussianCopula returns the same sample for every seed, so its 100
replications collapse to one draw and its Monte Carlo standard deviation
is zero.

Results reproduce from the master seed on any platform.  Small
floating-point differences arise between platforms, since torch and the
BLAS libraries order their reductions differently, so regenerated
outputs should be compared by value rather than by file hash.

## License

GNU Affero General Public License v3.0 (see `LICENSE`).
Copyright (c) 2026 Ruijie Li.
