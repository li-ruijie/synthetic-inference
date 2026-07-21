#!/usr/bin/env Rscript
# +------------------------------------------------------------------------+
# | 12-export-scenarios.r -- Export simulation data to CSV for Python       |
# +------------------------------------------------------------------------+
#
# Exports the 11 simulation scenarios (from 00-simulate.r) to CSV files
# consumed by the Python ML generators (13-gen-ml-default.py and
# 23-gen-ml-tuned.py). This is the bridge between the R-based DGP/copula
# pipeline and the SDV library [patki-syntheti-2016].
#
# Why this is necessary:
# Section 6 ("Methods compared") compares the copula/projection generators
# against five ML generators: CTGAN [xu-modeling-2019], TVAE, GaussianCopula,
# CopulaGAN, and DataSynthesizer [ping-datasynt-2017a]. These are
# implemented in Python (SDV 1.35, Python 3.13). The comparison tests
# Proposition 2 (no-guarantee): no distributional divergence criterion
# guarantees inferential fidelity. Exporting to CSV ensures all methods
# operate on identical input data.
#
# Input:  d.sim.scenarios.RData (from 00-simulate.r)
# Output: sim_scenario_{01..11}_orig.csv
.shared.loaded <- TRUE
source("lib/00-shared.r")
load("results/d.sim.scenarios.RData")

Map(function(i) {
    d <- l.scenarios[[i]][["data"]]
    s.path <- sprintf("results/sim_scenario_%02d_orig.csv", i)
    write.csv(d, s.path, row.names = FALSE)
    cat(sprintf("Scenario %d: %d rows -> %s\n", i, nrow(d), s.path))
}, seq_along(l.scenarios))
