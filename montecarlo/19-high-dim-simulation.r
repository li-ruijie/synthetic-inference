#!/usr/bin/env r
# +------------------------------------------------------------+
# |  High-dimensional simulation: S proportional to n          |
# +------------------------------------------------------------+
#
# This script addresses the high-dimensional strata discussion,
# which identifies two asymptotic regimes:
#   (a) S = o(n)  -- strata grow slower than sample size;
#       per-stratum counts n_s -> infinity, copula/calibrated
#       generators remain consistent.
#   (b) S proportional to n  -- strata grow with sample size;
#       per-stratum counts stay O(1), stratum-mean estimators
#       do not converge, and copula-based generators fail.
#
# The optimal shrinkage rate (supplement S8) gives the crossover point:
#   n ~ S * sigma^2 / tau^2
# When n/S is small relative to sigma^2/tau^2, the shrinkage
# component of MSE dominates and copula/calibrated generators
# show large coefficient deviations (22--44% in the results).
#
# Theorem 7 (sufficient statistics projection) sidesteps S
# entirely: the projection operates on the p-dimensional
# sufficient statistic (X'X, X'y), not on the S-dimensional
# vector of stratum means.  The number of strata does not enter the
# projection.  Its residual submodel deviations (up to 8.6% in the
# sweep) trace to the bootstrap draw's Gram noise, and the top
# configuration is exact (0.0%), the fitted model there being the
# projected full model.
#
# +------------------------------------------------------------+
# |  Setup                                                     |
# +------------------------------------------------------------+
source("lib/00-shared.r")
source("lib/gen-syn-family.r")
# prevent 00-shared.r from clearing env on re-source from 07/09/11
.sobol.main <- TRUE
.shared.loaded <- TRUE
source("lib/07-calibrated-generator.r")
source("lib/09-symmetric-generator.r")
source("lib/11-projection-generator.r")
rm(.sobol.main, .shared.loaded)

# +------------------------------------------------------------+
# |  Generate high-dimensional data                            |
# +------------------------------------------------------------+

# gen.high.dim: create a dataset with many categorical predictors
# whose cross-product produces S = 233 occupied strata for n = 2000
# observations, giving S/n ~ 0.12.  This places the design in the
# second regime of the high-dimensional strata discussion, an
# unnumbered run-in heading carrying no environment: S proportional
# to n, median n_s ~ 4.
#
# The design uses 7 categorical predictors:
#   5 binary (female, married, smoking, diabetes, hypertension)
#     + 2 three-level (race, edu)  -->  2^5 * 3^2 = 288 possible strata
#   + 1 continuous (age)
# Not all 288 strata are occupied. At seed 6489, S = 233.
gen.high.dim <- function(n = 2000, seed = 6489) {
    set.seed(seed)
    # 5 binary predictors -> 2^5 = 32 combinations
    # + 2 three-level predictors -> 2^5 * 3^2 = 288 possible strata
    # + 1 continuous predictor
    # = 233 occupied strata, median 4 obs each (S/n ~ 0.12)
    v.age    <- round(rnorm(n, 50, 12))
    v.female <- factor(rbinom(n, 1, 0.5))
    v.race   <- factor(sample(1:3, n, replace = TRUE, prob = c(0.6, 0.25, 0.15)))
    v.edu    <- factor(sample(1:3, n, replace = TRUE, prob = c(0.3, 0.4, 0.3)))
    v.married <- factor(rbinom(n, 1, 0.55))
    v.smoking <- factor(rbinom(n, 1, 0.25))
    v.diabetes <- factor(rbinom(n, 1, 0.12))
    v.hypertension <- factor(rbinom(n, 1, 0.30))

    # Response: linear model with true effects from all predictors.
    # sigma = 8 (noise SD), tau^2 varies by predictor (effect sizes
    # range from 2 to 6).  The optimal shrinkage rate crossover at
    # n ~ S*sigma^2/tau^2 = 233*64/~12 ~ 1240 means that for our
    # n=2000 we are near the crossover -- enough for the projection
    # to work well, but copula/calibrated methods already struggle.
    v.eta <- 80 +
        -0.3 * v.age +
        -5 * as.numeric(v.female == "1") +
        3 * as.numeric(v.race == "2") +
        6 * as.numeric(v.race == "3") +
        2 * as.numeric(v.edu == "2") +
        4 * as.numeric(v.edu == "3") +
        -2 * as.numeric(v.married == "1") +
        -3 * as.numeric(v.smoking == "1") +
        -4 * as.numeric(v.diabetes == "1") +
        -3 * as.numeric(v.hypertension == "1") +
        rnorm(n, 0, 8)

    d <- data.frame(
        response = v.eta,
        age = v.age,
        female = v.female,
        race = v.race,
        edu = v.edu,
        married = v.married,
        smoking = v.smoking,
        diabetes = v.diabetes,
        hypertension = v.hypertension)

    # Count strata formed by crossing all categorical predictors.
    # This S is what enters the regime classification of the high-dimensional
    # strata discussion and the optimal shrinkage rate's crossover formula.
    v.strata <- interaction(d[c("female", "race", "edu", "married",
                                "smoking", "diabetes", "hypertension")],
                            drop = TRUE)
    s.S <- length(unique(v.strata))
    v.ns <- table(v.strata)
    cat(sprintf("n=%d, S=%d, S/n=%.2f, median n_s=%d, min n_s=%d\n",
                n, s.S, s.S / n, median(v.ns), min(v.ns)))
    d
}
# +------------------------------------------------------------+
# |  Helpers                                                   |
# +------------------------------------------------------------+
pct <- function(a, b) round(100 * (a - b) / abs(b), 1)
maxpct <- function(a, b) max(abs(pct(a, b)), na.rm = TRUE)
# check: compare synthetic model to original, reporting the same
# five diagnostics used in the MC study (script 18).
# Under Theorem 7, the projection row should show near-zero
# deviations regardless of S/n.  Copula and calibrated rows
# degrade as S/n increases.
check <- function(mo.orig, mo.syn, label) {
    v.o <- coef(mo.orig)
    v.s <- coef(mo.syn)
    s.p <- min(length(v.o), length(v.s))
    v.o <- v.o[seq_len(s.p)]
    v.s <- v.s[seq_len(s.p)]
    s.o <- summary(mo.orig)[["coefficients"]]
    s.s <- summary(mo.syn)[["coefficients"]]
    v.pct <- abs(100 * (v.s - v.o) / abs(v.o))
    v.se.pct <- abs(100 * (s.s[seq_len(s.p), 2] -
                            s.o[seq_len(s.p), 2]) /
                    s.o[seq_len(s.p), 2])
    cat(sprintf("  %-20s  sig=%3d%%  dir=%3d%%  <5%%=%3d%%  max|%%|=%5.1f%%  SE=%5.1f%%\n",
                label,
                round(100 * mean((s.o[seq_len(s.p), 4] < 0.05) ==
                                 (s.s[seq_len(s.p), 4] < 0.05))),
                round(100 * mean(sign(v.o) == sign(v.s))),
                round(100 * mean(v.pct < 5, na.rm = TRUE)),
                max(v.pct, na.rm = TRUE),
                max(v.se.pct, na.rm = TRUE)))
}
# +------------------------------------------------------------+
# |  Run: full model with all predictors                      |
# +------------------------------------------------------------+
#
# First comparison: all 7 categorical predictors active, giving
# the maximum S/n ~ 0.1.  This is the headline result for the
# high-dimensional strata discussion: copula and calibrated
# generators show 22--44% max coefficient deviation while the
# projection is exact (0.0%), the fitted model being the one the
# response was projected for.
cat("=== HIGH-DIMENSIONAL SIMULATION (S proportional to n) ===\n\n")
d <- gen.high.dim()

fml <- response ~ age + female + race + edu + married +
       smoking + diabetes + hypertension
mo.orig <- lm(fml, d)
cat(sprintf("\nOriginal fit: %d coefficients\n", length(coef(mo.orig))))
cat("Significant (p<0.05):",
    sum(summary(mo.orig)[["coefficients"]][, 4] < 0.05),
    "of", nrow(summary(mo.orig)[["coefficients"]]), "\n\n")

# Support: response and age are real-valued.  All other predictors
# are categorical and handled automatically by stratum assignment.
sup <- c(response = "real", age = "real")

# Generate synthetic data from all three generators and compare.
# copula: uses S-dimensional stratum means (vulnerable to high S)
d.cop <- gen.syn.family(d, order = 4, supports = sup)
mo.cop <- lm(fml, d.cop)

# calibrated: applies shrinkage to stratum means (helps, but the
# high-dimensional strata discussion puts the shrinkage
# Contribution to the MSE at O(1) when S/n ~ const)
d.cal <- gen.syn.calibrated(d, order.cf = 4, order.cal = 2,
                            supports = sup)
mo.cal <- lm(fml, d.cal)

# projection: Theorem 7 -- uses p-dimensional sufficient statistic
# (X'X, X'y), not the S-dimensional stratum means.  The number of
# strata S does not enter the projection at all.
d.proj <- gen.syn.projected(d, y.nm = "response",
                            family = "gaussian", link = "identity",
                            supports = sup)
mo.proj <- lm(fml, d.proj)

check(mo.orig, mo.cop, "copula")
check(mo.orig, mo.cal, "calibrated")
check(mo.orig, mo.proj, "projection")

cat("\nCoefficients:\n")
v.o <- coef(mo.orig)
v.c <- coef(mo.cop)
v.l <- coef(mo.cal)
v.p <- coef(mo.proj)
s.p <- min(length(v.o), length(v.c), length(v.l), length(v.p))
d.out <- data.frame(
    orig  = sprintf("%7.2f", v.o[seq_len(s.p)]),
    copul = sprintf("%7.2f (%+.1f%%)", v.c[seq_len(s.p)],
                    pct(v.c[seq_len(s.p)], v.o[seq_len(s.p)])),
    calib = sprintf("%7.2f (%+.1f%%)", v.l[seq_len(s.p)],
                    pct(v.l[seq_len(s.p)], v.o[seq_len(s.p)])),
    proj  = sprintf("%7.2f (%+.1f%%)", v.p[seq_len(s.p)],
                    pct(v.p[seq_len(s.p)], v.o[seq_len(s.p)])))
rownames(d.out) <- names(v.o)[seq_len(s.p)]
print(d.out)

# +------------------------------------------------------------+
# |  Varying S/n ratio                                        |
# +------------------------------------------------------------+
#
# Systematically vary S/n by progressively including more
# categorical predictors in the stratum definition.  This traces
# the transition from the S=o(n) regime to the S ~ n regime
# described in the high-dimensional strata discussion.
#
# The sweep runs S/n from 0.001 to 0.117 (the table below and the
# paper's supplement):
#   - Copula max|%| runs 24.2--39.7% and worsens as S/n grows
#   - Calibrated max|%| runs 22.0--43.6%, worst at the highest S/n
#   - Projection stays below 9% on every submodel row and is exact
#     at the highest S/n, where the fitted model is the projected full
#     model (score pinned, so the refit is exact on any draw)
#
# This directly illustrates Theorem 7's advantage: by projecting
# onto sufficient statistics of dimension p (fixed), the generator
# is invariant to the number of strata S.
cat("\n\n=== VARYING S/n RATIO ===\n\n")
cat(sprintf("%-10s  %-10s  %8s  %8s  %8s\n",
            "S/n", "median_ns", "copula", "calibr", "project"))
cat(strrep("-", 60), "\n")

# Each configuration adds one more categorical predictor to the
# stratum definition.  S counts the occupied cells, since the code
# takes interaction(..., drop = TRUE), so it tracks the nominal
# product only while every cell is filled and falls behind it once
# the strata outrun the sample.  The S/n sequence is:
#   "1 bin"      -> S/n=0.001
#   "+race"      -> S/n=0.003
#   "+edu"       -> S/n=0.009
#   "+married"   -> S/n=0.018
#   "+smoking"   -> S/n=0.036
#   "+diabetes"  -> S/n=0.066  (approaching crossover)
#   "+hypert"    -> S/n=0.117  (past crossover)
l.configs <- list(
    list(cats = c("female"), label = "1 bin"),
    list(cats = c("female", "race"), label = "1+3lev"),
    list(cats = c("female", "race", "edu"), label = "+edu"),
    list(cats = c("female", "race", "edu", "married"), label = "+married"),
    list(cats = c("female", "race", "edu", "married", "smoking"),
         label = "+smoking"),
    list(cats = c("female", "race", "edu", "married", "smoking",
                  "diabetes"), label = "+diabetes"),
    list(cats = c("female", "race", "edu", "married", "smoking",
                  "diabetes", "hypertension"), label = "+hypert"))

# For each configuration, count strata, generate synthetic data from
# all three generators, and report max|%| coefficient deviation.
# The copula and calibrated columns degrade with S/n (as predicted
# by the high-dimensional strata discussion and the optimal
# shrinkage rate, neither of them numbered).  The projection stays flat
# (as predicted by Theorem 7).
Map(function(cfg) {
    v.strata <- interaction(d[cfg[["cats"]]], drop = TRUE)
    s.S <- length(unique(v.strata))
    v.ns <- table(v.strata)

    # fit original
    fml.sub <- reformulate(c("age", cfg[["cats"]]), response = "response")
    mo.o <- lm(fml.sub, d)

    sup.sub <- c(response = "real", age = "real")
    d.c <- gen.syn.family(d, order = 4, supports = sup.sub)
    d.l <- gen.syn.calibrated(d, order.cf = 4, order.cal = 2,
                              supports = sup.sub)
    d.p <- gen.syn.projected(d, y.nm = "response",
                             family = "gaussian", link = "identity",
                             supports = sup.sub)
    mo.c <- lm(fml.sub, d.c)
    mo.l <- lm(fml.sub, d.l)
    mo.p <- lm(fml.sub, d.p)

    s.pc <- min(length(coef(mo.o)), length(coef(mo.c)))
    s.pl <- min(length(coef(mo.o)), length(coef(mo.l)))
    s.pp <- min(length(coef(mo.o)), length(coef(mo.p)))

    cat(sprintf("%-10s  ns=%4d     %6.1f%%   %6.1f%%   %6.1f%%\n",
                sprintf("%.3f", s.S / nrow(d)),
                median(v.ns),
                maxpct(coef(mo.c)[seq_len(s.pc)], coef(mo.o)[seq_len(s.pc)]),
                maxpct(coef(mo.l)[seq_len(s.pl)], coef(mo.o)[seq_len(s.pl)]),
                maxpct(coef(mo.p)[seq_len(s.pp)], coef(mo.o)[seq_len(s.pp)])))
}, l.configs)
# +------------------------------------------------------------+
