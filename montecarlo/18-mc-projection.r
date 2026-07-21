#!/usr/bin/env r
# +------------------------------------------------------------+
# |  Monte Carlo replications of the projection generator      |
# +------------------------------------------------------------+
#
# Section 6 of the paper describes the Monte Carlo design:
#   "Monte Carlo summaries use 100 independent replications,
#    reproducible from master seed 6489."
#
# Purpose: empirically verify Theorem 7 (sufficient statistics
# projection).  For continuous GLMs the projection generator
# reproduces the original sufficient statistics exactly, so every
# replication should return identical coefficient estimates
# (sd = 0.0 across 100 runs).  The SE column is a different matter.
# This script synthesises the design as well, so the refit
# information is formed on the synthetic design and the standard
# errors move on every row, by 1.6 to 6.9% on this grid and 19.2%
# on PHI (18b-fixed-design-se.r separates the two regimes).  The
# only non-zero coefficient deviation arises from Poisson/count
# responses where probabilistic rounding injects randomness
# (~2.4% max coefficient deviation, sd ~1.1%).
#
# Expected results (this projection-only MC, master seed 6489):
#   Scenario 1  gaussian(identity)   0.0 (0.0)
#   Scenario 3  gamma(log)           0.0 (0.0)
#   Scenario 6  poisson(log)         2.4 (1.1)
# Continuous families achieve exact zero in the coefficient column
# because the projection pins the response-side sufficient statistic
# that determines the MLE, X'y under a canonical link and the
# weighted score identity otherwise, leaving no stochastic component
# in the estimate.
#
# +------------------------------------------------------------+
# |  Setup                                                     |
# +------------------------------------------------------------+
.sobol.main <- TRUE
source("lib/00-shared.r")
s.wd <- tryCatch(
    normalizePath(dirname(sys.frame(1)[["ofile"]]), mustWork = TRUE),
    error = function(e) getwd())
setwd(s.wd)
source("lib/gen-syn-family.r")
source("lib/07-calibrated-generator.r")
source("lib/09-symmetric-generator.r")
source("lib/11-projection-generator.r")
rm(.sobol.main)

# +------------------------------------------------------------+
# |  Helpers                                                   |
# +------------------------------------------------------------+

# fit.model: fit the analysis model on either original or synthetic
# data.  The family/link combinations mirror the simulation scenarios:
# gaussian(identity), gamma(log), quasi-gamma(log), poisson(log).
# For gamma with zeros in the response, rows with y=0 are dropped
# and a quasi(log, mu^2) family is used with explicit start values
# to aid convergence (quasi-gamma approximation).
fit.model <- function(d, family, link) {
    v.pred <- intersect(c("age", "female", "race", "income",
                           "edu", "housing"), colnames(d))
    v.resp <- setdiff(names(d), v.pred)[1]
    fml <- reformulate(v.pred, response = v.resp)
    tryCatch(
        `if`(family == "gaussian", lm(fml, d),
        `if`(family == "poisson", glm(fml, d, family = poisson()), {
            d2 <- `if`(link == "log" && any(d[[v.resp]] == 0),
                       d[d[[v.resp]] > 0, ], d)
            s.np <- ncol(model.matrix(fml, d2))
            glm(fml, d2, family = quasi(link = "log", variance = "mu^2"),
                start = c(1, rep(0, s.np - 1)),
                control = glm.control(maxit = 200))
        })),
        error = function(e) NULL)
}

# mc.one: run a single MC replication of the projection generator.
# Generates one synthetic dataset via gen.syn.projected, fits the
# analysis model, and returns diagnostics comparing synthetic to
# original:
#   maxpct -- max absolute % deviation across all coefficients
#   se     -- max absolute % deviation in standard errors
#   sig    -- proportion of coefficients with matching significance
#   dir    -- proportion of coefficients with matching sign
#   lt5    -- proportion of coefficients with <5% deviation
#
# Under Theorem 7, for continuous GLMs the coefficient metrics
# should be exactly zero (maxpct=0, sig=100%, dir=100%, lt5=100%)
# across all 100 seeds, because the sufficient statistics are
# preserved exactly.  se is not among them.  The theorem holds the
# design fixed and synthesises only the response, while
# gen.syn.projected draws X* from the copula too, so the refit
# information is formed on the synthetic design and the standard
# errors move, by 1.6 to 6.9% on this grid and 19.2% on PHI.
# 18b-fixed-design-se.r separates the two regimes.  For Poisson,
# probabilistic rounding makes maxpct non-zero as well.
mc.one <- function(d, y.nm, family, link, supports, seed, mo.orig) {
    d.p <- gen.syn.projected(d, y.nm = y.nm, family = family,
                             link = link, supports = supports,
                             seed = seed)
    mo.p <- fit.model(d.p, family, link)
    if(is.null(mo.p)) return(NULL)

    # Align every comparison by coefficient name.  A bootstrapped design
    # can lose a rare factor level (rank deficiency), in which case the
    # synthetic fit carries an NA coefficient and summary() drops the
    # aliased row, so positional indexing would misalign the remaining
    # rows.  Name indexing keeps the estimable coefficients paired and
    # leaves NA for the lost level, which na.rm then excludes.
    v.o <- coef(mo.orig)
    v.p <- coef(mo.p)[names(v.o)]

    s.o <- summary(mo.orig)[["coefficients"]]
    s.p.s <- summary(mo.p)[["coefficients"]]
    v.se.o <- s.o[, 2][names(v.o)]
    v.se.p <- s.p.s[, 2][names(v.o)]
    v.p4.o <- s.o[, 4][names(v.o)]
    v.p4.p <- s.p.s[, 4][names(v.o)]

    v.pct <- abs(100 * (v.p - v.o) / abs(v.o))
    v.se.pct <- abs(100 * (v.se.p - v.se.o) / v.se.o)

    list(maxpct = max(v.pct, na.rm = TRUE),
         se     = max(v.se.pct, na.rm = TRUE),
         sig    = mean((v.p4.o < 0.05) == (v.p4.p < 0.05), na.rm = TRUE),
         dir    = mean(sign(v.o) == sign(v.p), na.rm = TRUE),
         lt5    = mean(v.pct < 5, na.rm = TRUE),
         fullrank = !anyNA(v.p))
}
# +------------------------------------------------------------+
# |  Run Monte Carlo                                          |
# +------------------------------------------------------------+

# Load the pre-generated simulation scenarios (created by the
# scenario setup script).  Each scenario encodes a (family, link)
# combination with a simulated dataset of n=2000.
load("results/d.sim.scenarios.RData")

# B=100 replications, as stated in Section 6.  This is sufficient
# to estimate mean and SD of the deviation metrics.  For continuous
# GLMs we expect sd=0 exactly, so even B=2 would suffice -- 100
# replications are used to also characterise the Poisson rounding
# distribution.
s.B <- 100

# Support specifications for each scenario's response and continuous
# predictors.  These tell the projection generator the domain of each
# variable (real, positive, count, zeroinfl) so it can apply the
# correct rounding/clamping after projection.  NULL means all
# variables are categorical (no continuous support needed).
l.sup <- list(
    c(response = "real",     age = "real", income = "real"),
    NULL,
    c(response = "positive", age = "real", income = "positive"),
    c(response = "zeroinfl", age = "real", income = "positive"),
    NULL,
    c(response = "count",    age = "real", income = "positive"),
    c(response = "positive", age = "real", income = "positive"),
    c(response = "zeroinfl", age = "real", income = "positive"),
    NULL, NULL)

# Subset of scenarios tested in the MC study.  These correspond to
# the rows of the results table in Section 6.  Scenarios 7 and 8
# follow the current numbering, where 7 is a clean Gamma(log)
# response with correlated predictors, 8 is the all-combined
# quasi-gamma design with covariate-dependent zero-inflation.
#   1 = gaussian(identity)   -- Theorem 7 exact (sd=0)
#   3 = gamma(log)           -- exact (sd=0), non-canonical route
#   4 = quasi-gamma(log)     -- zero-inflated gamma variant
#   6 = poisson(log)         -- rounding error ~2.4% (1.1%)
#   7 = gamma(log)           -- correlated predictors, clean Gamma
#   8 = quasi-gamma(log)     -- all combined (ZI + ceiling)
#   9 = gaussian(identity)   -- bounded response [0, 100]
# Theorem 7 is stated for a canonical link and the Gamma canonical
# link is the inverse, so the log-link rows owe their exactness to the
# Newton iteration on the log scale (the numerical implementation on
# positive support), which solves the score equation for any link,
# rather than to the theorem's X'y* = X'y condition.
v.test <- c(1, 3, 4, 6, 7, 8, 9)
l.fam  <- c("gaussian", "gamma", "quasi-gamma",
             "poisson", "gamma", "quasi-gamma", "gaussian")
l.lnk  <- c("identity", "log", "log", "log", "log", "log", "identity")

cat(sprintf("MONTE CARLO PROJECTION RESULTS (%d replications)\n", s.B))
cat(strrep("=", 80), "\n\n")

cat(sprintf("%-28s  %10s  %10s  %10s  %10s\n",
            "", "max|%|", "SE max", "sig%", "<5%"))
cat(strrep("-", 80), "\n")

# Main MC loop: for each scenario, fit the original model once,
# then generate B=100 synthetic datasets with different seeds and
# compare each to the original.  Seeds are 6489+1, ..., 6489+100,
# derived from the master seed 6489 (Section 6).
#
# For continuous GLMs the projection matches the response-side
# sufficient statistics Theorem 7 asks for, X'y* = X'y and the
# saturated kernel under a canonical link, and the weighted score
# identity under a non-canonical one.  Therefore the MLE from every
# synthetic dataset is identical to the original MLE, yielding sd=0
# across all 100 replications.  The information is a separate object,
# formed on the synthetic design, and is not preserved, so the SE
# column is non-zero.
Map(function(i, fam, lnk) {
    s   <- l.scenarios[[i]]
    d   <- s[["data"]]
    sup <- l.sup[[i]]
    mo.o <- fit.model(d, s[["family"]], s[["link"]])
    if(is.null(mo.o)) return(NULL)

    v.seeds <- 6489 + seq_len(s.B)
    l.res <- Map(function(seed) mc.one(d, "response", fam, lnk, sup,
                                       seed, mo.o),
                 v.seeds)
    l.res <- Filter(Negate(is.null), l.res)

    v.maxpct <- mapply(function(r) r[["maxpct"]], l.res)
    v.se     <- mapply(function(r) r[["se"]], l.res)
    v.sig    <- mapply(function(r) r[["sig"]], l.res)
    v.lt5    <- mapply(function(r) r[["lt5"]], l.res)

    cat(sprintf("Scenario %d %-17s  %4.1f (%.1f)  %4.1f (%.1f)  %5.1f%%       %5.1f%%\n",
                i, sprintf("%s(%s)", s[["family"]], s[["link"]]),
                mean(v.maxpct), sd(v.maxpct),
                mean(v.se), sd(v.se),
                100 * mean(v.sig),
                100 * mean(v.lt5)))
}, v.test, l.fam, l.lnk)

# +------------------------------------------------------------+
# |  Real data: PHI waist circumference                       |
# +------------------------------------------------------------+
#
# Apply the same MC protocol to the real PHI baseline dataset.
# The analysis model is OLS (gaussian, identity) for waist
# circumference, so Theorem 7 guarantees exact reproduction of
# the estimable coefficients.  A row-bootstrap of the n=1927
# design is rank-deficient in about 37% of draws ((1-1/n)^n:
# the lone observation in the rarest housing level is not
# resampled).  project.ols absorbs the deficiency via QR, so
# those draws still reproduce the estimable coefficients exactly
# but lose the missing level.  The full-rank count below records
# the empirical rate.  This confirms on real data that the
# sufficient statistics projection holds exactly in finite
# samples wherever the design retains rank.
cat("\n")
# The PHI dataset is restricted.  Everything above runs without it.  The
# application section below needs it, so stop cleanly when it is absent.
if (!file.exists("phi_phase1_baseline.dta")) {
    cat("\nPHI data not available (restricted). Skipping the application section.\n")
    quit(save = "no", status = 0)
}

library(haven)
d.raw <- read_dta("phi_phase1_baseline.dta")
d.app <- data.frame(
    waist   = as.numeric(d.raw[["waist_circum"]]),
    age     = 2016 - as.numeric(d.raw[["birthyear"]]),
    female  = factor(as.numeric(d.raw[["female"]])),
    race    = factor(as.numeric(d.raw[["race"]])),
    edu     = factor(findInterval(as.numeric(d.raw[["education"]]),
                                  c(5, 8)) + 1),
    housing = factor(as.numeric(d.raw[["housingtype"]])))
d.app <- d.app[complete.cases(d.app), ]

sup.app <- c(waist = "real", age = "real")
fml.app <- waist ~ age + female + race + edu + housing
mo.o.app <- lm(fml.app, d.app)

# Same seed sequence as above: 6490, ..., 6589.
v.seeds <- 6489 + seq_len(s.B)
l.res.app <- Map(function(seed) mc.one(d.app, "waist", "gaussian",
                                       "identity", sup.app, seed,
                                       mo.o.app),
                 v.seeds)
l.res.app <- Filter(Negate(is.null), l.res.app)

v.maxpct <- mapply(function(r) r[["maxpct"]], l.res.app)
v.se     <- mapply(function(r) r[["se"]], l.res.app)
v.sig    <- mapply(function(r) r[["sig"]], l.res.app)
v.lt5    <- mapply(function(r) r[["lt5"]], l.res.app)

cat(sprintf("PHI waist %-18s  %4.1f (%.1f)  %4.1f (%.1f)  %5.1f%%       %5.1f%%\n",
            "OLS (n=1927)",
            mean(v.maxpct), sd(v.maxpct),
            mean(v.se), sd(v.se),
            100 * mean(v.sig),
            100 * mean(v.lt5)))

v.fr <- mapply(function(r) r[["fullrank"]], l.res.app)
cat(sprintf(paste0("PHI fits: %d/%d replications returned a fit; ",
                   "%d/%d bootstrap designs were full rank\n",
                   "(expected full-rank rate 1 - (1-1/n)^n = 63%%; ",
                   "rank-deficient draws lose the rarest housing level).\n"),
            length(l.res.app), s.B, sum(v.fr), length(v.fr)))
# +------------------------------------------------------------+
