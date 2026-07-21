#!/usr/bin/env r
# +------------------------------------------------------------------------+
# | 14-report-inferential.r -- Inferential fidelity report                 |
# +------------------------------------------------------------------------+
#
# Generates the primary evaluation table for the paper (Section 6,
# "Results"). For each scenario, compares original vs copula vs
# projection on five criteria that operationalise Definition 2
# (inferential fidelity):
#
#   sig%   -- Significance preservation: the fraction of coefficients
#             retaining the same significance status (p < 0.05).
#             (Proposition 1(d): significance preserved when moment-
#             truncation error < margin from critical value)
#   dir%   -- Direction preservation: whether all coefficient signs match.
#   <5%    -- Fraction of coefficients with |% deviation| < 5%
#   max|%| -- Maximum absolute percentage deviation across all
#             coefficients. This is the primary metric in Section 6.
#             Theorem 7 predicts 0.0% for projection on canonical GLMs.
#   SE%    -- Maximum SE deviation. Proposition 1(b)-(c) bound the
#             dispersion and SE error by O(nu_(k+1)/(k+1)!) when k >= 2
#             on the h-scale.  Neither clause makes it zero.
#
# The copula column validates Theorem 2 (partial moment matching):
# bias should decay with moment order. The projection column validates
# Theorem 7 (sufficient statistics projection): exact preservation.
#
# References:
#   Proposition (sufficient conditions for summary statistics)
#     -- k=2 controls the full summary
#   Theorem (sufficient statistics projection)
#     -- exact fidelity across all submodels
#   Yap et al. (2020) [yap-developm-2020] -- PHI dataset (real-data test)
#
# +------------------------------------------------------------------------+
# | Setup                                                                  |
# +------------------------------------------------------------------------+
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
# +------------------------------------------------------------------------+
# | Helpers                                                                |
# +------------------------------------------------------------------------+
# pct: percentage deviation of synthetic coefficient from original.
# check.row: computes the five fidelity criteria for one scenario.
# fit.model: fits the appropriate GLM family for each scenario.
# For quasi-gamma (Scenarios 4-5, 7-8), uses Wedderburn's quasi-likelihood
# [wedderbu-quasi-li-1974] with variance function mu^2 (gamma-like).
# For zero-inflated data, fits on positive subset only (two-part model).
pct <- function(a, b)
    round(100 * (a - b) / abs(b), 1)

check.row <- function(mo.orig, mo.syn) {
    v.o  <- coef(mo.orig)
    v.s  <- coef(mo.syn)
    s.o  <- summary(mo.orig)[["coefficients"]]
    s.s  <- summary(mo.syn)[["coefficients"]]
    v.pct    <- abs(100 * (v.s - v.o) / abs(v.o))
    v.se.pct <- abs(100 * (s.s[, 2] - s.o[, 2]) / s.o[, 2])
    list(
        sig    = round(100 * mean((s.o[, 4] < 0.05) == (s.s[, 4] < 0.05))),
        dir    = round(100 * mean(sign(v.o) == sign(v.s))),
        lt5    = round(100 * mean(v.pct < 5, na.rm = TRUE)),
        maxpct = round(max(v.pct, na.rm = TRUE), 1),
        se     = round(max(v.se.pct, na.rm = TRUE), 1))
}

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
# +------------------------------------------------------------------------+
# | Simulation scenarios                                                   |
# +------------------------------------------------------------------------+
# Test on the scenarios the paper reports:
# Sc1 (Gaussian/identity), Sc3 (Gamma/log), Sc4 (ZI quasi-gamma),
# Sc6 (Poisson/log), Sc7 (clean Gamma, correlated predictors),
# Sc8 (ZI all combined), Sc9 (bounded [0,100]).
# Each scenario is generated with copula (gen.syn.family) and
# projection (gen.syn.projected), then fitted with the same GLM.
# The l.sup list specifies the support type for each variable,
# which determines the support transform h (Definition 3).
load("results/d.sim.scenarios.RData")

l.sup <- list(
    c(response = "real",     age = "real", income = "real"),        # 1
    NULL,                                                            # 2
    c(response = "positive", age = "real", income = "positive"),   # 3
    c(response = "zeroinfl", age = "real", income = "positive"),   # 4
    NULL,                                                            # 5
    c(response = "count",    age = "real", income = "positive"),   # 6
    c(response = "positive", age = "real", income = "positive"),   # 7
    c(response = "zeroinfl", age = "real", income = "positive"),   # 8
    c(response = "bounded",  age = "real", income = "positive"),   # 9
    NULL)

v.test <- c(1, 3, 4, 6, 7, 8, 9)
l.fam  <- c("gaussian", "gamma", "quasi-gamma",
             "poisson", "gamma", "quasi-gamma", "gaussian")
l.lnk  <- c("identity", "log", "log", "log", "log", "log", "identity")

cat("INFERENTIAL FIDELITY REPORT\n")
cat(strrep("=", 80), "\n\n")

cat(sprintf("%-25s  %5s %5s %5s %7s %7s\n",
            "", "sig%", "dir%", "<5%", "max|%|", "SE%"))
cat(strrep("-", 80), "\n")

Map(function(i, fam, lnk) {
    s   <- l.scenarios[[i]]
    d   <- s[["data"]]
    sup <- l.sup[[i]]
    mo.o <- fit.model(d, s[["family"]], s[["link"]])
    if(is.null(mo.o)) return(NULL)

    d.c <- gen.syn.family(d, order = 4, supports = sup)
    d.p <- gen.syn.projected(d, y.nm = "response", family = fam,
                             link = lnk, supports = sup)
    mo.c <- fit.model(d.c, s[["family"]], s[["link"]])
    mo.p <- fit.model(d.p, s[["family"]], s[["link"]])
    if(is.null(mo.c) || is.null(mo.p)) return(NULL)

    r.c <- check.row(mo.o, mo.c)
    r.p <- check.row(mo.o, mo.p)

    s.label <- sprintf("Sc%d %s(%s)", i, s[["family"]], s[["link"]])
    cat(sprintf("%-25s  %4d%% %4d%% %4d%% %6.1f%% %6.1f%%  copula\n",
                s.label, r.c[["sig"]], r.c[["dir"]], r.c[["lt5"]],
                r.c[["maxpct"]], r.c[["se"]]))
    cat(sprintf("%-25s  %4d%% %4d%% %4d%% %6.1f%% %6.1f%%  projection\n",
                "", r.p[["sig"]], r.p[["dir"]], r.p[["lt5"]],
                r.p[["maxpct"]], r.p[["se"]]))
}, v.test, l.fam, l.lnk)

# +------------------------------------------------------------------------+
# | Real-data application: PHI waist circumference                        |
# +------------------------------------------------------------------------+
# Section 6 ("Real-data application"): "As a test on real data, we apply
# the generators to the Population Health Index (PHI) study
# [yap-developm-2020]: n = 1927 complete cases of the 1942 survey
# participants, p = 16, and waist circumference modelled on age, sex,
# race, education, and housing type".
#
# The cross-classification female(2) x race(4) x edu(3) x housing(9) gives
# 216 possible cross-strata, of which 123 are occupied: 26 singletons and 15
# doubletons.  This is the sparse strata extreme case (Theorem 5):
# the copula fails (1084% deviation on housing) while the projection achieves
# 0.0% (Theorem 7 guarantee).
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
d.c <- gen.syn.family(d.app, order = 4, supports = sup.app)
d.p <- gen.syn.projected(d.app, y.nm = "waist", family = "gaussian",
                         link = "identity", supports = sup.app)

fml.app <- waist ~ age + female + race + edu + housing
mo.o <- lm(fml.app, d.app)
mo.c <- lm(fml.app, d.c)
mo.p <- lm(fml.app, d.p)

r.c <- check.row(mo.o, mo.c)
r.p <- check.row(mo.o, mo.p)

cat(sprintf("%-25s  %4d%% %4d%% %4d%% %6.1f%% %6.1f%%  copula\n",
            "PHI waist OLS (n=1927)", r.c[["sig"]], r.c[["dir"]],
            r.c[["lt5"]], r.c[["maxpct"]], r.c[["se"]]))
cat(sprintf("%-25s  %4d%% %4d%% %4d%% %6.1f%% %6.1f%%  projection\n",
            "", r.p[["sig"]], r.p[["dir"]], r.p[["lt5"]],
            r.p[["maxpct"]], r.p[["se"]]))
# +------------------------------------------------------------------------+
# | End of 14-report-inferential.r                                         |
# +------------------------------------------------------------------------+
