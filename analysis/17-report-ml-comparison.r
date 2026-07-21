#!/usr/bin/env r
# +------------------------------------------------------------------------+
# | 17-report-ml-comparison.r -- Six-method comparison                     |
# +------------------------------------------------------------------------+
#
# The central empirical test of Proposition 2 (no-guarantee): compares
# the copula and projection generators against four ML generators.
# Section 6 ("Methods compared") explains: "The machine learning
# generators are selected to cover the major categories of learning
# objectives, so that the comparison tests whether the failure is
# structural (as Proposition 2 predicts) rather than architecture-
# specific."
#
# Methods:
#   1. copula     -- CF quantile generator with EB shrinkage (ours)
#   2. projection -- sufficient statistics projection (ours, Theorem 7)
#   3. ctgan      -- CTGAN [xu-modeling-2019]: adversarial, Wasserstein
#   4. tvae       -- TVAE: variational, ELBO (the VAE ELBO discussion)
#   5. gauscop    -- SDV GaussianCopula [patki-syntheti-2016]: pairwise corr
#   6. copgan     -- CopulaGAN: copula marginals + GAN dependence
#
# Key results (single seeded run):
#   - Projection: 0.0% max deviation (Sc1,3,4,8,9), 2.0% (Sc6, Poisson rounding)
#   - Copula: 20-48% (rare-level attenuation, Theorem 5)
#   - CTGAN: 98-210% (mode averaging, Section 4(i))
#   - TVAE: 29.8-935% (high variance, Section 4 VAE objective)
#   - GaussianCopula: 86-234% (wrong-sign coefficients, Theorem 4)
#   - CopulaGAN: 94-338%
#
# The comparison measures inferential fidelity (max|%|, sig%, dir%),
# distributional fidelity (avg KS), and privacy (DCR) -- all three
# to show that distributional closeness does not predict inferential
# closeness (the paper's central thesis).
#
# References:
#   Proposition (no inferential guarantee from distributional objectives)
#   Section (failure of distributional objectives for explanatory
#     modelling) -- why ML generators fail
#   GAN: Goodfellow et al. (2020) [goodfell-generati-2020]
#   VAE: Kingma & Welling (2014) [kingma-auto-enc-2014]
#   SDV: Patki et al. (2016) [patki-syntheti-2016]
#   CTGAN: Xu et al. (2019) [xu-modeling-2019]
#
# +------------------------------------------------------------------------+
# | Setup                                                                  |
# +------------------------------------------------------------------------+
.sobol.main <- TRUE
source("lib/00-shared.r")
source("lib/gen-syn-family.r")
source("lib/07-calibrated-generator.r")
source("lib/09-symmetric-generator.r")
source("lib/11-projection-generator.r")
rm(.sobol.main)
# +------------------------------------------------------------------------+
# | Helpers                                                                |
# +------------------------------------------------------------------------+
# check.row: computes inferential fidelity metrics (same as 14-report).
# fit.model: fits appropriate GLM per scenario (quasi-gamma for ZI).
# load.ml: reads CSV output from Python ML generators (13-gen-ml-default.py).
# dcr.sample: subsampled Gower DCR for privacy comparison.
pct <- function(a, b)
    `if`(abs(b) < 1e-10, `if`(abs(a) < 1e-10, 0, Inf),
         round(100 * (a - b) / abs(b), 1))

check.row <- function(mo.orig, mo.syn) {
    v.o <- coef(mo.orig)
    v.s <- coef(mo.syn)
    s.o <- summary(mo.orig)[["coefficients"]]
    s.s <- summary(mo.syn)[["coefficients"]]
    # align lengths (ML may drop levels)
    s.p <- min(nrow(s.o), nrow(s.s))
    v.o <- v.o[seq_len(s.p)]
    v.s <- v.s[seq_len(s.p)]
    v.pct <- abs(100 * (v.s - v.o) / abs(v.o))
    v.se.pct <- abs(100 * (s.s[seq_len(s.p), 2] -
                            s.o[seq_len(s.p), 2]) /
                    s.o[seq_len(s.p), 2])
    list(
        sig    = round(100 * mean((s.o[seq_len(s.p), 4] < 0.05) ==
                                  (s.s[seq_len(s.p), 4] < 0.05))),
        dir    = round(100 * mean(sign(v.o) == sign(v.s))),
        lt5    = round(100 * mean(v.pct < 5, na.rm = TRUE)),
        maxpct = round(max(v.pct, na.rm = TRUE), 1),
        se     = round(max(v.se.pct, na.rm = TRUE), 1))
}

fit.model <- function(d, family, link) {
    v.pred <- intersect(c("age", "female", "race", "income"), colnames(d))
    fml <- reformulate(v.pred, response = "response")
    tryCatch(
        `if`(family == "gaussian", lm(fml, d),
        `if`(family == "poisson", glm(fml, d, family = poisson()), {
            d2 <- `if`(link == "log" && any(d[["response"]] == 0),
                       d[d[["response"]] > 0, ], d)
            s.np <- ncol(model.matrix(fml, d2))
            glm(fml, d2, family = quasi(link = "log", variance = "mu^2"),
                start = c(1, rep(0, s.np - 1)),
                control = glm.control(maxit = 200))
        })),
        error = function(e) NULL)
}

load.ml <- function(i, method) {
    s.path <- sprintf("results/sim_scenario_%02d_%s.csv", i, method)
    if(!file.exists(s.path)) return(NULL)
    d <- read.csv(s.path)
    # ensure factors match
    if("female" %in% names(d)) d[["female"]] <- factor(d[["female"]])
    if("race" %in% names(d)) d[["race"]] <- factor(d[["race"]])
    d
}

dcr.sample <- function(d.syn, d.orig, n.sample = 200) {
    s.ns <- min(n.sample, nrow(d.syn))
    v.idx <- sample(nrow(d.syn), s.ns)
    v.types <- mapply(id.cat, d.orig)
    v.cnt <- names(v.types[v.types == "cnt"])
    v.cat <- names(v.types[v.types != "cnt"])
    l.range <- Map(function(nm) {
        s.r <- diff(range(d.orig[[nm]], na.rm = TRUE))
        `if`(s.r > 0, s.r, 1)
    }, v.cnt)
    s.p <- length(v.cnt) + length(v.cat)

    mapply(function(i) {
        s.min <- Inf
        Map(function(j) {
            s.d <- 0
            Map(function(nm) {
                s.d <<- s.d + abs(as.numeric(d.syn[[nm]][i]) -
                                  as.numeric(d.orig[[nm]][j])) /
                        l.range[[nm]]
            }, v.cnt)
            Map(function(nm) {
                s.d <<- s.d + as.integer(
                    as.character(d.syn[[nm]][i]) !=
                    as.character(d.orig[[nm]][j]))
            }, v.cat)
            s.min <<- min(s.min, s.d / s.p)
        }, seq_len(nrow(d.orig)))
        s.min
    }, v.idx)
}
# +------------------------------------------------------------------------+
# | Run comparison                                                         |
# +------------------------------------------------------------------------+
# Tests Scenarios 1, 3, 4, 6, 8, 9 -- the six the paper's ML tables report
# (race3, age, S2 six-method).  All six come from one seeded run, and the
# families and supports match 25-mc-allmethod.r, the MC driver.
# For each scenario and method, reports:
#   - Inferential: sig%, dir%, <5%, max|%|, SE% (Definition 2)
#   - Distributional: avg KS across continuous variables
#   - Privacy: median DCR (Gower distance to closest record)
load("results/d.sim.scenarios.RData")

# indexed by scenario number.  NULL where the scenario is not tested here.
# Sc.7 is the clean Gamma/log design (correlated predictors).  The zero-inflated
# design is Sc.8.
l.sup <- list(
    c(response = "real",     age = "real", income = "real"),
    NULL,
    c(response = "positive", age = "real", income = "positive"),
    c(response = "zeroinfl", age = "real", income = "positive"),
    NULL,
    c(response = "count",    age = "real", income = "positive"),
    c(response = "positive", age = "real", income = "positive"),
    c(response = "zeroinfl", age = "real", income = "positive"),
    c(response = "bounded",  age = "real", income = "positive"),
    NULL)

v.test <- c(1, 3, 4, 6, 8, 9)
l.fam  <- c("gaussian", "gamma", "quasi-gamma", "poisson", "quasi-gamma", "gaussian")
l.lnk  <- c("identity", "log", "log", "log", "log", "identity")
v.ml   <- c("ctgan", "tvae", "gauscop", "copgan")

cat("ML COMPARISON REPORT\n")
cat(strrep("=", 90), "\n\n")

Map(function(i, fam, lnk) {
    s   <- l.scenarios[[i]]
    d   <- s[["data"]]
    sup <- l.sup[[i]]
    mo.o <- fit.model(d, s[["family"]], s[["link"]])
    if(is.null(mo.o)) return(NULL)

    cat(sprintf("=== Scenario %d: %s(%s) ===\n", i, s[["family"]], s[["link"]]))
    cat(sprintf("%-14s  %5s %5s %5s %7s %7s   %6s  %7s\n",
                "method", "sig%", "dir%", "<5%", "max|%|", "SE%",
                "avgKS", "DCR.med"))
    cat(strrep("-", 90), "\n")

    # copula
    d.c <- gen.syn.family(d, order = 4, supports = sup)
    mo.c <- fit.model(d.c, s[["family"]], s[["link"]])
    if(!is.null(mo.c)) {
        r <- check.row(mo.o, mo.c)
        v.ks <- mapply(function(nm) {
            if(is.factor(d[[nm]])) return(NA)
            suppressWarnings(ks.test(d[[nm]], d.c[[nm]])[["statistic"]])
        }, names(d))
        s.ks <- round(mean(v.ks, na.rm = TRUE), 3)
        v.dcr <- dcr.sample(d.c, d, 100)
        cat(sprintf("%-14s  %4d%% %4d%% %4d%% %6.1f%% %6.1f%%   %.3f  %.4f\n",
                    "copula", r[["sig"]], r[["dir"]], r[["lt5"]],
                    r[["maxpct"]], r[["se"]], s.ks, median(v.dcr)))
    }

    # projection
    d.p <- gen.syn.projected(d, y.nm = "response", family = fam,
                             link = lnk, supports = sup)
    mo.p <- fit.model(d.p, s[["family"]], s[["link"]])
    if(!is.null(mo.p)) {
        r <- check.row(mo.o, mo.p)
        v.ks <- mapply(function(nm) {
            if(is.factor(d[[nm]])) return(NA)
            suppressWarnings(ks.test(d[[nm]], d.p[[nm]])[["statistic"]])
        }, names(d))
        s.ks <- round(mean(v.ks, na.rm = TRUE), 3)
        v.dcr <- dcr.sample(d.p, d, 100)
        cat(sprintf("%-14s  %4d%% %4d%% %4d%% %6.1f%% %6.1f%%   %.3f  %.4f\n",
                    "projection", r[["sig"]], r[["dir"]], r[["lt5"]],
                    r[["maxpct"]], r[["se"]], s.ks, median(v.dcr)))
    }

    # ML methods
    Map(function(ml.nm) {
        d.ml <- load.ml(i, ml.nm)
        if(is.null(d.ml)) {
            cat(sprintf("%-14s  (not available)\n", ml.nm))
            return(NULL)
        }
        mo.ml <- fit.model(d.ml, s[["family"]], s[["link"]])
        if(is.null(mo.ml)) {
            cat(sprintf("%-14s  (fit failed)\n", ml.nm))
            return(NULL)
        }
        r <- check.row(mo.o, mo.ml)
        v.ks <- mapply(function(nm) {
            if(is.factor(d[[nm]]) || is.factor(d.ml[[nm]])) return(NA)
            suppressWarnings(
                ks.test(as.numeric(d[[nm]]),
                        as.numeric(d.ml[[nm]]))[["statistic"]])
        }, intersect(names(d), names(d.ml)))
        s.ks <- round(mean(v.ks, na.rm = TRUE), 3)
        v.dcr <- dcr.sample(d.ml, d, 100)
        cat(sprintf("%-14s  %4d%% %4d%% %4d%% %6.1f%% %6.1f%%   %.3f  %.4f\n",
                    ml.nm, r[["sig"]], r[["dir"]], r[["lt5"]],
                    r[["maxpct"]], r[["se"]], s.ks, median(v.dcr)))
    }, v.ml)

    cat("\n")
}, v.test, l.fam, l.lnk)
# +------------------------------------------------------------+
