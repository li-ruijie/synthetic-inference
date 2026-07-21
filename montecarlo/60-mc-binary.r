#!/usr/bin/env r
# +------------------------------------------------------------------------+
# | 60-mc-binary.r -- Monte Carlo validation of the binomial projection    |
# +------------------------------------------------------------------------+
#
# Scenario 11 (binomial/logit, defined in 00-simulate.r) is the first
# scenario to exercise Theorem 7's binomial branch (path B in
# 11-projection-generator.r, the Stage-2 candidate P2). Path B
# matches the logistic score X'(y* - mu_hat) = 0 and shrinks the score-orthogonal
# copula residual to keep y* in [0, 1], giving an exact, valid fractional
# binomial response in one closed-form step.
#
# This script confirms the projection meets the same 0.0% deviation bar the
# continuous families meet, and contrasts it with the copula baseline, which
# distorts the logistic coefficients like every other scenario. B = 100
# replications, master seed 6489, deterministic under OPENBLAS_NUM_THREADS=1.
#
# Output: mc-binary.csv and a console table.
#
# +------------------------------------------------------------------------+
# | Setup (mirror 18-mc-projection.r: dotted .sobol.main survives the       |
# | globalenv clear in 00-shared.r, suppressing the generator harness)      |
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
# | Deviation metrics (aligned by coefficient name, as 18-mc-projection.r)  |
# +------------------------------------------------------------------------+
dev.metrics <- function(mo.o, mo.s) {
    if(is.null(mo.s)) return(NULL)
    v.o   <- coef(mo.o)
    v.s   <- coef(mo.s)[names(v.o)]
    s.o   <- summary(mo.o)[["coefficients"]]
    s.s   <- summary(mo.s)[["coefficients"]]
    v.se.o <- s.o[, 2][names(v.o)]
    v.se.s <- s.s[, 2][names(v.o)]
    v.p4.o <- s.o[, 4][names(v.o)]
    v.p4.s <- s.s[, 4][names(v.o)]
    v.pct  <- abs(100 * (v.s - v.o) / abs(v.o))
    list(maxpct = max(v.pct, na.rm = TRUE),
         se     = max(abs(100 * (v.se.s - v.se.o) / abs(v.se.o)), na.rm = TRUE),
         sig    = mean((v.p4.o < 0.05) == (v.p4.s < 0.05), na.rm = TRUE),
         dir    = mean(sign(v.o) == sign(v.s), na.rm = TRUE),
         lt5    = mean(v.pct < 5, na.rm = TRUE))
}

agg.method <- function(l.res, s.method) {
    l.res <- Filter(Negate(is.null), l.res)
    v.maxpct <- vapply(l.res, `[[`, numeric(1), "maxpct")
    v.se     <- vapply(l.res, `[[`, numeric(1), "se")
    v.sig    <- vapply(l.res, `[[`, numeric(1), "sig")
    v.dir    <- vapply(l.res, `[[`, numeric(1), "dir")
    v.lt5    <- vapply(l.res, `[[`, numeric(1), "lt5")
    data.frame(method = s.method, n = length(l.res),
               maxpct_mean = mean(v.maxpct), maxpct_sd = sd(v.maxpct),
               se_mean = mean(v.se), se_sd = sd(v.se),
               sig = 100 * mean(v.sig), dir = 100 * mean(v.dir),
               lt5 = 100 * mean(v.lt5), stringsAsFactors = FALSE)
}

# +------------------------------------------------------------------------+
# | Run Monte Carlo on Scenario 11                                          |
# +------------------------------------------------------------------------+
load("results/d.sim.scenarios.RData")
s.B <- 100
s   <- l.scenarios[[11]]
d   <- s[["data"]]
stopifnot(s[["family"]] == "binomial", s[["link"]] == "logit")

# support of the two continuous predictors (income is the skewed log-normal for
# scenario >= 2). The binary response is categorical (bootstrapped by the
# copula), so it needs no support token.
sup <- c(age = "real", income = "positive")
fml <- response ~ age + female + race + income
mo.o <- glm(fml, d, family = binomial())

v.seeds <- 6489 + seq_len(s.B)

# projection: exact logistic refit expected (path B)
l.proj <- Map(function(seed) {
    d.p <- tryCatch(
        gen.syn.projected(d, y.nm = "response", family = "binomial",
                          link = "logit", supports = sup, seed = seed),
        error = function(e) NULL)
    if(is.null(d.p)) return(NULL)
    # glm(binomial) warns "non-integer #successes" on the fractional [0, 1]
    # projected release, which is expected and harmless (weights = 1).  Suppress
    # it as the package test-projection.R case does.
    dev.metrics(mo.o, tryCatch(suppressWarnings(glm(fml, d.p, family = binomial())),
                               error = function(e) NULL))
}, v.seeds)

# copula baseline: distributional match only, distorts logistic coefficients
l.cop <- Map(function(seed) {
    d.c <- tryCatch(
        gen.syn.family(d, order = 4, supports = sup, shrink = TRUE, seed = seed),
        error = function(e) NULL)
    if(is.null(d.c)) return(NULL)
    dev.metrics(mo.o, tryCatch(suppressWarnings(glm(fml, d.c, family = binomial())),
                               error = function(e) NULL))
}, v.seeds)

d.out <- rbind(agg.method(l.proj, "projection"),
               agg.method(l.cop,  "copula"))
write.csv(d.out, "results/mc-binary.csv", row.names = FALSE)

cat(sprintf("MONTE CARLO BINARY SCENARIO 11 binomial(logit) (%d replications)\n", s.B))
cat(strrep("=", 78), "\n\n")
cat(sprintf("%-12s %6s %14s %14s %7s %7s %7s\n",
            "method", "n", "max|%| (sd)", "SE% (sd)", "sig%", "dir%", "<5%"))
cat(strrep("-", 78), "\n")
Map(function(k) {
    cat(sprintf("%-12s %6d %6.2f (%5.2f) %6.2f (%5.2f) %6.1f %6.1f %6.1f\n",
                d.out[["method"]][k], d.out[["n"]][k],
                d.out[["maxpct_mean"]][k], d.out[["maxpct_sd"]][k],
                d.out[["se_mean"]][k], d.out[["se_sd"]][k],
                d.out[["sig"]][k], d.out[["dir"]][k], d.out[["lt5"]][k]))
}, seq_len(nrow(d.out)))
cat("\nWrote mc-binary.csv\n")
# +------------------------------------------------------------------------+
# | End of 60-mc-binary.r                                                   |
# +------------------------------------------------------------------------+
