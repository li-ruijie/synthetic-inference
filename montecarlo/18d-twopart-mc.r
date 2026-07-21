# +------------------------------------------------------------------------+
# | 18d-twopart-mc.r -- 100-replication Monte Carlo of the two-part        |
# |                     projection's standard-error exactness              |
# +------------------------------------------------------------------------+
#
# Why this script exists.
#
# 18c-twopart-se.r shows on a single seed that the two-part (zero-pattern
# preserving) projection returns the positive-block standard errors to
# exactness on the two zero-inflated scenarios (Sc.4, Sc.8).  This script
# confirms the result is not a one-seed artefact by repeating it over 100
# independent copula draws (seeds 6489+1..100 from the master seed 6489 of
# Section 6) and reporting the positive-block coefficient and SE deviation as a
# Monte Carlo mean and standard deviation, matching the basis of the S11
# all-method tables.
#
# The original data are fixed, as in 25-mc-allmethod.r, and each replication
# varies the copula seed, so the spread measures the generator's stochasticity
# under the two-part release.  Proposition S5 (supplement S13) makes the
# positive-block standard errors exact whenever the release preserves the zero
# pattern, so the mean deviation sits at the projection's numerical tolerance
# with a near-zero standard deviation, and the zero block is exact by
# construction on every replication.
#
# Output:  twopart-se-mc.csv
# Usage:   Rscript 18d-twopart-mc.r
#
# Both source-guard flags are required (see 18c-, 25-mc-allmethod.r).
.shared.loaded <- TRUE
.sobol.main <- TRUE
source("lib/00-shared.r")

s.wd <- tryCatch(
    normalizePath(dirname(sys.frame(1)[["ofile"]]), mustWork = TRUE),
    error = function(e) getwd())
setwd(s.wd)

source("lib/gen-syn-family.r")
source("lib/11-projection-generator.r")
rm(.sobol.main)

s.B <- 100L
s.master.seed <- 6489L

fit.pos <- function(d) {
    dp <- d[d[["response"]] > 0, ]
    glm(response ~ age + female + race + income, dp,
        family = quasi(link = "log", variance = "mu^2"),
        start = c(1, rep(0, 5)), control = glm.control(maxit = 100L))
}
fit.zero <- function(d) {
    d[["z"]] <- as.integer(d[["response"]] > 0)
    glm(z ~ age + female + race + income, d, family = binomial())
}
maxpct <- function(a, b) max(abs(100 * (a - b) / b), na.rm = TRUE)

sup <- c(response = "zeroinfl", age = "real", income = "positive")

load("results/d.sim.scenarios.RData")

# One replication: two-part release on the fixed original data at copula seed
# `seed`, refit both blocks, and return the coefficient and SE deviations
# against the original fits passed in.
mc.one <- function(d, seed, coef.pos.o, se.pos.o, coef.zero.o, se.zero.o) {
    d.fix <- gen.syn.projected.twopart(d, y.nm = "response",
                                       family = "quasi-gamma", link = "log",
                                       supports = sup, seed = seed)
    mo.pos.x  <- fit.pos(d.fix)
    mo.zero.x <- fit.zero(d.fix)
    list(
        pos_coef_dev  = maxpct(coef(mo.pos.x), coef.pos.o),
        pos_se_dev    = maxpct(summary(mo.pos.x)[["coefficients"]][, 2], se.pos.o),
        zero_coef_dev = maxpct(coef(mo.zero.x), coef.zero.o),
        zero_se_dev   = maxpct(summary(mo.zero.x)[["coefficients"]][, 2], se.zero.o))
}

l.rows <- Map(function(sc) {
    d <- l.scenarios[[sc]][["data"]]

    mo.pos.o  <- fit.pos(d)
    mo.zero.o <- fit.zero(d)
    coef.pos.o  <- coef(mo.pos.o)
    se.pos.o    <- summary(mo.pos.o)[["coefficients"]][, 2]
    coef.zero.o <- coef(mo.zero.o)
    se.zero.o   <- summary(mo.zero.o)[["coefficients"]][, 2]

    v.seeds <- s.master.seed + seq_len(s.B)
    l.res <- Map(function(seed) mc.one(d, seed, coef.pos.o, se.pos.o,
                                       coef.zero.o, se.zero.o),
                 v.seeds)

    v.pos.coef  <- mapply(function(r) r[["pos_coef_dev"]],  l.res)
    v.pos.se    <- mapply(function(r) r[["pos_se_dev"]],    l.res)
    v.zero.coef <- mapply(function(r) r[["zero_coef_dev"]], l.res)
    v.zero.se   <- mapply(function(r) r[["zero_se_dev"]],   l.res)

    m.pos.coef <- mean(v.pos.coef); d.pos.coef <- sd(v.pos.coef)
    m.pos.se   <- mean(v.pos.se);   d.pos.se   <- sd(v.pos.se)
    m.zero.coef <- mean(v.zero.coef)
    m.zero.se   <- mean(v.zero.se)

    cat(sprintf(paste0("sc%d (%d reps)  positive block: coef %.6f%% (sd %.6f)  ",
                       "SE %.6f%% (sd %.6f)   zero block: coef %.6f%% SE %.6f%%\n"),
                sc, length(l.res),
                m.pos.coef, d.pos.coef, m.pos.se, d.pos.se,
                m.zero.coef, m.zero.se))

    data.frame(
        scenario = sc, reps = length(l.res),
        pos_coef_dev_mean = m.pos.coef, pos_coef_dev_sd = d.pos.coef,
        pos_se_dev_mean   = m.pos.se,   pos_se_dev_sd   = d.pos.se,
        zero_coef_dev_mean = m.zero.coef, zero_se_dev_mean = m.zero.se)
}, c(4L, 8L))

d.out <- do.call(rbind, l.rows)
write.csv(d.out, "results/twopart-se-mc.csv", row.names = FALSE)
cat(sprintf("\nWritten to twopart-se-mc.csv (%d scenarios, %d replications)\n",
            nrow(d.out), s.B))
