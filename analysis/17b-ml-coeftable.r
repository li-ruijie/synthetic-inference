# +------------------------------------------------------------------------+
# | 17b-ml-coeftable.r -- per-coefficient ML comparison, emitted as a CSV   |
# +------------------------------------------------------------------------+
#
# Why this script exists.
#
# The race, age and six-method coefficient tables of the supplement all report
# Per-coefficient values for the ML generators.
# 17-report-ml-comparison.r prints only the summary metrics (sig%, dir%, <5%,
# max|%|, SE%), so the per-coefficient values need a driver of their own.
#
# This script is that driver.  One seeded run over the six scenarios the paper
# tabulates, every method, every coefficient, written to a CSV the tables can
# be diffed against.  Per-coefficient ML values should be cited from here and
# from nowhere else.
#
# Generator calls mirror 25-mc-allmethod.r, the MC driver, so the single run and
# the 100-replication run describe the same estimators.
#
# Input:   d.sim.scenarios.RData, sim_scenario_{01,03,04,06,08,09}_{ml}.csv
#          (the ML CSVs come from 13-gen-ml-default.py, which is seeded)
# Output:  ml-coeftable.csv
#
# Usage:   Rscript 17b-ml-coeftable.r
# +------------------------------------------------------------------------+
# Both flags are needed, as 25-mc-allmethod.r sets them.  .shared.loaded stops
# 00-shared.r re-initialising, and its re-init clears the global environment
# (ut.clear.globalwspc).  07- and 11- each source 00-shared.r themselves, so
# without the flag, sourcing 11- wipes the gen.syn.calibrated that 07- had just
# defined, and the calibrated arm silently vanishes from the output.
# .sobol.main suppresses the demo blocks in 07- and 11-.
.shared.loaded <- TRUE
.sobol.main <- TRUE
source("lib/00-shared.r")

s.wd <- tryCatch(
    normalizePath(dirname(sys.frame(1)[["ofile"]]), mustWork = TRUE),
    error = function(e) getwd())
setwd(s.wd)

source("lib/gen-syn-family.r")
source("lib/07-calibrated-generator.r")
source("lib/11-projection-generator.r")
rm(.sobol.main)

stopifnot(is.function(gen.syn.family), is.function(gen.syn.calibrated),
          is.function(gen.syn.projected))

SEED <- 6489L

v.scenarios <- c(1L, 3L, 4L, 6L, 8L, 9L)
l.fam <- list(`1` = "gaussian", `3` = "gamma",      `4` = "quasi-gamma",
              `6` = "poisson",  `8` = "quasi-gamma", `9` = "gaussian")
l.lnk <- list(`1` = "identity", `3` = "log", `4` = "log",
              `6` = "log",      `8` = "log", `9` = "identity")
# Scenario 9 is the bounded [0, 100] response (generalised logit), matching
# 07-calibrated-generator.r and the DGP.  A "real" support here would make the
# Monte Carlo skip the transform the paper tests.
l.sup <- list(
    `1` = c(response = "real",     age = "real", income = "real"),
    `3` = c(response = "positive", age = "real", income = "positive"),
    `4` = c(response = "zeroinfl", age = "real", income = "positive"),
    `6` = c(response = "count",    age = "real", income = "positive"),
    `8` = c(response = "zeroinfl", age = "real", income = "positive"),
    `9` = c(response = "bounded",  age = "real", income = "positive"))

v.ml <- c("ctgan", "tvae", "gauscop", "copgan")

fit.model <- function(d, family, link) {
    v.pred <- intersect(c("age", "female", "race", "income"), colnames(d))
    fml <- reformulate(v.pred, response = "response")
    # The zero-drop applies to the two-part quasi families only.  A zero
    # count is a legitimate Poisson outcome, and the projected count draw
    # can carry probabilistic-rounding zeros, so filtering the Poisson
    # refit would fit a truncated sample against an untruncated original.
    d.fit <- `if`(link == "log" && any(d[["response"]] == 0),
                  d[d[["response"]] > 0, ], d)
    if(nrow(d.fit) < 20L) return(NULL)
    tryCatch(
        `if`(family == "gaussian",
             lm(fml, d.fit),
        `if`(family == "binomial",
             suppressWarnings(glm(fml, d.fit, family = binomial())),
        `if`(family == "poisson",
             glm(fml, d, family = poisson()),
             glm(fml, d.fit,
                 family = quasi(link = "log", variance = "mu^2"),
                 start = c(1, rep(0, length(v.pred) + 1L)),
                 control = glm.control(maxit = 100L))))),
        error = function(e) NULL)
}

load.ml <- function(sc, method) {
    s.path <- sprintf("results/sim_scenario_%02d_%s.csv", sc, method)
    if(!file.exists(s.path)) return(NULL)
    d <- read.csv(s.path)
    if("female" %in% names(d)) d[["female"]] <- factor(d[["female"]])
    if("race"   %in% names(d)) d[["race"]]   <- factor(d[["race"]])
    d
}

load("results/d.sim.scenarios.RData")

l.rows <- Map(function(sc) {
    k   <- as.character(sc)
    s   <- l.scenarios[[sc]]
    d   <- s[["data"]]
    sup <- l.sup[[k]]
    fam <- l.fam[[k]]
    lnk <- l.lnk[[k]]

    mo.o <- fit.model(d, s[["family"]], s[["link"]])
    if(is.null(mo.o)) stop(sprintf("scenario %d: original fit failed", sc))

    # R generators, same calls as 25-mc-allmethod.r.  Errors are reported, not
    # swallowed: a tryCatch that returns NULL silently is how this project has
    # shipped an empty CSV with exit code 0 before.
    warn.null <- function(nm, e) {
        cat(sprintf("  !! %s generator FAILED on scenario %d: %s\n",
                    nm, sc, conditionMessage(e)))
        NULL
    }
    set.seed(SEED)
    d.cop <- tryCatch(gen.syn.family(d, order = 4, supports = sup, seed = SEED),
                      error = function(e) warn.null("copula", e))
    d.cal <- tryCatch(gen.syn.calibrated(d, order.cf = 4, order.cal = 2,
                                         supports = sup, seed = SEED),
                      error = function(e) warn.null("calibrated", e))
    d.prj <- tryCatch(gen.syn.projected(d, y.nm = "response", family = fam,
                                        link = lnk, supports = sup, seed = SEED),
                      error = function(e) warn.null("projection", e))

    l.d <- c(list(copula = d.cop, calibrated = d.cal, projection = d.prj),
             setNames(Map(function(m) load.ml(sc, m), v.ml), v.ml))

    v.o <- coef(mo.o)
    l.out <- Map(function(nm) {
        d.s <- l.d[[nm]]
        if(is.null(d.s)) {
            cat(sprintf("  !! %s: no data on scenario %d\n", nm, sc))
            return(NULL)
        }
        mo.s <- fit.model(d.s, s[["family"]], s[["link"]])
        if(is.null(mo.s)) {
            cat(sprintf("  !! %s: FIT FAILED on scenario %d\n", nm, sc))
            return(NULL)
        }
        v.s <- coef(mo.s)
        v.common <- intersect(names(v.o), names(v.s))
        data.frame(scenario = sc, method = nm, coefficient = v.common,
                   original = as.numeric(v.o[v.common]),
                   estimate = as.numeric(v.s[v.common]),
                   stringsAsFactors = FALSE)
    }, names(l.d))

    d.sc <- do.call(rbind, Filter(Negate(is.null), l.out))
    cat(sprintf("scenario %2d (%s/%s): %d rows\n", sc, fam, lnk, nrow(d.sc)))
    d.sc
}, v.scenarios)

d.all <- do.call(rbind, l.rows)
d.all[["pct.dev"]] <- ifelse(
    abs(d.all[["original"]]) < 1e-10, NA_real_,
    100 * (d.all[["estimate"]] - d.all[["original"]]) / abs(d.all[["original"]]))

write.csv(d.all, "results/ml-coeftable.csv", row.names = FALSE)
cat(sprintf("\nWritten to results/ml-coeftable.csv (%d rows, %d scenarios, %d methods)\n",
            nrow(d.all), length(unique(d.all[["scenario"]])),
            length(unique(d.all[["method"]]))))
