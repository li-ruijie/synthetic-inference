#!/usr/bin/env Rscript
# +------------------------------------------------------------+
# |  Monte Carlo replications: All methods (parallel)          |
# +------------------------------------------------------------+
#
# Phase 2 of the MC pipeline.  Phase 1 (26-mc-ml-generate.py)
# generates ML synthetic datasets in parallel via Python/SDV.
# This script:
#   1. Generates copula/calibrated/projection datasets in parallel
#      via par_map (fork or PSOCK by platform).
#   2. Reads the pre-generated ML CSVs from results/mc/.
#   3. Computes diagnostics for all 7 methods.
#
# Output:
#   mc-allmethod-results.rds  -- full results list
#   mc-allmethod-summary.csv  -- mean (sd) per method per scenario
#
# Usage:
#   Rscript 25-mc-allmethod.r
#
# +------------------------------------------------------------+
# |  Setup                                                     |
# +------------------------------------------------------------+
options(width = 128)

.shared.loaded <- TRUE
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
source("lib/parallel_setup.R")
rm(.sobol.main)

# +------------------------------------------------------------+
# |  Configuration                                             |
# +------------------------------------------------------------+
s.B <- 100L
s.master.seed <- 6489L
v.scenarios <- c(1L, 3L, 4L, 6L, 8L, 9L)
v.ml.methods <- c("ctgan", "tvae", "gauscop", "copgan")

l.fam <- list(`1` = "gaussian", `3` = "gamma",     `4` = "quasi-gamma",
              `6` = "poisson",  `8` = "quasi-gamma", `9` = "gaussian")
l.lnk <- list(`1` = "identity", `3` = "log", `4` = "log",
              `6` = "log",      `8` = "log", `9` = "identity")
# Scenario 9 is the bounded [0, 100] response, and the paper says so
# ("Scenario 9 tests the bounded support transform").  A "real" support here
# would stop the Monte Carlo applying the generalised-logit transform the paper
# tests, and would treat a skewed income (skewness 1.96) as Gaussian.  This
# matches 07-calibrated-generator.r and the DGP.
l.sup <- list(
    `1` = c(response = "real",     age = "real", income = "real"),
    `3` = c(response = "positive", age = "real", income = "positive"),
    `4` = c(response = "zeroinfl", age = "real", income = "positive"),
    `6` = c(response = "count",    age = "real", income = "positive"),
    `8` = c(response = "zeroinfl", age = "real", income = "positive"),
    `9` = c(response = "bounded",  age = "real", income = "positive"))

# Seeds are discovered from results/mc/ filenames (written by
# 26-mc-ml-generate.py).  No cross-language seed coordination needed.

# +------------------------------------------------------------+
# |  Helpers                                                   |
# +------------------------------------------------------------+
fit.model <- function(d, family, link) {
    v.pred <- intersect(c("age", "female", "race", "income"),
                        colnames(d))
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

diagnostics <- function(mo.orig, mo.syn) {
    if (is.null(mo.syn)) return(NULL)
    v.o <- coef(mo.orig)
    v.s <- coef(mo.syn)
    s.p <- min(length(v.o), length(v.s))
    v.o <- v.o[seq_len(s.p)]
    v.s <- v.s[seq_len(s.p)]
    s.o <- summary(mo.orig)[["coefficients"]]
    s.s <- summary(mo.syn)[["coefficients"]]
    v.pct <- abs(100 * (v.s - v.o) / abs(v.o))
    v.se  <- abs(100 * (s.s[seq_len(s.p), 2] -
                         s.o[seq_len(s.p), 2]) / s.o[seq_len(s.p), 2])
    list(maxpct = max(v.pct, na.rm = TRUE),
         se     = max(v.se,  na.rm = TRUE),
         sig    = mean((s.o[seq_len(s.p), 4] < 0.05) ==
                       (s.s[seq_len(s.p), 4] < 0.05)),
         dir    = mean(sign(v.o) == sign(v.s)),
         lt5    = mean(v.pct < 5, na.rm = TRUE),
         coefs  = v.s)
}

# Read ML CSV from results/mc/.  Returns data.frame or NULL.
read.ml.csv <- function(sc, method, seed) {
    f <- sprintf("results/mc/sc%02d_seed%d_%s.csv", sc, seed, method)
    if (file.exists(f)) {
        tryCatch(read.csv(f), error = function(e) NULL)
    } else NULL
}

# Discover which seeds the Python script generated for a scenario.
discover.seeds <- function(sc) {
    v.files <- list.files("results/mc", pattern = sprintf("^sc%02d_seed", sc))
    v.seeds <- sub(sprintf("^sc%02d_seed(\\d+)_.*\\.csv$", sc),
                   "\\1", v.files)
    sort(unique(as.integer(v.seeds)))
}

# One replication: generate R methods + read ML CSVs + diagnostics.
mc.one.rep <- function(seed, d, sc, fam, lnk, sup, mo.orig) {
    set.seed(seed)

    # R generators
    d.cop <- tryCatch(
        gen.syn.family(d, order = 4, supports = sup, seed = seed),
        error = function(e) NULL)
    d.cal <- tryCatch(
        gen.syn.calibrated(d, order.cf = 4, order.cal = 2,
                           supports = sup, seed = seed),
        error = function(e) NULL)
    d.prj <- tryCatch(
        gen.syn.projected(d, y.nm = "response", family = fam,
                          link = lnk, supports = sup, seed = seed),
        error = function(e) NULL)

    # ML generators (read pre-generated CSVs)
    d.ctgan   <- read.ml.csv(sc, "ctgan",   seed)
    d.tvae    <- read.ml.csv(sc, "tvae",    seed)
    d.gauscop <- read.ml.csv(sc, "gauscop", seed)
    d.copgan  <- read.ml.csv(sc, "copgan",  seed)

    list(
        scenario   = sc,
        seed       = seed,
        copula     = diagnostics(mo.orig, fit.model(d.cop, fam, lnk)),
        calibrated = diagnostics(mo.orig, fit.model(d.cal, fam, lnk)),
        projection = diagnostics(mo.orig, fit.model(d.prj, fam, lnk)),
        ctgan      = diagnostics(mo.orig, fit.model(d.ctgan, fam, lnk)),
        tvae       = diagnostics(mo.orig, fit.model(d.tvae, fam, lnk)),
        gauscop    = diagnostics(mo.orig, fit.model(d.gauscop, fam, lnk)),
        copgan     = diagnostics(mo.orig, fit.model(d.copgan, fam, lnk)))
}

# +------------------------------------------------------------+
# |  Load simulation data                                      |
# +------------------------------------------------------------+
if (!file.exists("results/d.sim.scenarios.RData")) source("simulate/00-simulate.r")
load("results/d.sim.scenarios.RData")

# +------------------------------------------------------------+
# |  Parallel setup                                            |
# +------------------------------------------------------------+
cfg <- setup_parallel(
    n_cores = max(1L, parallel::detectCores() - 2L),
    export = c("fit.model", "diagnostics", "read.ml.csv"),
    init_worker = function() {
        Sys.setenv(OPENBLAS_NUM_THREADS = "1",
                   MKL_NUM_THREADS      = "1",
                   OMP_NUM_THREADS      = "1")
        # PSOCK workers start with an empty global environment, so they
        # source the generator stack the way the script head does.  The
        # flags skip the globalenv clear in 00-shared.r, and the guard
        # makes this a no-op on the fork and sequential backends, whose
        # parent process already loaded the stack.
        if (!exists("gen.syn.family", envir = globalenv())) {
            assign(".shared.loaded", TRUE, envir = globalenv())
            assign(".sobol.main", TRUE, envir = globalenv())
            source("lib/00-shared.r")
            source("lib/gen-syn-family.r")
            source("lib/07-calibrated-generator.r")
            source("lib/09-symmetric-generator.r")
            source("lib/11-projection-generator.r")
            rm(list = ".sobol.main", envir = globalenv())
        }
    })
on.exit(cleanup_parallel(cfg), add = TRUE)

# +------------------------------------------------------------+
# |  Main MC loop                                              |
# +------------------------------------------------------------+
cat(sprintf("=== MC ALL-METHOD: %d scenarios (%s, %d cores) ===\n\n",
            length(v.scenarios), cfg[["type"]], cfg[["n_cores"]]))

# Check results/mc/ exists
if (!dir.exists("results/mc")) {
    stop("results/mc/ not found. Run 26-mc-ml-generate.py first.")
}

l.all <- Map(function(sc) {
    s   <- l.scenarios[[sc]]
    d   <- s[["data"]]
    fam <- l.fam[[as.character(sc)]]
    lnk <- l.lnk[[as.character(sc)]]
    sup <- l.sup[[as.character(sc)]]

    mo.orig <- fit.model(d, fam, lnk)
    if (is.null(mo.orig)) return(list())

    # Discover seeds from Python-generated CSVs
    v.seeds <- discover.seeds(sc)
    s.n <- length(v.seeds)
    cat(sprintf("Scenario %d %s(%s): %d seeds found in results/mc/\n",
                sc, fam, lnk, s.n))

    if (s.n == 0L) return(list())

    l.res <- par_map(cfg, mc.one.rep,
                     v.seeds,
                     rep(list(d), s.n),
                     rep(sc, s.n),
                     rep(fam, s.n),
                     rep(lnk, s.n),
                     rep(list(sup), s.n),
                     rep(list(mo.orig), s.n),
                     .progress = TRUE)

    # Unwrap safe_fun wrappers
    l.res <- Map(function(r) {
        if (isTRUE(r[["error"]])) {
            cat(sprintf("  FAILED: %s\n", r[["message"]]))
            NULL
        } else r[["value"]]
    }, l.res) |> Filter(f = Negate(is.null))

    cat(sprintf("  %d/%d succeeded\n\n", length(l.res), s.n))
    l.res
}, v.scenarios) |> Reduce(f = c)

# +------------------------------------------------------------+
# |  Save results                                              |
# +------------------------------------------------------------+
saveRDS(l.all, "results/mc-allmethod-results.rds")

v.methods <- c("copula", "calibrated", "projection",
               "ctgan", "tvae", "gauscop", "copgan")

d.summary <- Map(function(sc) {
    l.sc <- Filter(function(r) r[["scenario"]] == sc, l.all)
    Map(function(m) {
        v.vals <- Map(function(r) r[[m]], l.sc) |>
            Filter(f = Negate(is.null))
        if (length(v.vals) == 0L) {
            return(data.frame(
                scenario = sc, method = m,
                maxpct_mean = NA, maxpct_sd = NA,
                se_mean = NA, se_sd = NA,
                sig_mean = NA, lt5_mean = NA, n_ok = 0L))
        }
        data.frame(
            scenario    = sc,
            method      = m,
            maxpct_mean = mean(mapply(`[[`, v.vals, "maxpct"), na.rm = TRUE),
            maxpct_sd   = sd(mapply(`[[`, v.vals, "maxpct"), na.rm = TRUE),
            se_mean     = mean(mapply(`[[`, v.vals, "se"), na.rm = TRUE),
            se_sd       = sd(mapply(`[[`, v.vals, "se"), na.rm = TRUE),
            sig_mean    = mean(mapply(`[[`, v.vals, "sig"), na.rm = TRUE),
            lt5_mean    = mean(mapply(`[[`, v.vals, "lt5"), na.rm = TRUE),
            n_ok        = length(v.vals))
    }, v.methods) |> Reduce(f = rbind)
}, v.scenarios) |> Reduce(f = rbind)

# A summary in which every replication failed is a broken run, not a
# result.  Stop loudly rather than write an all-NA table that exits 0.
if (all(d.summary[["n_ok"]] == 0L)) {
    stop("every replication failed on every scenario, refusing to ",
         "write an all-NA summary")
}

write.csv(d.summary, "results/mc-allmethod-summary.csv", row.names = FALSE)

cat("\n=== SUMMARY ===\n\n")
print(d.summary, digits = 3)
cat("\nSaved: mc-allmethod-results.rds, mc-allmethod-summary.csv\n")
