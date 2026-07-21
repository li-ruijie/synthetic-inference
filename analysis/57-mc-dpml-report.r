#!/usr/bin/env Rscript
# +------------------------------------------------------------+
# |  Monte Carlo diagnostics: DP-ML generators (parallel)      |
# +------------------------------------------------------------+
#
# Phase 2 of the DP-ML pipeline. Phase 1 (56-mc-dpml-generate.py)
# trains DP-SGD CTGAN and PATE-CTGAN at eps = 10 and writes synthetic
# CSVs to results/mc-dpml/. This script fits the true model on each
# synthetic dataset and computes the same diagnostics as
# 25-mc-allmethod.r (fit.model and diagnostics are verbatim copies
# from there), supporting the DP machine learning comparison sentence
# in the main text's DP simulation subsection.
#
# The run of record used R 4.6.1 on Windows (PSOCK backend).  A
# Linux fork run under R 4.5.0 reproduces it to 2.0e-12.
#
# Output:
#   results/mc-dpml-results.rds  -- full results list
#   results/mc-dpml-summary.csv  -- mean (sd) per method per scenario
#
# Usage:
#   nohup Rscript 57-mc-dpml-report.r > mc-dpml-report.log 2>&1 &
#
# +------------------------------------------------------------+
# |  Setup                                                     |
# +------------------------------------------------------------+
options(width = 128)

.shared.loaded <- TRUE
source("lib/00-shared.r")

s.wd <- tryCatch(
    normalizePath(dirname(sys.frame(1)[["ofile"]]), mustWork = TRUE),
    error = function(e) getwd())
setwd(s.wd)

source("lib/parallel_setup.R")

# +------------------------------------------------------------+
# |  Configuration                                             |
# +------------------------------------------------------------+
v.scenarios <- c(1L, 3L, 4L, 6L, 8L, 9L)
v.methods   <- c("dpctgan", "patectgan")

l.fam <- list(`1` = "gaussian", `3` = "gamma",       `4` = "quasi-gamma",
              `6` = "poisson",  `8` = "quasi-gamma", `9` = "gaussian")
l.lnk <- list(`1` = "identity", `3` = "log", `4` = "log",
              `6` = "log",      `8` = "log", `9` = "identity")

# +------------------------------------------------------------+
# |  Helpers (verbatim from 25-mc-allmethod.r)                 |
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

# Categorical columns are returned by the DP generators as plain
# integers. Coerce to factors with the original levels so that
# model.matrix builds the same dummy columns as the original fit.
align.types <- function(d.syn, d.orig) {
    Reduce(function(d, s.col) {
        if (is.factor(d.orig[[s.col]]) && !is.factor(d[[s.col]]))
            d <- `[[<-`(d, s.col,
                        value = factor(d[[s.col]],
                                       levels = levels(d.orig[[s.col]])))
        d
    }, intersect(colnames(d.syn), colnames(d.orig)), d.syn)
}

read.dpml.csv <- function(sc, method, seed) {
    f <- sprintf("results/mc-dpml/sc%02d_seed%d_%s.csv", sc, seed, method)
    if (file.exists(f)) {
        tryCatch(read.csv(f), error = function(e) NULL)
    } else NULL
}

discover.seeds <- function(sc) {
    v.files <- list.files("results/mc-dpml",
                          pattern = sprintf("^sc%02d_seed", sc))
    v.seeds <- sub(sprintf("^sc%02d_seed(\\d+)_.*\\.csv$", sc),
                   "\\1", v.files)
    sort(unique(as.integer(v.seeds)))
}

mc.one.rep <- function(seed, d, sc, fam, lnk, mo.orig) {
    l.diag <- Map(function(m) {
        d.syn <- read.dpml.csv(sc, m, seed)
        if (is.null(d.syn)) return(NULL)
        diagnostics(mo.orig, fit.model(align.types(d.syn, d), fam, lnk))
    }, v.methods)
    c(list(scenario = sc, seed = seed),
      `names<-`(l.diag, v.methods))
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
    export = c("fit.model", "diagnostics", "align.types", "read.dpml.csv",
               "v.methods"),
    init_worker = function() {
        Sys.setenv(OPENBLAS_NUM_THREADS = "1",
                   MKL_NUM_THREADS      = "1",
                   OMP_NUM_THREADS      = "1")
        # PSOCK workers start with an empty global environment, so they
        # source 00-shared.r the way the script head does.  The flag
        # skips the globalenv clear, and the guard makes this a no-op
        # on the fork and sequential backends.
        if (!exists(".shared.loaded", envir = globalenv())) {
            assign(".shared.loaded", TRUE, envir = globalenv())
            source("lib/00-shared.r")
        }
    })
on.exit(cleanup_parallel(cfg), add = TRUE)

# +------------------------------------------------------------+
# |  Main MC loop                                              |
# +------------------------------------------------------------+
cat(sprintf("=== MC DP-ML: %d scenarios (%s, %d cores) ===\n\n",
            length(v.scenarios), cfg[["type"]], cfg[["n_cores"]]))

if (!dir.exists("results/mc-dpml")) {
    stop("results/mc-dpml/ not found. Run 56-mc-dpml-generate.py first.")
}

l.all <- Map(function(sc) {
    s   <- l.scenarios[[sc]]
    d   <- s[["data"]]
    fam <- l.fam[[as.character(sc)]]
    lnk <- l.lnk[[as.character(sc)]]

    mo.orig <- fit.model(d, fam, lnk)
    if (is.null(mo.orig)) return(list())

    v.seeds <- discover.seeds(sc)
    s.n <- length(v.seeds)
    cat(sprintf("Scenario %d %s(%s): %d seeds found in results/mc-dpml/\n",
                sc, fam, lnk, s.n))

    if (s.n == 0L) return(list())

    l.res <- par_map(cfg, mc.one.rep,
                     v.seeds,
                     rep(list(d), s.n),
                     rep(sc, s.n),
                     rep(fam, s.n),
                     rep(lnk, s.n),
                     rep(list(mo.orig), s.n),
                     .progress = TRUE)

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
saveRDS(l.all, "results/mc-dpml-results.rds")

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
                sig_mean = NA, dir_mean = NA, lt5_mean = NA, n_ok = 0L))
        }
        data.frame(
            scenario    = sc,
            method      = m,
            maxpct_mean = mean(mapply(`[[`, v.vals, "maxpct"), na.rm = TRUE),
            maxpct_sd   = sd(mapply(`[[`, v.vals, "maxpct"), na.rm = TRUE),
            se_mean     = mean(mapply(`[[`, v.vals, "se"), na.rm = TRUE),
            se_sd       = sd(mapply(`[[`, v.vals, "se"), na.rm = TRUE),
            sig_mean    = mean(mapply(`[[`, v.vals, "sig"), na.rm = TRUE),
            dir_mean    = mean(mapply(`[[`, v.vals, "dir"), na.rm = TRUE),
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

write.csv(d.summary, "results/mc-dpml-summary.csv", row.names = FALSE)

cat("\n=== SUMMARY ===\n\n")
print(d.summary, digits = 3)
cat("\nSaved: results/mc-dpml-results.rds, results/mc-dpml-summary.csv\n")
