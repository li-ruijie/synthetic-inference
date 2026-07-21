#!/usr/bin/env r
# +------------------------------------------------------------------------+
# | 59-mc-cross-family.r -- Cross-family/link robustness of exact refit     |
# +------------------------------------------------------------------------+
#
# Candidate P1 of the paper's replication record. Theorem 7
# (11-projection-generator.r, and the manuscript) proves that refitting the
# Same family/link the generator matched reproduces beta_hat exactly, over
# every predictor submodel. This script tests the untested case: what a
# Different analyst family/link recovers from the same release.
#
# Two release modes (see the paper's replication record):
#
#   synthetic-X: the full generator. Covariates are synthesised too, so the
#     projection reconstructs beta_hat_A on a new design X_syn. The same-family
#     refit stays exact, but a different analyst family reads y* through a
#     different estimating equation and incurs a deviation. Empirically this
#     deviation is modest for a canonical-link generator (its projection
#     reconstructs the actual conditional mean g_A^{-1}(X beta_hat_A), which a
#     neighbouring family also fits well) and large for a non-canonical log-link
#     generator (only a variance-weighted mean is reconstructed).
#
#   public-X: covariates are released unchanged (X = X_orig) and only the
#     response is synthesised. A canonical-A projection then enforces
#     X' y* = X' y on the original design, so any canonical analyst family B
#     solves X' y = X' mu_B(X beta) exactly and refits with 0.0% deviation
#     (Proposition, public-covariate cross-family invariance). Non-canonical B
#     still deviate, because their score uses X' W_B y, not X' y.
#
# Method: Monte Carlo protocol of 18-mc-projection.r (B=100 seeds 6489+1..100).
# The generator draw depends only on (scenario, seed), not the analyst family,
# so y*_A is generated once per (scenario, seed) and the whole analyst panel is
# refitted on it. Output: mc-cross-family.csv plus a console report.
# Deterministic under OPENBLAS_NUM_THREADS=1.
#
# +------------------------------------------------------------------------+
# | Setup (mirror 18-mc-projection.r: dotted .sobol.main survives the       |
# | globalenv clear in 00-shared.r:103, suppressing the generator harness)  |
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
# | Analyst panel (family B, link B)                                        |
# +------------------------------------------------------------------------+
# The families a downstream analyst might fit on the release. "pos" flags a
# positive-support family (log/inverse link or Gamma): rows with y <= 0 are
# dropped before the fit, exactly as a real analyst would. "canon" is whether
# (family, link) is the canonical pair (drives the public-X exactness split).
l.analyst <- list(
    list(tag = "gaussian(identity)", family = "gaussian", link = "identity", pos = FALSE, canon = TRUE),
    list(tag = "gaussian(log)",      family = "gaussian", link = "log",      pos = TRUE,  canon = FALSE),
    list(tag = "poisson(log)",       family = "poisson",  link = "log",      pos = FALSE, canon = TRUE),
    list(tag = "poisson(identity)",  family = "poisson",  link = "identity", pos = FALSE, canon = FALSE),
    list(tag = "gamma(log)",         family = "gamma",    link = "log",      pos = TRUE,  canon = FALSE),
    list(tag = "gamma(inverse)",     family = "gamma",    link = "inverse",  pos = TRUE,  canon = TRUE))

# fit.analyst: robust GLM/LM fit under an arbitrary analyst spec. Applies the
# same support filtering to original and synthetic data so the two fits are
# comparable. Tries default glm initialisation first, then an eta-based start
# (intercept = g(y_bar), slopes 0), mirroring project.glm's fallback. Returns
# NULL when the fit is infeasible (e.g. positive family on all-zero rows) or
# fails to converge, which the caller records as an infeasible cell.
fit.analyst <- function(d, y.nm, x.nms, spec) {
    fml <- reformulate(x.nms, response = y.nm)
    v.y <- d[[y.nm]]
    d.fit <- `if`(spec[["pos"]] && any(v.y <= 0), d[v.y > 0, , drop = FALSE], d)
    if(nrow(d.fit) < length(x.nms) + 2) return(NULL)
    if(spec[["family"]] == "gaussian" && spec[["link"]] == "identity")
        return(tryCatch(lm(fml, d.fit), error = function(e) NULL))
    fam.obj <- make.glm.family(spec[["family"]], spec[["link"]])
    suppressWarnings(tryCatch(
        glm(fml, d.fit, family = fam.obj, control = glm.control(maxit = 200)),
        error = function(e) tryCatch({
            s.np   <- ncol(model.matrix(fml, d.fit))
            s.eta0 <- fam.obj[["linkfun"]](max(mean(d.fit[[y.nm]]), 1e-3))
            glm(fml, d.fit, family = fam.obj,
                start = c(s.eta0, rep(0, s.np - 1)),
                control = glm.control(maxit = 200))
        }, error = function(e) NULL)))
}

# compare.fit: deviation of a synthetic-data refit from the original-data refit
# under the same analyst family. Aligned by coefficient name so a bootstrapped
# design that loses a rare factor level (rank deficiency) leaves an NA rather
# than misaligning rows (same guard as 18-mc-projection.r:mc.one).
compare.fit <- function(mo.o, mo.s) {
    if(is.null(mo.o) || is.null(mo.s)) return(NULL)
    v.o   <- coef(mo.o)
    v.s   <- coef(mo.s)[names(v.o)]
    s.o   <- summary(mo.o)[["coefficients"]]
    s.s   <- summary(mo.s)[["coefficients"]]
    v.se.o <- s.o[, 2][names(v.o)]
    v.se.s <- s.s[, 2][names(v.o)]
    v.p4.o <- s.o[, 4][names(v.o)]
    v.p4.s <- s.s[, 4][names(v.o)]
    v.pct    <- abs(100 * (v.s - v.o) / abs(v.o))
    v.se.pct <- abs(100 * (v.se.s - v.se.o) / abs(v.se.o))
    list(maxpct   = max(v.pct, na.rm = TRUE),
         se       = max(v.se.pct, na.rm = TRUE),
         sig      = mean((v.p4.o < 0.05) == (v.p4.s < 0.05), na.rm = TRUE),
         dir      = mean(sign(v.o) == sign(v.s), na.rm = TRUE),
         lt5      = mean(v.pct < 5, na.rm = TRUE),
         fullrank = !anyNA(v.s))
}

# gen.public.projected: the "public covariate" release. Covariates are kept
# unchanged (X = X_orig).  Only the response is synthesised. A permuted response
# serves as the pre-projection draw (its role is residual novelty only -- the
# projected coefficients are invariant to it), then project.glm imposes the
# family-A score on the original design. X being original, a canonical-A
# projection enforces X' y* = X' y exactly.
gen.public.projected <- function(d, y.nm, x.nms, family, link, seed) {
    set.seed(seed)
    d.cop <- d
    d.cop[[y.nm]] <- sample(d[[y.nm]])
    d.out <- d
    d.out[[y.nm]] <- project.glm(d.cop, d, y.nm, x.nms, family, link)
    d.out
}

# +------------------------------------------------------------------------+
# | Run Monte Carlo                                                         |
# +------------------------------------------------------------------------+
load("results/d.sim.scenarios.RData")
s.B <- 100

# Generator scenarios (family A) and their support declarations, matching the
# rows of 18-mc-projection.r. Sc 1 and 6 are canonical.  Sc 3 and 4 are not.
v.gen       <- c(1, 3, 4, 6)
v.canon.gen <- c(1, 6)
l.sup <- list(
    `1` = c(response = "real",     age = "real", income = "real"),
    `3` = c(response = "positive", age = "real", income = "positive"),
    `4` = c(response = "zeroinfl", age = "real", income = "positive"),
    `6` = c(response = "count",    age = "real", income = "positive"))

# run.scenario: for each seed, produce the release under `mode`, refit the whole
# analyst panel on it, and aggregate the per-analyst deviation over seeds.
run.scenario <- function(i, mode = "synthetic") {
    s     <- l.scenarios[[i]]
    d     <- s[["data"]]
    y.nm  <- "response"
    x.nms <- setdiff(names(d), y.nm)
    a.fam <- s[["family"]]
    a.lnk <- s[["link"]]
    sup   <- l.sup[[as.character(i)]]

    # original-data fit under each analyst family (the target of comparison)
    l.mo.o <- lapply(l.analyst, function(sp) fit.analyst(d, y.nm, x.nms, sp))

    v.seeds <- 6489 + seq_len(s.B)
    l.per.seed <- Map(function(seed) {
        d.syn <- tryCatch(
            `if`(mode == "public",
                 gen.public.projected(d, y.nm, x.nms, a.fam, a.lnk, seed),
                 gen.syn.projected(d, y.nm = y.nm, family = a.fam, link = a.lnk,
                                   supports = sup, seed = seed)),
            error = function(e) NULL)
        if(is.null(d.syn)) return(NULL)
        lapply(seq_along(l.analyst), function(j)
            compare.fit(l.mo.o[[j]],
                        fit.analyst(d.syn, y.nm, x.nms, l.analyst[[j]])))
    }, v.seeds)
    l.per.seed <- Filter(Negate(is.null), l.per.seed)

    # aggregate per analyst family j
    l.rows <- lapply(seq_along(l.analyst), function(j) {
        sp   <- l.analyst[[j]]
        l.rj <- Filter(Negate(is.null), lapply(l.per.seed, function(x) x[[j]]))
        s.n  <- length(l.rj)
        d.base <- data.frame(
            mode          = mode,
            scenario      = i,
            gen           = sprintf("%s(%s)", a.fam, a.lnk),
            gen_canon     = is.canonical(a.fam, a.lnk),
            analyst       = sp[["tag"]],
            analyst_canon = sp[["canon"]],
            same_model    = sp[["tag"]] == sprintf("%s(%s)", a.fam, a.lnk),
            n_feasible    = s.n,
            stringsAsFactors = FALSE)
        if(s.n == 0)
            return(cbind(d.base, maxpct_mean = NA, maxpct_sd = NA,
                         se_mean = NA, se_sd = NA, sig = NA, dir = NA, lt5 = NA))
        v.maxpct <- vapply(l.rj, `[[`, numeric(1), "maxpct")
        v.se     <- vapply(l.rj, `[[`, numeric(1), "se")
        v.sig    <- vapply(l.rj, `[[`, numeric(1), "sig")
        v.dir    <- vapply(l.rj, `[[`, numeric(1), "dir")
        v.lt5    <- vapply(l.rj, `[[`, numeric(1), "lt5")
        cbind(d.base,
              maxpct_mean = mean(v.maxpct), maxpct_sd = sd(v.maxpct),
              se_mean = mean(v.se), se_sd = sd(v.se),
              sig = 100 * mean(v.sig), dir = 100 * mean(v.dir),
              lt5 = 100 * mean(v.lt5))
    })
    do.call(rbind, l.rows)
}

# print.grid: console table of one mode's rows, grouped by generator scenario.
print.grid <- function(d.all, v.scen, s.title) {
    cat(s.title, "\n")
    cat(strrep("=", 92), "\n\n")
    Map(function(i) {
        d.i <- d.all[d.all[["scenario"]] == i, ]
        cat(sprintf("Generator: Scenario %d %s  [%s]\n", i, d.i[["gen"]][1],
                    `if`(d.i[["gen_canon"]][1], "canonical", "non-canonical")))
        cat(sprintf("  %-20s %-9s %13s %13s %7s %5s\n",
                    "analyst refit", "canon", "max|%| (sd)", "SE% (sd)", "sig%", "n"))
        cat("  ", strrep("-", 74), "\n", sep = "")
        Map(function(k) {
            cat(sprintf("  %-20s %-9s %6.2f (%5.2f) %6.2f (%5.2f) %6.1f %4d%s\n",
                        d.i[["analyst"]][k],
                        `if`(d.i[["analyst_canon"]][k], "canon", "non"),
                        d.i[["maxpct_mean"]][k], d.i[["maxpct_sd"]][k],
                        d.i[["se_mean"]][k], d.i[["se_sd"]][k],
                        d.i[["sig"]][k], d.i[["n_feasible"]][k],
                        `if`(isTRUE(d.i[["same_model"]][k]), "  <- same as generator", "")))
        }, seq_len(nrow(d.i)))
        cat("\n")
    }, v.scen)
    invisible(NULL)
}

d.syn.grid <- do.call(rbind, lapply(v.gen,       function(i) run.scenario(i, "synthetic")))
d.pub.grid <- do.call(rbind, lapply(v.canon.gen, function(i) run.scenario(i, "public")))
d.all <- rbind(d.syn.grid, d.pub.grid)
write.csv(d.all, "results/mc-cross-family.csv", row.names = FALSE)

print.grid(d.syn.grid, v.gen,
           sprintf("MODE synthetic-X: full generator (X synthesised) -- %d reps", s.B))
print.grid(d.pub.grid, v.canon.gen,
           sprintf("MODE public-X: covariates released unchanged -- %d reps", s.B))

cat("Wrote mc-cross-family.csv\n")
# +------------------------------------------------------------------------+
# | End of 59-mc-cross-family.r                                            |
# +------------------------------------------------------------------------+
