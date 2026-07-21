#!/usr/bin/env Rscript
# +--------------------------------------------------------------------------+
# | 68-kwise-fidelity.r -- k-wise distributional fidelity of the projection  |
# +--------------------------------------------------------------------------+
#
# Companion instrument. The coupling bound states that the
# projection is a coupling of the projected sample to the base copula draw,
# so every k-margin's empirical W1 distance is bounded by the mean absolute
# per-observation displacement.  (The inheritance result is a separate
# Corollary, that the projected generator keeps the base generator's k-wise
# distributional fidelity.  It is not what this script measures.)  This
# script measures both sides of that inequality on every draw and stops if
# the bound is violated (a violation is an implementation bug, not noise),
# and reports energy distances (Szekely and Rizzo) of base and projected
# draws against the original data, per variable pair and on the full joint.
#
# Displacement convention: d_i = y_released,i - y_base,i, the total
# displacement including Poisson rounding, since project.glm() does not
# expose the pre-rounding response. W1 <= mean |d_i| holds for the total
# displacement by the same identity coupling.
#
# Output: kwise-fidelity.csv, one row per (design, seed). Smoke mode:
# KWISE_SMOKE=1 restricts to Scenario 1 with B = 3.
#
# Replication-tree copy (flat layout). The PHI design of the analysis tree
# is skipped here unless the PHI data file sits beside this script, since
# the published tree does not carry it.

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

# +--------------------------------------------------------------------------+
# | Metric helpers                                                           |
# +--------------------------------------------------------------------------+

# Exact empirical W1 between two univariate samples of equal size: the
# optimal coupling on the line is the sorted matching.
w1.1d <- function(a, b) mean(abs(sort(a) - sort(b)))

# Standardise a data frame onto the numeric indicator scale of a reference
# frame: continuous columns are centred and scaled by the reference mean and
# sd, and each categorical column expands to its indicator columns at the
# Reference levels. Returns a numeric matrix. Column grouping by source
# variable is carried in attr "var".
std.rows <- function(d, d.ref) {
    l.cols <- Map(function(nm) {
        v <- d[[nm]]
        v.ref <- d.ref[[nm]]
        # Categorical means factor or character only. A numeric binary
        # column stays on the continuous scale: scenario 11's projected
        # response is fractional in [0, 1] while the original is {0, 1},
        # and an indicator encoding would collapse every fractional value
        # to zero and fake a perfect margin.
        if(is.factor(v.ref) || is.character(v.ref)) {
            v.lev <- `if`(is.factor(v.ref), levels(v.ref),
                          sort(unique(v.ref)))
            m <- sapply(v.lev[-1], function(lv) as.numeric(v == lv))
            m <- matrix(m, nrow = length(v))
            colnames(m) <- paste0(nm, v.lev[-1])
            m
        } else {
            m <- matrix((v - mean(v.ref)) / sd(v.ref), ncol = 1)
            colnames(m) <- nm
            m
        }
    }, names(d.ref))
    m.out <- do.call(cbind, l.cols)
    attr(m.out, "var") <- rep(names(d.ref),
                              mapply(function(m) ncol(m), l.cols))
    m.out
}

# Energy distance between two samples given as numeric matrices:
#   E = 2 mean ||x - y|| - mean ||x - x'|| - mean ||y - y'||
# with the within-sample means taken over all ordered pairs (diagonal zero),
# so mean ||x - x'|| = mean(dist(x)) * (n - 1) / n.
energy.dist <- function(m.x, m.y) {
    s.nx <- nrow(m.x)
    s.ny <- nrow(m.y)
    m.cross <- sqrt(pmax(0, outer(rowSums(m.x^2), rowSums(m.y^2), `+`) -
                            2 * tcrossprod(m.x, m.y)))
    2 * mean(m.cross) -
        mean(dist(m.x)) * (s.nx - 1) / s.nx -
        mean(dist(m.y)) * (s.ny - 1) / s.ny
}

# Energy distances of a synthetic frame against the original, on the
# standardised indicator scale: the full joint, the mean over variable
# pairs involving the response, and the mean over pairs not involving it
# (the X block, identical between the base and projected arms).
pair.energy <- function(d.syn, d.orig, y.nm) {
    m.o <- std.rows(d.orig, d.orig)
    m.s <- std.rows(d.syn,  d.orig)
    v.var <- attr(m.o, "var")
    v.nms <- unique(v.var)
    d.pairs <- t(combn(v.nms, 2))
    v.e <- apply(d.pairs, 1, function(pr) {
        v.keep <- v.var %in% pr
        energy.dist(m.s[, v.keep, drop = FALSE],
                    m.o[, v.keep, drop = FALSE])
    })
    b.resp <- apply(d.pairs, 1, function(pr) y.nm %in% pr)
    c(joint = energy.dist(m.s, m.o),
      rpair = mean(v.e[b.resp]),
      xpair = mean(v.e[!b.resp]))
}

# +--------------------------------------------------------------------------+
# | Self-checks (hand-computable cases, stop on any failure)                 |
# +--------------------------------------------------------------------------+

stopifnot(abs(w1.1d(c(0, 1), c(3, 0)) - 1) < 1e-12)
stopifnot(abs(w1.1d(c(5, 2, 9), c(2, 5, 9))) < 1e-12)
d.chk <- data.frame(a = c(0, 1, 2, 3), g = factor(c("u", "u", "v", "v")))
m.chk <- std.rows(d.chk, d.chk)
stopifnot(identical(dim(m.chk), c(4L, 2L)),
          abs(mean(m.chk[, "a"])) < 1e-12,
          abs(sd(m.chk[, "a"]) - 1) < 1e-12,
          identical(as.numeric(m.chk[, "gv"]), c(0, 0, 1, 1)))
m.e1 <- matrix(c(0, 0), ncol = 1)
m.e2 <- matrix(c(1, 1), ncol = 1)
# E = 2 * 1 - 0 - 0 = 2 for two point masses at distance 1
stopifnot(abs(energy.dist(m.e1, m.e2) - 2) < 1e-12)
stopifnot(abs(energy.dist(m.e1, m.e1)) < 1e-12)
cat("self-checks passed\n")

# +--------------------------------------------------------------------------+
# | One replication: reproduce the generator's own coupling                  |
# +--------------------------------------------------------------------------+

# Reproduces gen.syn.projected() internals (11-projection-generator.r):
# gen.syn.family() then project.glm(), same defaults, same seed, so the RNG
# stream matches the shipped generator exactly. Asserted per design in
# run.design() below.
mc.one.fid <- function(d, y.nm, family, link, sup, seed) {
    d.c <- gen.syn.family(d, order = 4, supports = sup, shrink = TRUE,
                          seed = seed)
    v.x.nms <- setdiff(names(d), y.nm)
    v.y.proj <- project.glm(d.c, d, y.nm, v.x.nms, family, link)
    d.p <- d.c
    d.p[[y.nm]] <- v.y.proj

    s.sd.y <- sd(d[[y.nm]])
    v.d <- (v.y.proj - d.c[[y.nm]]) / s.sd.y
    s.w1 <- w1.1d(v.y.proj, d.c[[y.nm]]) / s.sd.y

    # The guard: the coupling inequality holds deterministically for every
    # draw. A violation is an implementation bug and stops the run.
    stopifnot(s.w1 <= mean(abs(v.d)) + 1e-9)

    v.e.b <- pair.energy(d.c, d, y.nm)
    v.e.p <- pair.energy(d.p, d, y.nm)

    data.frame(n = nrow(d), mean.disp = mean(abs(v.d)), w1.resp = s.w1,
               e.joint.base = v.e.b[["joint"]],
               e.joint.proj = v.e.p[["joint"]],
               e.rpair.base = v.e.b[["rpair"]],
               e.rpair.proj = v.e.p[["rpair"]],
               e.xpair      = v.e.b[["xpair"]])
}

# +--------------------------------------------------------------------------+
# | Designs                                                                  |
# +--------------------------------------------------------------------------+

load("results/d.sim.scenarios.RData")

# Simulation designs: scenario id, family and link tokens as accepted by
# project.glm(), and the support vectors of the analysis tree's
# 18-mc-projection.r l.sup, plus scenario 11 per 60-mc-binary.r.
l.designs <- list(
    list(nm = "sc1",  i = 1L,  fam = "gaussian",    lnk = "identity",
         sup = c(response = "real",     age = "real", income = "real")),
    list(nm = "sc3",  i = 3L,  fam = "gamma",       lnk = "log",
         sup = c(response = "positive", age = "real", income = "positive")),
    list(nm = "sc4",  i = 4L,  fam = "quasi-gamma", lnk = "log",
         sup = c(response = "zeroinfl", age = "real", income = "positive")),
    list(nm = "sc6",  i = 6L,  fam = "poisson",     lnk = "log",
         sup = c(response = "count",    age = "real", income = "positive")),
    list(nm = "sc8",  i = 8L,  fam = "quasi-gamma", lnk = "log",
         sup = c(response = "zeroinfl", age = "real", income = "positive")),
    list(nm = "sc9",  i = 9L,  fam = "gaussian",    lnk = "identity",
         sup = NULL),
    list(nm = "sc11", i = 11L, fam = "binomial",    lnk = "logit",
         sup = c(age = "real", income = "positive")))

b.smoke <- nzchar(Sys.getenv("KWISE_SMOKE"))
s.B <- `if`(b.smoke, 3L, 100L)
if(b.smoke) l.designs <- l.designs[1]
v.seeds <- 6489 + seq_len(s.B)

# +--------------------------------------------------------------------------+
# | Run                                                                      |
# +--------------------------------------------------------------------------+

run.design <- function(nm, d, y.nm, fam, lnk, sup) {
    # Coupling-fidelity guard against the shipped generator: the paired
    # reconstruction must equal gen.syn.projected() bit for bit at the
    # first seed (catches any future drift in the generator's internals).
    d.direct <- gen.syn.projected(d, y.nm = y.nm, family = fam, link = lnk,
                                  supports = sup, seed = v.seeds[1])
    d.c1 <- gen.syn.family(d, order = 4, supports = sup, shrink = TRUE,
                           seed = v.seeds[1])
    v.y1 <- project.glm(d.c1, d, y.nm, setdiff(names(d), y.nm), fam, lnk)
    d.p1 <- d.c1
    d.p1[[y.nm]] <- v.y1
    stopifnot(isTRUE(all.equal(d.direct, d.p1)))

    l.rows <- Map(function(seed) {
        d.row <- tryCatch(mc.one.fid(d, y.nm, fam, lnk, sup, seed),
                          error = function(e) {
                              # A guard violation must stop the run: only
                              # generator-side failures (rare copula draws)
                              # are skippable.
                              if(grepl("s.w1 <=", conditionMessage(e)))
                                  stop(e)
                              NULL
                          })
        if(!is.null(d.row)) { d.row[["design"]] <- nm
                              d.row[["seed"]]   <- seed }
        d.row
    }, v.seeds)
    do.call(rbind, Filter(Negate(is.null), l.rows))
}

l.out <- Map(function(cfg) {
    s <- l.scenarios[[cfg[["i"]]]]
    cat(sprintf("design %s (%s/%s) ...\n", cfg[["nm"]],
                cfg[["fam"]], cfg[["lnk"]]))
    run.design(cfg[["nm"]], s[["data"]], "response",
               cfg[["fam"]], cfg[["lnk"]], cfg[["sup"]])
}, l.designs)

# PHI waist application: runs only when the data file sits beside this
# script. The published tree does not carry it, so the design is skipped
# there with a notice, and the simulation designs above stay byte-identical
# to the analysis tree's run.
if(!b.smoke && file.exists("phi_phase1_baseline.dta")) {
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
    cat("design phi (gaussian/identity) ...\n")
    l.out <- c(l.out, list(run.design("phi", d.app, "waist", "gaussian",
                                      "identity",
                                      c(waist = "real", age = "real"))))
} else if(!b.smoke) {
    cat("design phi skipped: phi_phase1_baseline.dta not in this tree\n")
}

d.out <- do.call(rbind, l.out)
d.out <- d.out[, c("design", "seed", "n", "mean.disp", "w1.resp",
                   "e.joint.base", "e.joint.proj", "e.rpair.base",
                   "e.rpair.proj", "e.xpair")]
write.csv(d.out, "results/kwise-fidelity.csv", row.names = FALSE)

# +--------------------------------------------------------------------------+
# | Summary table and LaTeX body                                             |
# +--------------------------------------------------------------------------+

agg <- function(v) sprintf("%.4f (%.4f)", mean(v), sd(v))
cat("\nkwise fidelity, mean (sd) over", s.B, "replications\n")
cat(sprintf("%-6s %6s %18s %18s %18s %18s %22s\n", "design", "n",
            "mean.disp", "w1.resp", "e.joint base", "e.joint proj",
            "e.rpair base -> proj"))
for(nm in unique(d.out[["design"]])) {
    d.g <- d.out[d.out[["design"]] == nm, ]
    cat(sprintf("%-6s %6d %18s %18s %18s %18s %9s -> %9s\n", nm,
                d.g[["n"]][1], agg(d.g[["mean.disp"]]),
                agg(d.g[["w1.resp"]]), agg(d.g[["e.joint.base"]]),
                agg(d.g[["e.joint.proj"]]),
                sprintf("%.4f", mean(d.g[["e.rpair.base"]])),
                sprintf("%.4f", mean(d.g[["e.rpair.proj"]]))))
}
cat("\nLaTeX body (booktabs, copy verbatim into the supplement table):\n")
for(nm in unique(d.out[["design"]])) {
    d.g <- d.out[d.out[["design"]] == nm, ]
    cat(sprintf("%s & %d & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\\n",
                nm, d.g[["n"]][1], mean(d.g[["mean.disp"]]),
                mean(d.g[["w1.resp"]]), mean(d.g[["e.joint.base"]]),
                mean(d.g[["e.joint.proj"]]), mean(d.g[["e.rpair.base"]]),
                mean(d.g[["e.rpair.proj"]])))
}
cat(sprintf("\nguard: W1 <= mean displacement held on every one of %d draws\n",
            nrow(d.out)))
