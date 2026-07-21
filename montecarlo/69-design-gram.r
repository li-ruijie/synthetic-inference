#!/usr/bin/env r
# +--------------------------------------------------------------------------+
# | 69-design-gram.r -- design-Gram projection under full synthesis          |
# +--------------------------------------------------------------------------+
#
# Companion instrument to the design-Gram proposition and the two-part
# full-synthesis corollary. Two paired arms per (design, seed): "free" is the
# shipped full-synthesis generator (design.exact = FALSE), "pinned" adds the
# design-Gram projection. The free arm reproduces the section-5 SE drift.  The
# pinned arm's refit triple is exact for OLS channels, within the fixed-point
# tolerance for GLM channels, and per block for the two-part channels.
#
# The guard: on every pinned draw the driver re-asserts the Gram identity
# at analysis time (raw Gram for OLS and two-part channels, Fisher-weighted
# Gram for GLM channels), on top of the generator's own postconditions. A
# violation is an implementation bug, not noise, and stops the run.
#
# Columns of design-gram.csv, one row per (design, seed, arm):
#   coef.dev, se.dev   max |%| deviation of full-model coefficients / SEs,
#                      aligned by name (two-part: the positive block)
#   sub.coef.dev, sub.se.dev  the same on a per-seed random strict submodel
#                      (OLS channels only, design-Gram clause (a), NA otherwise)
#   z.coef.dev         zero-block max absolute coefficient deviation
#                      (two-part channels only, the two-part corollary, NA
#                      else.  Absolute since the zero model's true coefficients
#                      include near-zeros and a relative metric explodes)
#   viol.rate, viol.max  support violations of continuous predictors
#                      against their token bounds (rate in %, max distance)
#   e.joint            energy distance to the original on the standardised
#                      indicator scale (full joint)
#   full.rank          refit design full rank (PHI only, NA otherwise)
#
# Smoke mode: DGRAM_SMOKE=1 restricts to Scenario 1 with B = 3.
# Output tables: the summary table and a LaTeX
# table body are printed at the end of the run.
#
# Replication-tree copy (flat layout). The PHI design of the analysis
# tree is skipped here unless the PHI data file sits beside this script,
# since the published tree does not carry it.

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
# | Metric helpers (std.rows and energy.dist copied from 68-kwise-fidelity.r |
# | which is a driver, not a library.  Sourcing it would run its MC)        |
# +--------------------------------------------------------------------------+

std.rows <- function(d, d.ref) {
    l.cols <- Map(function(nm) {
        v <- d[[nm]]
        v.ref <- d.ref[[nm]]
        if(is.factor(v.ref) || is.character(v.ref)) {
            v.lev <- `if`(is.factor(v.ref), levels(v.ref),
                          sort(unique(v.ref)))
            m <- mapply(function(lv) as.numeric(v == lv), v.lev[-1])
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

energy.dist <- function(m.x, m.y) {
    s.nx <- nrow(m.x)
    s.ny <- nrow(m.y)
    m.cross <- sqrt(pmax(0, outer(rowSums(m.x^2), rowSums(m.y^2), `+`) -
                            2 * tcrossprod(m.x, m.y)))
    2 * mean(m.cross) -
        mean(dist(m.x)) * (s.nx - 1) / s.nx -
        mean(dist(m.y)) * (s.ny - 1) / s.ny
}

# Max |%| deviation between two named coefficient vectors, aligned by name,
# NA-dropped (rank-deficient PHI draws lose inestimable levels).
dev.pct <- function(v.s, v.o) {
    v.common <- intersect(names(v.s)[!is.na(v.s)], names(v.o)[!is.na(v.o)])
    max(abs((v.s[v.common] - v.o[v.common]) / v.o[v.common])) * 100
}

# Support-violation rate (%) and max distance of the continuous predictors
# against their token bounds. Token per predictor from detect.support() of
# the original column: positive -> (0, Inf), zeroinfl -> [0, Inf),
# bounded -> [min, max] of the original, real / count -> unbounded here.
support.viol <- function(d.syn, d.orig, v.x.cnt) {
    l.v <- Map(function(nm) {
        s.tok <- detect.support(d.orig[[nm]])
        v <- d.syn[[nm]]
        v.dist <- `if`(s.tok %in% c("positive", "zeroinfl"),
                       pmax(0, -v),
                       `if`(s.tok == "bounded",
                            pmax(0, min(d.orig[[nm]]) - v) +
                            pmax(0, v - max(d.orig[[nm]])),
                            rep(0, length(v))))
        v.dist
    }, v.x.cnt)
    v.all <- unlist(l.v)
    c(rate = mean(v.all > 0) * 100, max = max(v.all, 0))
}

# +--------------------------------------------------------------------------+
# | Refits                                                                   |
# +--------------------------------------------------------------------------+

# Family-appropriate full-model refit.  Returns list(coef, se) or NULL.
# Quasi families use the project.glm() start idiom.
refit <- function(d, y.nm, fam, lnk) {
    fml <- reformulate(setdiff(names(d), y.nm), response = y.nm)
    mo <- tryCatch({
        if(fam == "gaussian") {
            lm(fml, d)
        } else if(fam == "quasi-gamma") {
            s.np <- ncol(model.matrix(fml, d))
            glm(fml, d[d[[y.nm]] > 0, ],
                family = quasi(link = "log", variance = "mu^2"),
                start = c(1, rep(0, s.np - 1)),
                control = glm.control(maxit = 200))
        } else if(fam == "gamma") {
            glm(fml, d, family = Gamma(link = "log"),
                control = glm.control(maxit = 200))
        } else if(fam == "poisson") {
            glm(fml, d, family = poisson(),
                control = glm.control(maxit = 200))
        } else stop("refit: unknown family ", fam)
    }, error = function(e) NULL)
    if(is.null(mo)) return(NULL)
    m.s <- coef(summary(mo))
    list(coef = coef(mo), se = `names<-`(m.s[, 2], rownames(m.s)))
}

refit.zero <- function(d, y.nm) {
    fml <- reformulate(setdiff(names(d), y.nm),
                       response = sprintf("I(%s > 0)", y.nm))
    mo <- glm(fml, d, family = binomial(),
              control = glm.control(maxit = 200))
    coef(mo)
}

# +--------------------------------------------------------------------------+
# | One replication, one arm                                                 |
# +--------------------------------------------------------------------------+

mc.one.dg <- function(d, y.nm, fam, lnk, sup, seed, b.pinned, b.phi) {
    b.twopart <- fam == "quasi-gamma"
    d.s <- `if`(b.twopart && b.pinned,
        gen.syn.projected.twopart(d, y.nm = y.nm, family = fam, link = lnk,
                                  supports = sup, seed = seed,
                                  design.exact = TRUE),
        gen.syn.projected(d, y.nm = y.nm, family = fam, link = lnk,
                          supports = sup, seed = seed,
                          design.exact = b.pinned))

    v.x.nms <- setdiff(names(d), y.nm)
    v.types <- mapply(id.cat, d[v.x.nms])
    v.x.cnt <- names(v.types[v.types == "cnt"])
    fml.x   <- reformulate(v.x.nms)

    # The guard (pinned arm): the identity of the design-Gram proposition,
    # re-checked at analysis time. Raw Gram for OLS and two-part channels (and
    # the positive block for two-part), Fisher-weighted Gram for GLM channels.
    if(b.pinned) {
        m.a <- model.matrix(fml.x, d)
        m.b <- model.matrix(fml.x, d.s)
        s.scl <- max(abs(crossprod(m.a)))
        if(fam %in% c("gaussian", "quasi-gamma")) {
            stopifnot(max(abs(crossprod(m.b) - crossprod(m.a)))
                      < 1e-7 * s.scl)
            if(b.twopart) {
                m.a.p <- m.a[d[[y.nm]] > 0, , drop = FALSE]
                m.b.p <- m.b[d.s[[y.nm]] > 0, , drop = FALSE]
                stopifnot(max(abs(crossprod(m.b.p) - crossprod(m.a.p)))
                          < 1e-7 * s.scl)
            }
        } else {
            fml.y <- reformulate(v.x.nms, response = y.nm)
            fam.obj <- `if`(fam == "gamma", Gamma(link = "log"), poisson())
            mo.w <- glm(fml.y, d, family = fam.obj,
                        control = glm.control(maxit = 200))
            w.chk <- function(m.X) {
                v.eta <- as.numeric(m.X %*% coef(mo.w))
                fam.obj[["mu.eta"]](v.eta)^2 /
                    fam.obj[["variance"]](fam.obj[["linkinv"]](v.eta))
            }
            m.I.a <- crossprod(m.a * sqrt(w.chk(m.a)))
            m.I.b <- crossprod(m.b * sqrt(w.chk(m.b)))
            stopifnot(max(abs(m.I.b - m.I.a)) < 1e-6 * max(abs(m.I.a)))
        }
    }

    l.o <- refit(d,   y.nm, fam, lnk)
    l.s <- refit(d.s, y.nm, fam, lnk)
    if(is.null(l.o) || is.null(l.s)) return(NULL)
    s.coef <- dev.pct(l.s[["coef"]], l.o[["coef"]])
    s.se   <- dev.pct(l.s[["se"]],   l.o[["se"]])

    # full rank means the refit reproduces the original coefficient name
    # set with no NA. A dropped factor level shows as a missing column,
    # not an NA (the synthetic factor's levels attribute keeps only the
    # drawn levels), so an any-NA check alone is vacuously true.
    b.fr <- !any(is.na(l.s[["coef"]])) &&
            setequal(names(l.s[["coef"]]), names(l.o[["coef"]]))

    # zero block for the two-part channels. Absolute deviation: the zero
    # model's true coefficients include near-zeros by design, so a
    # relative metric explodes on them. coef.dev stays the positive
    # block's own relative measure.
    s.z <- NA_real_
    if(b.twopart) {
        v.z.o <- refit.zero(d,   y.nm)
        v.z.s <- refit.zero(d.s, y.nm)
        v.common <- intersect(names(v.z.s)[!is.na(v.z.s)],
                              names(v.z.o)[!is.na(v.z.o)])
        s.z <- max(abs(v.z.s[v.common] - v.z.o[v.common]))
    }

    # random strict submodel, OLS channels only (design-Gram clause (a))
    s.sub.c <- NA_real_
    s.sub.s <- NA_real_
    if(fam == "gaussian" && length(v.x.nms) > 1) {
        set.seed(seed + 500000L)
        v.sub <- sample(v.x.nms, sample(seq_len(length(v.x.nms) - 1), 1))
        fml.s <- reformulate(v.sub, response = y.nm)
        mo.so <- lm(fml.s, d)
        mo.ss <- lm(fml.s, d.s)
        m.so <- coef(summary(mo.so))
        m.ss <- coef(summary(mo.ss))
        s.sub.c <- dev.pct(coef(mo.ss), coef(mo.so))
        s.sub.s <- dev.pct(`names<-`(m.ss[, 2], rownames(m.ss)),
                           `names<-`(m.so[, 2], rownames(m.so)))
    }

    v.viol <- `if`(length(v.x.cnt) > 0, support.viol(d.s, d, v.x.cnt),
                   c(rate = 0, max = 0))
    s.e <- energy.dist(std.rows(d.s, d), std.rows(d, d))

    data.frame(n = nrow(d), arm = `if`(b.pinned, "pinned", "free"),
               coef.dev = s.coef, se.dev = s.se,
               sub.coef.dev = s.sub.c, sub.se.dev = s.sub.s,
               z.coef.dev = s.z,
               viol.rate = v.viol[["rate"]], viol.max = v.viol[["max"]],
               e.joint = s.e,
               full.rank = `if`(b.phi, b.fr, NA))
}

# +--------------------------------------------------------------------------+
# | Designs (the section-5 grid: scenarios 1, 3, 4, 6, 7, 8, 9, and PHI)     |
# +--------------------------------------------------------------------------+

load("results/d.sim.scenarios.RData")

l.designs <- list(
    list(nm = "sc1", i = 1L, fam = "gaussian",    lnk = "identity",
         sup = c(response = "real",     age = "real", income = "real")),
    list(nm = "sc3", i = 3L, fam = "gamma",       lnk = "log",
         sup = c(response = "positive", age = "real", income = "positive")),
    list(nm = "sc4", i = 4L, fam = "quasi-gamma", lnk = "log",
         sup = c(response = "zeroinfl", age = "real", income = "positive")),
    list(nm = "sc6", i = 6L, fam = "poisson",     lnk = "log",
         sup = c(response = "count",    age = "real", income = "positive")),
    list(nm = "sc7", i = 7L, fam = "gamma",       lnk = "log",
         sup = c(response = "positive", age = "real", income = "positive")),
    list(nm = "sc8", i = 8L, fam = "quasi-gamma", lnk = "log",
         sup = c(response = "zeroinfl", age = "real", income = "positive")),
    list(nm = "sc9", i = 9L, fam = "gaussian",    lnk = "identity",
         sup = NULL))

b.smoke <- nzchar(Sys.getenv("DGRAM_SMOKE"))
s.B <- `if`(b.smoke, 3L, 100L)
if(b.smoke) l.designs <- l.designs[1]
v.seeds <- 6489 + seq_len(s.B)

run.design <- function(nm, d, y.nm, fam, lnk, sup, b.phi = FALSE) {
    l.rows <- Map(function(seed) {
        l.arm <- Map(function(b.pinned) {
            d.row <- mc.one.dg(d, y.nm, fam, lnk, sup, seed, b.pinned,
                               b.phi)
            if(!is.null(d.row)) { d.row[["design"]] <- nm
                                  d.row[["seed"]]   <- seed }
            d.row
        }, c(FALSE, TRUE))
        do.call(rbind, Filter(Negate(is.null), l.arm))
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

if(!b.smoke && file.exists("phi_phase1_baseline.dta")) {
    cat("design phi (gaussian/identity) ...\n")
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
    l.out <- c(l.out, list(run.design("phi", d.app, "waist", "gaussian",
                                      "identity",
                                      c(waist = "real", age = "real"),
                                      b.phi = TRUE)))
} else if(!b.smoke) {
    cat("design phi skipped: phi_phase1_baseline.dta not in this tree\n")
}

d.out <- do.call(rbind, l.out)
d.out <- d.out[, c("design", "seed", "n", "arm", "coef.dev", "se.dev",
                   "sub.coef.dev", "sub.se.dev", "z.coef.dev",
                   "viol.rate", "viol.max", "e.joint", "full.rank")]
write.csv(d.out, "results/design-gram.csv", row.names = FALSE)

# +--------------------------------------------------------------------------+
# | Summary table and LaTeX body, both emitted here                         |
# +--------------------------------------------------------------------------+

agg <- function(v) {
    v <- v[!is.na(v)]
    `if`(length(v) == 0, "--", sprintf("%.4g (%.4g)", mean(v), sd(v)))
}
cat("\ndesign-gram, mean (sd) over", s.B, "replications, per arm\n")
cat(sprintf("%-5s %-7s %14s %14s %14s %12s %12s %10s\n", "des", "arm",
            "coef.dev%", "se.dev%", "z.coef.abs", "viol.rate%",
            "e.joint", "full.rank"))
for(nm in unique(d.out[["design"]])) {
    for(arm in c("free", "pinned")) {
        d.g <- d.out[d.out[["design"]] == nm & d.out[["arm"]] == arm, ]
        s.fr <- `if`(all(is.na(d.g[["full.rank"]])), "--",
                     sprintf("%d/%d", sum(d.g[["full.rank"]]), nrow(d.g)))
        cat(sprintf("%-5s %-7s %14s %14s %14s %12s %12s %10s\n", nm, arm,
                    agg(d.g[["coef.dev"]]), agg(d.g[["se.dev"]]),
                    agg(d.g[["z.coef.dev"]]), agg(d.g[["viol.rate"]]),
                    agg(d.g[["e.joint"]]), s.fr))
    }
}
cat("\nsubmodel spot-check (OLS channels), mean (sd)\n")
for(nm in unique(d.out[["design"]])) {
    d.g <- d.out[d.out[["design"]] == nm & d.out[["arm"]] == "pinned", ]
    if(all(is.na(d.g[["sub.coef.dev"]]))) next
    cat(sprintf("%-5s sub.coef.dev %s  sub.se.dev %s\n", nm,
                agg(d.g[["sub.coef.dev"]]), agg(d.g[["sub.se.dev"]])))
}

cat("\nLaTeX body (design-Gram table): design & free se & pinned se",
    "& viol & e.joint(pinned) \\\\\n")
for(nm in unique(d.out[["design"]])) {
    d.f <- d.out[d.out[["design"]] == nm & d.out[["arm"]] == "free", ]
    d.p <- d.out[d.out[["design"]] == nm & d.out[["arm"]] == "pinned", ]
    cat(sprintf("%s & %s & %s & %s & %s \\\\\n", nm,
                agg(d.f[["se.dev"]]), agg(d.p[["se.dev"]]),
                agg(d.p[["viol.rate"]]), agg(d.p[["e.joint"]])))
}
cat("\ndone.\n")
