#!/usr/bin/env r
# +------------------------------------------------------------+
# |  DP-Projection for zero-inflated scenarios (4, 8)          |
# +------------------------------------------------------------+
#
# Extension of DP-Projection to ZI quasi-gamma scenarios.
#
# Two arms, both on the ZI quasi-gamma scenarios 4 and 8:
#   Reference arm (gen.syn.dp.projected): the positive-only projection
#     fidelity engine.  It uses local (DFBETA) sensitivity, which is not
#     sound DP under Remedy B, so it is reported only as a non-private
#     fidelity reference (the deviation the projection achieves before
#     sound DP noise), and the zero indicator leaks through the copula.
#   Full DP-ZI arm (gen.syn.dp.projected.zi): a private zero-model
#     gamma_DP (mu_Z-GDP) plus the private positive coefficient beta_DP
#     (mu_P-GDP), both by subsample-and-aggregate, composed
#     sqrt(mu_Z^2 + mu_P^2)-GDP by mu-quadrature.  Z_syn is sampled from
#     gamma_DP, so the whole release is mu-GDP, equivalent to (eps, delta)
#     by the exact dual.
#
# +------------------------------------------------------------+
# |  Setup                                                     |
# +------------------------------------------------------------+
.sobol.main <- TRUE
source("lib/00-shared.r")
source("lib/gen-syn-family.r")
source("lib/07-calibrated-generator.r")
source("lib/09-symmetric-generator.r")
source("lib/11-projection-generator.r")
source("lib/28-dp-projection.r")
source("lib/dp-public-bounds.r")
rm(.sobol.main)

fit.positive.model <- function(d, y.nm) {
    v.pred <- intersect(c("age", "female", "race", "income",
                           "edu", "housing"), colnames(d))
    fml <- reformulate(v.pred, response = y.nm)
    d.pos <- d[d[[y.nm]] > 0, , drop = FALSE]
    s.np <- ncol(model.matrix(fml, d.pos))
    tryCatch(
        glm(fml, d.pos,
            family = quasi(link = "log", variance = "mu^2"),
            start = c(1, rep(0, s.np - 1)),
            control = glm.control(maxit = 200)),
        error = function(e) NULL)
}

# Predictor set shared by the positive and zero models (mirrors
# fit.positive.model's intersect.  Keep the two lists identical so the
# DP-ZI beta_DP and gamma_DP are comparable to the original refits).
zi.pred.names <- function(d)
    intersect(c("age", "female", "race", "income", "edu", "housing"),
              colnames(d))

# Original (non-private) zero-model gamma_hat: logistic I(y > 0) ~ X.
# The full DP-ZI arm compares the synthetic zero-model refit against this
# to measure whether the private zero indicator preserves the analyst's
# sign and significance decisions.
fit.zero.model <- function(d, y.nm) {
    v.pred <- zi.pred.names(d)
    d.z <- d[v.pred]
    d.z[[".Z"]] <- as.integer(d[[y.nm]] > 0)
    tryCatch(
        glm(reformulate(v.pred, response = ".Z"), d.z,
            family = binomial(link = "logit")),
        error = function(e) NULL)
}

load("results/d.sim.scenarios.RData")

v.test <- c(4, 8)
l.sup <- list(
    `4` = c(response = "zeroinfl", age = "real", income = "positive"),
    `8` = c(response = "zeroinfl", age = "real", income = "positive"))

s.B <- 20L
v.eps <- c(1, 2, 5, 10)
s.delta <- 1e-6

# +------------------------------------------------------------------------+
# |  Public boxes for the DP-ZI release                                    |
# +------------------------------------------------------------------------+
# The box width sets s_j = (u_j - l_j)/m, hence the noise scale, so the box must
# be public.  A generator that built its own box
# internally as [beta_hat +/- 4 SE] from the full-data fit and called it "public
# prior knowledge".  It was a statistic of the sensitive sample: a single
# adversarial row moved the noise scale by 5.7%.  This arm drives the
# zero-inflated DP table in the supplement, so that box reached the manuscript.
#
# The boxes below are +/- kappa * |beta_design| centred on zero, with beta_design
# the population coefficient under the published DGP (64-dp-design-constants.r) --
# a functional of the DGP's own parameters, computed from an independent draw at a
# public seed, and not the n=2000 seed-6489 analysis sample.  Centring on zero
# costs nothing (the width, hence the noise, is the same) and assumes strictly
# less: only a magnitude scale, not a sign or a location.
KAPPA   <- 3
l.const <- readRDS("results/dp-design-constants.rds")

zi.public.box <- function(sc, x.nms) {
    k <- l.const[[paste0("sc", sc)]]
    keep <- function(v) v[intersect(names(v), c("(Intercept)", x.nms,
                             "female1", "race2", "race3"))]
    v.b <- keep(k[["beta"]]); v.g <- keep(k[["gamma"]])
    list(lo   = -KAPPA * abs(v.b), hi   = KAPPA * abs(v.b),
         z.lo = -KAPPA * abs(v.g), z.hi = KAPPA * abs(v.g),
         disp = dp.dispersion.public(k[["sigma2"]],
                    source = paste("population dispersion under the published",
                                   "DGP; see 64-dp-design-constants.r")))
}

cat("DP-PROJECTION: ZERO-INFLATED SCENARIOS\n")
cat(strrep("=", 80), "\n")
cat(sprintf("B = %d reps, delta = %g\n", s.B, s.delta))
cat("Reference: positive-only projection (non-private fidelity engine,\n")
cat("           zero indicator leaks through the copula bootstrap).\n")
cat("Full DP-ZI: mu-GDP zero-model (mu_Z) + mu-GDP positive coefficient\n")
cat("           (mu_P), composed sqrt(mu_Z^2 + mu_P^2)-GDP, mu split 50/50.\n\n")

l.results <- list()
l.results.full <- list()
for(i in v.test) {
    s <- l.scenarios[[i]]; d <- s[["data"]]
    sup <- l.sup[[as.character(i)]]
    cat(sprintf("\n=== Scenario %d (quasi-gamma log, ZI) ===\n", i))
    cat(sprintf("n = %d, positive count = %d (%.1f%%)\n",
                nrow(d), sum(d[["response"]] > 0),
                100 * mean(d[["response"]] > 0)))

    mo.o <- fit.positive.model(d, "response")
    if(is.null(mo.o)) { cat("Original fit failed\n"); next }
    v.beta.o <- coef(mo.o)
    cat("beta_hat (positive only):", round(v.beta.o, 4), "\n")

    y.pos <- d[d[["response"]] > 0, "response"]

    cat(sprintf("\n%-6s  %8s  %8s  %8s  %6s  %10s\n",
                "eps", "max%.m", "max%.sd", "mean%", "n.ok", "sigma"))
    cat(strrep("-", 60), "\n")

    for(eps in v.eps) {
        v.seeds <- 6489L + seq_len(s.B)
        l.res <- Map(function(sd) {
            res <- tryCatch(
                gen.syn.dp.projected(d, y.nm = "response",
                                      family = "quasi-gamma",
                                      link = "log", supports = sup,
                                      eps = eps, delta = s.delta,
                                      seed = sd,
                                      clip.sigma = 3),
                error = function(e) NULL)
            if(is.null(res)) return(NULL)
            mo.s <- fit.positive.model(res[["d.syn"]], "response")
            if(is.null(mo.s)) return(NULL)
            v.b <- coef(mo.s)
            v.common <- intersect(names(v.beta.o), names(v.b))
            if(length(v.common) < 2) return(NULL)
            v.pct <- abs(100 * (v.b[v.common] - v.beta.o[v.common]) /
                          abs(v.beta.o[v.common]))
            list(maxpct = max(v.pct, na.rm = TRUE),
                 meanpct = mean(v.pct, na.rm = TRUE),
                 sigma = res[["sigma"]])
        }, v.seeds)
        l.res <- Filter(Negate(is.null), l.res)
        if(length(l.res) == 0) {
            cat(sprintf("%-6g  all failed\n", eps)); next
        }
        v.max <- mapply(function(r) r[["maxpct"]], l.res)
        v.mean <- mapply(function(r) r[["meanpct"]], l.res)
        v.sigma <- mapply(function(r) r[["sigma"]], l.res)
        cat(sprintf("%-6g  %7.1f%%  %7.1f%%  %6.1f%%  %3d/%d  %10.4f\n",
                    eps, mean(v.max), sd(v.max), mean(v.mean),
                    length(l.res), s.B, mean(v.sigma)))
        l.results[[length(l.results) + 1L]] <- data.frame(
            scenario = i, eps = eps, delta = s.delta,
            maxpct.mean = mean(v.max), maxpct.sd = sd(v.max),
            mean.mean = mean(v.mean), sigma = mean(v.sigma),
            n.success = length(l.res), n.total = s.B)
    }

    # ---- Full DP-ZI arm: private zero-model + private positive coef ----
    # Same positive metric as the reference arm (deviation of the synthetic
    # positive-only refit from the original beta_hat) so the arms are
    # comparable.  The full arm splits the mu-GDP budget, mu_P = mu_total /
    # sqrt(2) on the positive coefficient, so its positive deviation is
    # larger than the reference, and the gap is the cost of privatising the
    # zero-model at mu_Z.  The zero-model decision agreement (sign and
    # significance of the synthetic zero-model refit vs the original
    # gamma_hat) is the new quantity the release protects.
    x.nms <- zi.pred.names(d)
    l.box0 <- zi.public.box(i, x.nms)
    cat(sprintf("  subagg m: zero-model %d, positive %d (block-count rule, published rate)\n",
                dp.subagg.blocks(nrow(d), length(l.box0[["z.lo"]])),
                dp.subagg.blocks(nrow(d), length(l.box0[["lo"]]),
                                 rate = dp.public.zero.rate(i))))
    mo.zero.o <- fit.zero.model(d, "response")
    if(is.null(mo.zero.o)) {
        cat("Zero-model fit failed; skipping full arm\n"); next
    }
    v.gamma.o   <- coef(mo.zero.o)
    v.p.gamma.o <- summary(mo.zero.o)[["coefficients"]][, 4]

    cat("\n  full DP-ZI (mu split 50/50: mu_Z = mu_P = mu_total/sqrt2)\n")
    cat(sprintf("  %-6s  %8s  %8s  %7s  %7s  %6s\n",
                "eps", "pos.m%", "pos.sd%", "Zsign%", "Zsig%", "n.ok"))
    cat("  ", strrep("-", 56), "\n")

    for(eps in v.eps) {
        v.seeds <- 6489L + seq_len(s.B)
        l.rf <- Map(function(sd) {
            l.box <- zi.public.box(i, x.nms)
            res <- tryCatch(
                gen.syn.dp.projected.zi(d, y.nm = "response",
                    family = "quasi-gamma", link = "log", x.nms = x.nms,
                    eps = eps, delta = s.delta, seed = sd,
                    pos.rate = dp.public.zero.rate(i),
                    clip.lo   = l.box[["lo"]],   clip.hi   = l.box[["hi"]],
                    clip.z.lo = l.box[["z.lo"]], clip.z.hi = l.box[["z.hi"]],
                    dispersion = l.box[["disp"]]),
                error = function(e) NULL)
            if(is.null(res)) return(NULL)
            d.s <- res[["d.syn"]]
            # positive-component deviation vs original beta.o
            mo.ps <- fit.positive.model(d.s, "response")
            if(is.null(mo.ps)) return(NULL)
            v.b <- coef(mo.ps)
            v.cb <- intersect(names(v.beta.o), names(v.b))
            if(length(v.cb) < 2) return(NULL)
            v.pct <- abs(100 * (v.b[v.cb] - v.beta.o[v.cb]) /
                          abs(v.beta.o[v.cb]))
            # zero-model decision agreement vs original gamma.o
            d.zr <- d.s
            d.zr[[".Z"]] <- as.integer(d.s[["response"]] > 0)
            mo.zs <- tryCatch(
                glm(reformulate(x.nms, response = ".Z"), d.zr,
                    family = binomial(link = "logit")),
                error = function(e) NULL)
            if(is.null(mo.zs)) return(NULL)
            v.g   <- coef(mo.zs)
            v.p.g <- summary(mo.zs)[["coefficients"]][, 4]
            v.cg  <- intersect(names(v.gamma.o), names(v.g))
            list(maxpct  = max(v.pct, na.rm = TRUE),
                 meanpct = mean(v.pct, na.rm = TRUE),
                 signz = mean(sign(v.g[v.cg]) == sign(v.gamma.o[v.cg])),
                 sigz  = mean((v.p.g[v.cg] < 0.05) ==
                              (v.p.gamma.o[v.cg] < 0.05)),
                 sigma.beta  = res[["sigma.beta"]],
                 sigma.gamma = res[["sigma.gamma"]],
                 mu.Z = res[["mu.Z"]], mu.P = res[["mu.P"]],
                 mu.total = res[["mu.total"]])
        }, v.seeds)
        l.rf <- Filter(Negate(is.null), l.rf)
        if(length(l.rf) == 0) {
            cat(sprintf("  %-6g  all failed\n", eps)); next
        }
        v.max   <- mapply(function(r) r[["maxpct"]], l.rf)
        v.mean  <- mapply(function(r) r[["meanpct"]], l.rf)
        v.signz <- mapply(function(r) r[["signz"]], l.rf)
        v.sigz  <- mapply(function(r) r[["sigz"]], l.rf)
        cat(sprintf("  %-6g  %7.1f%%  %7.1f%%  %6.0f%%  %6.0f%%  %3d/%d\n",
                    eps, mean(v.max), sd(v.max),
                    100 * mean(v.signz), 100 * mean(v.sigz),
                    length(l.rf), s.B))
        l.results.full[[length(l.results.full) + 1L]] <- data.frame(
            scenario = i, eps = eps, delta = s.delta,
            mu.total = l.rf[[1]][["mu.total"]],
            mu.Z = l.rf[[1]][["mu.Z"]], mu.P = l.rf[[1]][["mu.P"]],
            maxpct.mean = mean(v.max), maxpct.sd = sd(v.max),
            mean.mean = mean(v.mean),
            zsign.mean = 100 * mean(v.signz),
            zsig.mean = 100 * mean(v.sigz),
            sigma.beta = mean(mapply(function(r) r[["sigma.beta"]], l.rf)),
            sigma.gamma = mean(mapply(function(r) r[["sigma.gamma"]], l.rf)),
            n.success = length(l.rf), n.total = s.B)
    }
}

d.zi.out <- do.call(rbind, l.results)
write.csv(d.zi.out, "results/mc-dp-projection-zi.csv", row.names = FALSE)
cat("\nResults written to results/mc-dp-projection-zi.csv\n")

if(length(l.results.full) > 0) {
    d.zi.full <- do.call(rbind, l.results.full)
    write.csv(d.zi.full, "results/mc-dp-projection-zi-full.csv",
              row.names = FALSE)
    cat("Full DP-ZI results written to",
        "results/mc-dp-projection-zi-full.csv\n")
}
