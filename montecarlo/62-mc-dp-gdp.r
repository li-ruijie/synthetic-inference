#!/usr/bin/env r
# +------------------------------------------------------------------------+
# |  62-mc-dp-gdp.r -- mu-Gaussian-DP re-run driver for both sound routes   |
# +------------------------------------------------------------------------+
#
# Remedy B, Task 6.  Re-runs the DP-projection utility study under the sound
# mu-GDP mechanisms (Dong, Roth and Su 2022), replacing the earlier
# local-sensitivity calibration.  Two routes:
#
#   Route 2 (primary, assumption-free): subsample-and-aggregate for
#     canonical-link GLMs, gen.syn.dp.projected.subagg (log-link scenarios
#     3, 6 and Gaussian-identity scenarios 1, 2, 9, 10).  Per-coordinate global
#     sensitivity s_j = (u_j - ell_j) / m from public coefficient bounds, so
#     the scaled aggregate is sqrt(p)/sigma.Q-GDP.  Fidelity is approximate
#     (block-aggregate estimand).
#   Route 1 (attained special case): public-covariate exact release for Gaussian
#     OLS, gen.syn.dp.ols.public.  Closed-form global sensitivity
#     Delta_beta = 2 B_y max_i ||(X'X)^{-1} x_i||, Gaussianised variance
#     release, exact beta_hat_DP recovery.
#
# Both report mu.total.  The equivalent (eps, delta) is recovered by the dual
# gdp.to.delta at delta = 1e-6.  Metrics per (route, scenario, eps): coefficient
# deviation (max% over substantive coefficients), decision agreement (sign,
# significance, both), and calibrated 95% CI coverage of the original beta_hat
# with the DP noise variance folded into the SE.
#
# Public bounds.  Every bound below is a functional
# of the published DGP, computed by 64-dp-design-constants.r from an independent
# draw at a public seed.  None is a statistic of the analysis sample.
#
#   Route 2 box   +/- kappa * |beta_design|, centred on zero, swept over
#                 kappa in {1.5, 3, 5} and reported as a function of kappa.
#   Route 2 disp  the population dispersion, declared public (no budget spent).
#   Route 1 B_y   the DGP's own response bound, with an explicit public clipping
#                 step where the response is unbounded (scenarios 1, 2, 10).
#
# What these replaced, and why it mattered.  The box was [beta_hat +/- 4 SE] from
# the full-data fit, described in this header as "public prior knowledge", and
# B.y was 1.5 * max|y|, the sample maximum.  Both are statistics of the sensitive
# data, and both set the noise scale, so the noise scale was a function of the
# data being protected and no DP theorem could be instantiated.  This driver backs
# the four differentially private release tables of the supplement, the baseline,
# the variants, the multi-link and the scaling tables, so those bounds are what
# the submission's DP numbers rested on.
#
# Expect the numbers to move, by more than a factor of ten on the coefficient
# noise.  A public box cannot see the standard errors, so it must be far wider
# than 4 SE.  The relative coefficient noise is 2 kappa sqrt(q) / (m mu), about
# 40% at kappa = 3, eps = 10, against approximately 2-3% before.  Decision
# agreement and CI coverage survive.  The scale-dependent max% does not, and it
# is dominated by `income` (beta ~ 0.013 on a predictor of order 1e3).  Lead the
# refreshed tables with agreement and coverage.
#
# The full run is gated: run the smoke pass first (GDP_SMOKE=1) to confirm the
# driver executes, then the full pass.
#
# +------------------------------------------------------------------------+
# |  Setup                                                                 |
# +------------------------------------------------------------------------+
.sobol.main <- TRUE
source("lib/00-shared.r")
source("lib/gen-syn-family.r")
source("lib/07-calibrated-generator.r")
source("lib/09-symmetric-generator.r")
source("lib/11-projection-generator.r")
source("lib/28-dp-projection.r")
.dp.public.no.main <- TRUE
source("lib/33-dp-projection-public-x.r")           # gen.syn.dp.ols.public (Route 1)
rm(.sobol.main, .dp.public.no.main)

load("results/d.sim.scenarios.RData")

# +------------------------------------------------------------------------+
# |  Configuration                                                         |
# +------------------------------------------------------------------------+
# Smoke trims the grid to prove the driver runs without the full re-run cost.
SMOKE <- nzchar(Sys.getenv("GDP_SMOKE"))

v.eps   <- c(1, 2, 5, 10)               # target epsilon points (tables' grid)
s.delta <- 1e-6                         # fixed delta for the (eps, delta) dual
s.m     <- NULL                         # Route 2 block count, NULL applies floor(n / (5 p))

v.route1 <- c(1, 2, 9, 10)              # Gaussian identity scenarios (OLS exact release)
v.route2 <- c(1, 2, 3, 6, 9, 10)        # subsample-aggregate on every covered scenario

v.kappa <- c(1.5, 3, 5)                 # PUBLIC box half-width, in |beta_design|
SUBSTANTIVE <- 0.01                     # maxpct is taken over |beta_hat| >= this

# The public design constants: population values under the published DGP,
# computed from an independent draw at a public seed (64-dp-design-constants.r).
l.const <- readRDS("results/dp-design-constants.rds")

# reps: dense for scenario 1 (headline), lighter for the broader set.
b.for <- function(i) if(i == 1L) 200L else 25L

if(SMOKE) {
    v.eps <- 10; v.route1 <- 1; v.route2 <- c(1, 3); v.kappa <- 3
    b.for <- function(i) 2L
    cat("[SMOKE] trimmed grid: eps = 10, kappa = 3, Route 1 = {1}, Route 2 = {1,3}, B = 2\n")
}

# +------------------------------------------------------------------------+
# |  Helpers                                                               |
# +------------------------------------------------------------------------+
# Public predictor sets (covariate metadata is public domain knowledge).
preds.ols     <- function(d) intersect(c("age", "female", "race", "income"),
                                       colnames(d))
preds.loglink <- function(d) intersect(c("age", "female", "race", "income",
                                         "edu", "housing"), colnames(d))

# Fit the analyst's model.  lm for gaussian/identity, glm otherwise.  Returns
# NULL on failure so a rank-deficient synthetic refit drops the rep, not the run.
fit.model <- function(d, y.nm, x.nms, family, link) {
    fml <- reformulate(x.nms, response = y.nm)
    tryCatch(switch(family,
        gaussian      = lm(fml, d),
        gamma         = glm(fml, d, family = Gamma(link = link),
                            control = glm.control(maxit = 200)),
        poisson       = glm(fml, d, family = poisson(link = link),
                            control = glm.control(maxit = 200)),
        `quasi-gamma` = {
            s.np <- ncol(model.matrix(fml, d))
            glm(fml, d, family = quasi(link = "log", variance = "mu^2"),
                start = c(1, rep(0, s.np - 1)),
                control = glm.control(maxit = 200))
        }),
        error = function(e) NULL)
}

# +------------------------------------------------------------------------+
# |  The public clip box                                                   |
# +------------------------------------------------------------------------+
# A clip box must not be built as:
#
#     clip.box <- function(mo) {
#         v.b  <- coef(mo)                                  # full-data fit
#         v.se <- summary(mo)[["coefficients"]][, 2]        # full-data SEs
#         list(lo = v.b - clip.K * v.se, hi = v.b + clip.K * v.se)
#     }
#
# and the file's own header called that box "public prior knowledge".  It is not.
# It is a statistic of the sensitive data, and the box width sets s_j, hence the
# noise scale.  The sensitivity argument binds one box and applies it to both
# neighbours, so a data-derived box cannot instantiate the
# theorem at all: there is no guarantee to state, not merely a weaker one.  This
# driver backs the four differentially private release tables of the supplement,
# so that box is what the submission's DP numbers rest on.
#
# The box is now +/- kappa * |beta_design| centred on zero, with beta_design the
# Population coefficient under the published DGP (64-dp-design-constants.r).
# Centring on zero costs nothing (same width, same noise) and assumes strictly
# less: a magnitude scale only, no sign and no location.
#
# The cost is real, and it is large.  The relative coefficient noise is
# 2 kappa sqrt(q) / (m mu), the same for every coordinate, which at kappa = 3,
# m = 20, eps = 10 was about 40%, against approximately 2-3% for an SE-derived box.
# (The block-count rule floor(n / (5 p)) sets m, which
# shrinks the noise in proportion to 20/m.)  A
# public box cannot peek at the standard errors, so it must be far wider.  The
# kappa sweep reports utility as A function of kappa rather than hiding the
# dependence, which is the correct answer to the obvious referee question of
# whether the box smuggles in the answer.
clip.box <- function(i, kappa, v.names) {
    v.b <- l.const[[paste0("sc", i)]][["beta"]]
    v.b <- v.b[intersect(names(v.b), v.names)]
    list(lo = -kappa * abs(v.b), hi = kappa * abs(v.b))
}

# +------------------------------------------------------------------------+
# |  Route 1's public response bound, and its clipping step                |
# +------------------------------------------------------------------------+
# A bound such as B.Y.MULT * max(abs(response)) is the sample maximum, which is
# unsound, since Delta_beta = 2 B_y max_i ||(X'X)^-1 x_i||
# is defined through B_y, so the noise scale was a function of the private data.
#
# Worse, it was circular.  Route 1 is evaluated on scenarios 1, 2, 9 and 10.
# Scenarios 1, 2 and 10 draw y from N(50, 15), which is unbounded, so no finite
# B_y exists for them: 1.5 * max|y| made the bounded-universe assumption true for
# the observed sample by construction.  OLS beta_hat has unbounded global
# sensitivity when y is unbounded, so no public B_y means no Route 1.
#
# Clip y to a public range as an explicit
# pre-processing step, taken from the DGP's own marginal (mean +/- 6 sd), not from
# the sample.  The clipping is a deterministic per-row map with public parameters,
# so it costs no budget, and it makes the bounded-response premise TRUE rather
# than assumed.  The estimand becomes the clipped-OLS one and must be declared as
# such in the public-covariate exact release proposition.  At 6 sd it binds with
# probability ~1e-9 per row, so on any realised sample it is the identity map.
#
# Scenario 9 needs no clipping: its DGP is round(pmax(pmin(y, 100), 0)), so
# y in [0, 100] structurally and B_y = 100 is exact, public, and spends no budget.
route1.bounds <- function(i) {
    k <- l.const[[paste0("sc", i)]]
    list(B.y = k[["B.y"]], y.clip = k[["y.clip"]], source = k[["y.source"]])
}

clip.response <- function(d, l.b) {
    if(is.null(l.b[["y.clip"]])) return(d)          # structurally bounded already
    `[[<-`(d, "response", value = pmin(pmax(d[["response"]], l.b[["y.clip"]][[1L]]),
                                       l.b[["y.clip"]][[2L]]))
}

# Utility and decision-agreement metrics vs the original full-data fit.
# v.sigma.dp: named per-coefficient DP noise SD on the native beta scale, used
# to calibrate the coverage SE (se.corr = sqrt(se.syn^2 + sigma.dp^2)).  maxpct
# is taken over substantive coefficients (|beta_o| >= substantive) so a near-zero
# coefficient does not make the relative deviation diverge.
#
# Caveat (reporting, revisit at table time).  maxpct is scale-dependent.  In the
# Gaussian scenarios the income coefficient is small on its native scale
# (beta ~ 0.013 for a predictor of order 1e3) yet sits just above substantive,
# so under the sound worst-case noise (much larger than the earlier
# local calibration) its relative deviation explodes and dominates maxpct.  For
# those scenarios decision agreement (sign, sig, both, cov) is the meaningful
# summary.  When refreshing the tables, either raise substantive, drop income
# from maxpct, or lead the table with agreement rather than maxpct.
gdp.metrics <- function(mo.real, mo.syn, v.sigma.dp) {
    if(is.null(mo.syn)) return(NULL)
    s.real <- summary(mo.real)[["coefficients"]]
    s.syn  <- summary(mo.syn)[["coefficients"]]
    v.nms  <- intersect(rownames(s.real), rownames(s.syn))
    if(length(v.nms) < 2) return(NULL)

    v.b.real <- s.real[v.nms, 1]; v.p.real <- s.real[v.nms, 4]
    v.b.syn  <- s.syn[v.nms, 1];  v.p.syn  <- s.syn[v.nms, 4]
    v.se.syn <- s.syn[v.nms, 2]

    v.sub <- v.nms[abs(v.b.real) >= SUBSTANTIVE]
    v.pct <- if(length(v.sub) > 0)
        abs(100 * (v.b.syn[v.sub] - v.b.real[v.sub]) / abs(v.b.real[v.sub]))
        else NA_real_

    b.sign <- sign(v.b.syn) == sign(v.b.real)
    b.sig  <- (v.p.syn < 0.05) == (v.p.real < 0.05)

    v.sd     <- v.sigma.dp[v.nms]
    v.se.cor <- sqrt(v.se.syn^2 + v.sd^2)
    b.cov    <- (v.b.real >= v.b.syn - 1.96 * v.se.cor) &
                (v.b.real <= v.b.syn + 1.96 * v.se.cor)

    list(maxpct = max(v.pct, na.rm = TRUE),
         sign = mean(b.sign, na.rm = TRUE),
         sig  = mean(b.sig,  na.rm = TRUE),
         both = mean(b.sign & b.sig, na.rm = TRUE),
         cov  = mean(b.cov, na.rm = TRUE))
}

# One (route, scenario, kappa, eps) cell: B reps, aggregate to a one-row frame.
run.cell <- function(route, i, kappa, eps, gen.rep) {
    s.B <- b.for(i)
    l.m <- Filter(Negate(is.null), Map(gen.rep, 6489L + seq_len(s.B)))
    if(length(l.m) == 0L) {
        cat(sprintf("  %-7s Sc %2d k %-4g eps %-4g : no successful reps\n",
                    route, i, kappa, eps))
        return(NULL)
    }
    grab <- function(k) vapply(l.m, function(r) r[[k]], numeric(1))
    v.mp <- grab("maxpct"); s.mu <- mean(grab("mu"))
    cat(sprintf("  %-7s Sc %2d k %-4g eps %-4g : mu %.3f  max%% %7.1f  sign %3.0f  sig %3.0f  both %3.0f  cov %3.0f  (%d/%d)\n",
                route, i, kappa, eps, s.mu, mean(v.mp), 100 * mean(grab("sign")),
                100 * mean(grab("sig")), 100 * mean(grab("both")),
                100 * mean(grab("cov")), length(l.m), s.B))
    data.frame(route = route, scenario = i, kappa = kappa, mu = s.mu, eps = eps,
               delta = s.delta, maxpct.mean = mean(v.mp), maxpct.sd = sd(v.mp),
               sign = mean(grab("sign")), sig = mean(grab("sig")),
               both = mean(grab("both")), cov = mean(grab("cov")), B = length(l.m))
}

# +------------------------------------------------------------------------+
# |  Run                                                                   |
# +------------------------------------------------------------------------+
cat("mu-GDP DP-PROJECTION RE-RUN (Route 2 primary, Route 1 restrictive special case)\n")
cat(strrep("=", 92), "\n")
cat(sprintf("delta = %g, blocks m = the block-count rule floor(n / (5 p)), PUBLIC box +/- kappa*|beta_design| in %s\n",
            s.delta, paste(v.kappa, collapse = ", ")))
cat("Boxes and response bounds are POPULATION values under the published DGP\n")
cat("(64-dp-design-constants.r), not statistics of the analysis sample.\n\n")

# ---- Route 2: subsample-and-aggregate (primary) --------------------------
# Route 2 assumes public covariates.  Its dispersion is taken as public too: on
# these designs it is a published constant of the DGP, and in the field it is the
# R^2 or residual SD a paper reports alongside its coefficient table.  That mode
# spends no budget, so beta_DP is exactly as noisy as before, and it makes the
# release a post-processing of (beta_DP, sigma^2_public).
l.r2 <- Map(function(i) {
    s <- l.scenarios[[i]]; d <- s[["data"]]
    fam <- s[["family"]]; lnk <- s[["link"]]
    x.nms <- preds.loglink(d)
    mo.o <- fit.model(d, "response", x.nms, fam, lnk)
    if(is.null(mo.o)) { cat(sprintf("  Sc %d original fit failed\n", i)); return(NULL) }
    cat(sprintf("  Sc %d subagg m = %d (block-count rule)\n", i,
                dp.subagg.blocks(nrow(d), length(coef(mo.o)))))
    s.s2 <- l.const[[paste0("sc", i)]][["sigma2"]]
    disp <- `if`(is.na(s.s2), NULL, dp.dispersion.public(s.s2,
        source = "population dispersion under the published DGP (64-dp-design-constants.r)"))

    Map(function(kappa) {
        l.clip <- clip.box(i, kappa, names(coef(mo.o)))
        Map(function(eps) {
            run.cell("subagg", i, kappa, eps, function(sd) {
                res <- tryCatch(gen.syn.dp.projected.subagg(d, y.nm = "response",
                            family = fam, link = lnk, x.nms = x.nms, eps = eps,
                            delta = s.delta, m = s.m, clip.lo = l.clip[["lo"]],
                            clip.hi = l.clip[["hi"]], dispersion = disp, seed = sd),
                        error = function(e) NULL)
                if(is.null(res)) return(NULL)
                mo.s <- fit.model(res[["d.syn"]], "response", x.nms, fam, lnk)
                m <- gdp.metrics(mo.o, mo.s, res[["sigma"]])
                if(is.null(m)) return(NULL)
                m[["mu"]] <- res[["mu.total"]]; m
            })
        }, v.eps)
    }, v.kappa)
}, v.route2)

# ---- Route 1: public-covariate exact release (attained) -----------------
# The response is clipped to a public range first (see route1.bounds).  That step
# makes B_y a bound rather than an assumption, and it moves the estimand to the
# clipped-OLS one, so mo.o is fitted on the clipped data: the target the analyst
# is trying to recover is what the mechanism actually releases.  On these
# designs the clip binds with probability ~1e-9 per row, so the two estimands
# coincide on any realised sample.
l.r1 <- Map(function(i) {
    s <- l.scenarios[[i]]
    l.b  <- route1.bounds(i)
    d    <- clip.response(s[["data"]], l.b)
    x.nms <- preds.ols(d)
    mo.o <- fit.model(d, "response", x.nms, "gaussian", "identity")
    if(is.null(mo.o)) { cat(sprintf("  Sc %d original fit failed\n", i)); return(NULL) }

    Map(function(eps) {
        run.cell("public", i, NA_real_, eps, function(sd) {
            res <- tryCatch(gen.syn.dp.ols.public(d, y.nm = "response",
                        x.nms = x.nms, eps = eps, delta = s.delta,
                        B.y = l.b[["B.y"]], seed = sd),
                    error = function(e) NULL)
            if(is.null(res)) return(NULL)
            mo.s <- fit.model(res[["d.syn"]], "response", x.nms,
                              "gaussian", "identity")
            if(is.null(mo.s)) return(NULL)
            v.sd <- `names<-`(rep(res[["sigma.beta"]], length(coef(mo.s))),
                              names(coef(mo.s)))
            m <- gdp.metrics(mo.o, mo.s, v.sd)
            if(is.null(m)) return(NULL)
            m[["mu"]] <- res[["mu.total"]]; m
        })
    }, v.eps)
}, v.route1)

# Flatten the nested Maps down to the data frames.  A plain unlist() would tear
# the data frames apart (a data.frame is a list), so recurse until one is found.
flatten.df <- function(x) {
    if(is.null(x)) return(list())
    if(is.data.frame(x)) return(list(x))
    unlist(Map(flatten.df, x), recursive = FALSE)
}
l.results <- c(flatten.df(l.r2), flatten.df(l.r1))

# +------------------------------------------------------------------------+
# |  Write                                                                 |
# +------------------------------------------------------------------------+
d.out <- do.call(rbind, Filter(Negate(is.null), l.results))

# Sanity: the reported mu maps back to (eps, delta) through the dual.
d.out[["delta.check"]] <- mapply(gdp.to.delta, d.out[["mu"]], d.out[["eps"]])
s.max.dev <- max(abs(d.out[["delta.check"]] - s.delta))
cat(sprintf("\nDual check: max |delta(mu, eps) - %g| = %.2e\n", s.delta, s.max.dev))
d.out[["delta.check"]] <- NULL

s.out <- if(SMOKE) "results/mc-dp-gdp-smoke.csv" else "results/mc-dp-gdp.csv"
write.csv(d.out, s.out, row.names = FALSE)
cat(sprintf("Results written to %s (%d rows)\n", sub("^\\.\\./", "", s.out), nrow(d.out)))
