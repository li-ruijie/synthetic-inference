#!/usr/bin/env r
# +------------------------------------------------------------------------+
# |  66b-mc-dp-suffstat.r -- MC driver for Route 3, the whitened            |
# |  sufficient-statistic release                                           |
# +------------------------------------------------------------------------+
#
# Runs gen.syn.dp.suffstats (66-dp-suffstat.r) on the 62-mc-dp-gdp.r
# protocol so its rows sit like-for-like beside the subagg and Route 1
# numbers the paper reports.  Scenarios 1, 2, 9, and
# 10 (the bounded-response set, where dp.public.y.std exists), eps in
# {1, 2, 5, 10} at delta = 1e-6, B = 200 replications on scenario 1 and 25
# elsewhere, per-rep seeds 6489 + 1..B as 62- uses.
#
# Two dispersion modes per cell.
#
#   private   dispersion = NULL.  The y'y channel is released and sigma^2
#             comes from the released statistics, the default three-channel
#             quadrature (vech 0.45, xty 0.45, yy 0.10).
#   public    the population dispersion declared public, exactly Route 2's
#             free mode.  The y'y channel is dropped and its share
#             renormalised (vech 0.5, xty 0.5), which reproduces the survey's
#             survey's winning split, so this mode is the like-for-like
#             footing against the subagg rows.
#
# The estimand is the clipped-OLS one, as for Route 1.  The frame is
# clipped into the public support box (dp.public.clip.x) and the response
# into the public range (dp.public.y.std's y.clip) before the original
# fit, since that is the target the mechanism actually releases.  On the
# realised samples both clips are the identity map (the 62- header records
# the per-row binding probability near 1e-9).
#
# Metrics per cell, aggregated over reps: maxpct (max relative coefficient
# deviation over the substantive filter |beta_hat| >= 0.01, the 62-
# filter), sign, significance, and joint decision agreement at the 5%
# level, and calibrated 95% CI coverage of the original beta_hat.
#
# The calibrated interval, recentred.  The corrected
# standard-error equation folds the DP noise variance into the refit SE as
# SE_corr^2 = SE_syn^2 + sigma_beta^2 with no column scale.  Route 3's
# release is a ridge solve, so its noise is neither isotropic nor centred on the
# OLS target, and the first live run of this driver read coverage
# as low as 66% at the strong-shrinkage cells (lambda = n) for the interval
# centred at beta.dp with the then-returned first-order fold lambda M beta_w
# folded into the SE: the fold reads half the true shrinkage at lambda = n, and
# it shrinks with the noisy beta_w exactly when the error is largest.  The
# calibrated interval is centred at the lambda-0 solve instead,
#
#     beta.infer_j +/- 1.96 sqrt(SE_syn,j^2 + Cov.beta.infer[j, j])
#
# the mechanism's returned inference pair (66-dp-suffstat.r), an exact
# post-processing de-shrink of the same released statistics whose
# covariance plug-in is evaluated at the de-shrunk coefficient.  beta.dp
# stays the point estimate behind maxpct and the decision columns.  The
# coverage column reads the recentred interval, and it is the first real
# test of the returned covariance, since Task 11 verifies the mechanism's
# accounting and leaves the analyst-side error law untested.
#
# The refit on d.syn recovers beta.dp exactly (the projection property,
# asserted to 1e-8 by verify-gdp-mechanisms.R Task 11), so maxpct measured
# through the refit equals maxpct on the release itself.  The refit exists
# to give the analyst-side SE_syn.  d.syn is emitted under x.public = TRUE
# (the premise Route 1 and Route 2 share, so the comparison is
# like-for-like), while the release's own calibration does not read X, which
# is what the whitened-scale sensitivities are for.
#
# The full run is gated: run the smoke pass first (SUFFSTAT_SMOKE=1) to
# confirm the driver executes, then the full pass.
#
# Output: mc-dp-suffstat.csv.
# +------------------------------------------------------------------------+

# +------------------------------------------------------------------------+
# |  Setup                                                                 |
# +------------------------------------------------------------------------+
source("lib/00-shared.r")
.shared.loaded    <- TRUE               # the 66- chain re-sources 00-shared.r, and
                                        # an unguarded re-init clears globalenv
.sobol.main       <- TRUE               # suppress the 07-/09-/11- demo blocks
.dp.suppress.demo <- TRUE               # suppress the 28- demo block
source("lib/dp-public-bounds.r")            # the ONLY permitted source of a DP bound
source("lib/66-dp-suffstat.r")              # sources 28-dp-projection.r (suppressed)
rm(.shared.loaded, .sobol.main, .dp.suppress.demo)

load("results/d.sim.scenarios.RData")
l.const <- readRDS("results/dp-design-constants.rds")

# +------------------------------------------------------------------------+
# |  Configuration                                                         |
# +------------------------------------------------------------------------+
# Smoke trims the grid to prove the driver runs without the full cost.
SMOKE <- nzchar(Sys.getenv("SUFFSTAT_SMOKE"))

v.eps   <- c(1, 2, 5, 10)               # target epsilon points (tables' grid)
s.delta <- 1e-6                         # fixed delta for the (eps, delta) dual
v.scen  <- c(1L, 2L, 9L, 10L)           # bounded-response set (public B_y)
V.PREDS <- c("age", "female", "race", "income")
SUBSTANTIVE <- 0.01                     # maxpct is taken over |beta_hat| >= this

# reps: dense for scenario 1 (headline), lighter for the broader set.
b.for <- function(i) if(i == 1L) 200L else 25L

if(SMOKE) {
    v.eps <- 10; v.scen <- 1L
    b.for <- function(i) 3L
    cat("[SMOKE] trimmed grid: eps = 10, scenario 1, both modes, B = 3\n")
}

# +------------------------------------------------------------------------+
# |  Helpers                                                               |
# +------------------------------------------------------------------------+
# Factor levels re-asserted defensively, the 65-dp-suffstat-survey.r
# pattern, so model.matrix emits the canonical six-column design the
# whitening map is built on.
prep <- function(d) {
    d <- `[[<-`(d, "female",
                value = factor(as.character(d[["female"]]), levels = c("0", "1")))
    `[[<-`(d, "race",
           value = factor(as.character(d[["race"]]), levels = c("1", "2", "3")))
}

# The clipped-OLS estimand frame: X into the public support box, y into the
# public range.  Identity on the realised samples, and the universe the
# sensitivities are computed over.
clip.frame <- function(d, l.y) {
    d <- dp.public.clip.x(d)
    `[[<-`(d, "response",
           value = pmin(pmax(d[["response"]], l.y[["y.clip"]][[1L]]),
                        l.y[["y.clip"]][[2L]]))
}

# Gaussian identity only on this route (the 66- scope).  NULL on failure so
# a degenerate refit drops the rep, not the run.
fit.ols <- function(d) tryCatch(lm(reformulate(V.PREDS, "response"), d),
                                error = function(e) NULL)

# Utility and decision-agreement metrics vs the original full-data fit, the
# 62- definitions with the recentred interval described in the header.
# v.b.inf and v.var.inf are named per-coefficient vectors on the native
# beta scale, beta.infer and diag(Cov.beta.infer) from the mechanism.
suffstat.metrics <- function(mo.real, mo.syn, v.b.inf, v.var.inf) {
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

    v.se.cor <- sqrt(v.se.syn^2 + v.var.inf[v.nms])
    b.cov    <- (v.b.real >= v.b.inf[v.nms] - 1.96 * v.se.cor) &
                (v.b.real <= v.b.inf[v.nms] + 1.96 * v.se.cor)

    list(maxpct = max(v.pct, na.rm = TRUE),
         sign = mean(b.sign, na.rm = TRUE),
         sig  = mean(b.sig,  na.rm = TRUE),
         both = mean(b.sign & b.sig, na.rm = TRUE),
         cov  = mean(b.cov, na.rm = TRUE))
}

# One (scenario, mode, eps) cell: B reps, aggregate to a one-row frame.
# lambda is constant per cell under the public curve (fixed public seed),
# lambda.eff varies by rep through the PD repair, so both are reported as
# cell means.
run.cell <- function(i, s.disp, eps, gen.rep) {
    s.B <- b.for(i)
    l.m <- Filter(Negate(is.null), Map(gen.rep, 6489L + seq_len(s.B)))
    if(length(l.m) == 0L) {
        cat(sprintf("  Sc %2d %-8s eps %-4g : no successful reps\n",
                    i, s.disp, eps))
        return(NULL)
    }
    grab <- function(k) vapply(l.m, function(r) r[[k]], numeric(1))
    v.mp <- grab("maxpct"); s.mu <- mean(grab("mu"))
    cat(sprintf(paste0("  Sc %2d %-8s eps %-4g : mu %.3f  lam.eff %8.1f  ",
                       "max%% %7.1f  sign %3.0f  sig %3.0f  both %3.0f  ",
                       "cov %3.0f  (%d/%d)\n"),
                i, s.disp, eps, s.mu, mean(grab("lambda.eff")), mean(v.mp),
                100 * mean(grab("sign")), 100 * mean(grab("sig")),
                100 * mean(grab("both")), 100 * mean(grab("cov")),
                length(l.m), s.B))
    data.frame(route = "suffstat", scenario = i, disp = s.disp, mu = s.mu,
               eps = eps, delta = s.delta,
               maxpct.mean = mean(v.mp), maxpct.sd = sd(v.mp),
               sign = mean(grab("sign")), sig = mean(grab("sig")),
               both = mean(grab("both")), cov = mean(grab("cov")),
               lambda = mean(grab("lambda")),
               lambda.eff = mean(grab("lambda.eff")), B = length(l.m))
}

# +------------------------------------------------------------------------+
# |  Run                                                                   |
# +------------------------------------------------------------------------+
cat("Route 3 MC -- the whitened sufficient-statistic release\n")
cat(strrep("=", 92), "\n")
cat(sprintf(paste0("delta = %g, mu.split (private) = (vech 0.45, xty 0.45, ",
                   "yy 0.10); the public mode\n"), s.delta))
cat("renormalises to (0.5, 0.5, 0).  Bounds are functionals of the published DGP\n")
cat("(64-dp-design-constants.r); the estimand is the clipped-OLS one.\n\n")

l.res <- Map(function(i) {
    l.w <- dp.public.whitening(i)
    l.y <- dp.public.y.std(i)
    l.s <- dp.public.suffstat.sens(l.w, l.y)
    cat(sprintf("--- Scenario %d ---\n", i))
    dp.bounds.announce(B.x.white = l.w, y.std = l.y, suffstat = l.s)

    d <- clip.frame(prep(l.scenarios[[i]][["data"]])[c(V.PREDS, "response")], l.y)
    mo.o <- fit.ols(d)
    if(is.null(mo.o)) { cat(sprintf("  Sc %d original fit failed\n", i)); return(NULL) }

    s.s2 <- l.const[[paste0("sc", i)]][["sigma2"]]
    l.disp <- list(
        private = NULL,
        public  = dp.dispersion.public(s.s2, source = paste0(
            "population dispersion under the published DGP ",
            "(64-dp-design-constants.r)")))

    Map(function(s.mode) {
        disp <- l.disp[[s.mode]]
        Map(function(eps) {
            run.cell(i, s.mode, eps, function(sd) {
                res <- tryCatch(gen.syn.dp.suffstats(d, "response", sc = i,
                            eps = eps, delta = s.delta, whitening = l.w,
                            y.std = l.y, sens = l.s, dispersion = disp,
                            x.public = TRUE, seed = sd),
                        error = function(e) NULL)
                if(is.null(res)) return(NULL)
                mo.s <- fit.ols(res[["d.syn"]])
                if(is.null(mo.s)) return(NULL)
                v.b.inf   <- res[["beta.infer"]]
                v.var.inf <- `names<-`(diag(res[["Cov.beta.infer"]]),
                                       names(res[["beta.dp"]]))
                m <- suffstat.metrics(mo.o, mo.s, v.b.inf, v.var.inf)
                if(is.null(m)) return(NULL)
                m[["mu"]]         <- res[["mu.total"]]
                m[["lambda"]]     <- res[["lambda"]]
                m[["lambda.eff"]] <- res[["lambda.eff"]]
                m
            })
        }, v.eps)
    }, names(l.disp))
}, v.scen)

# Flatten the nested Maps down to the data frames (the 62- helper: a plain
# unlist() would tear the data frames apart).
flatten.df <- function(x) {
    if(is.null(x)) return(list())
    if(is.data.frame(x)) return(list(x))
    unlist(Map(flatten.df, x), recursive = FALSE)
}
d.out <- do.call(rbind, Filter(Negate(is.null), flatten.df(l.res)))

# +------------------------------------------------------------------------+
# |  Write                                                                 |
# +------------------------------------------------------------------------+
# Sanity: the reported mu maps back to (eps, delta) through the dual.
v.delta.check <- mapply(gdp.to.delta, d.out[["mu"]], d.out[["eps"]])
cat(sprintf("\nDual check: max |delta(mu, eps) - %g| = %.2e\n",
            s.delta, max(abs(v.delta.check - s.delta))))

s.out <- `if`(SMOKE, "results/mc-dp-suffstat-smoke.csv",
              "results/mc-dp-suffstat.csv")
write.csv(d.out, s.out, row.names = FALSE)
cat(sprintf("Results written to %s (%d rows)\n",
            sub("^\\.\\./", "", s.out), nrow(d.out)))

cat("\nComparators, quoted from the recorded values (not recomputed here):\n")
cat("  subagg kappa = 3, sc 1, eps = 10     21.5%  (s25, B = 200, m = 66 rule)\n")
cat("  deprecated local path, sc 1, eps = 10  12.3%  (s14, B = 25, no guarantee)\n")
cat("Copy every reported number FROM THE CSV, with this script named.\n")
