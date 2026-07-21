#!/usr/bin/env r
# +------------------------------------------------------------------------+
# |  63-dp-audit-mugdp.r -- empirical DP audit of the mu-GDP routes         |
# +------------------------------------------------------------------------+
#
# Why this script exists.  The paper's supplement S6 audit claim had no
# script behind it.  (S6 states that empirical audits are not a substitute
# for the analytic guarantee but detect implementation errors the analytic
# composition cannot see, which is what this script supplies.)
# The existing audits (50-, 51-, 51b-) cover the deprecated local-sensitivity
# variants V1-V5 rather than the two mu-GDP routes, so this script answers:
#
#   whether the synthetic release of Route 2 delivers the mu it claims.
#
# The attack.  DP is a statement about neighbouring datasets, so the audit is a
# distinguishing game on a neighbouring pair (Jagielski et al. 2020; Nasr et al.
# 2021).  Fix n-1 background rows.  The two worlds differ in one row, a probe row:
#
#   world 0   probe row sits on the regression line   (residual 0)
#   world 1   probe row sits PROBE.R residual SDs off it
#
# The adversary sees only the synthetic release and must say which world.  Its
# statistic is the residual scale of the release, sigma(lm(y ~ x, d.syn)),
# because Route 2 reproduces the real residual scale exactly: the probe row inflates
# RSS by (PROBE.R * sigma)^2, and the release copies it straight through.
#
# The bound.  For a mu-GDP mechanism the optimal distinguisher's ROC is that of
# N(0,1) versus N(mu,1), so its AUC is at most
#
#       AUC <= Phi(mu / sqrt(2))                                   [attained]
#
# This is far sharper than the (e^eps + delta)/(1 + e^eps + delta) accuracy bound
# S6 quotes, which is 0.99995 at eps = 10 and therefore vacuous.  Both are
# reported.  An attack AUC above the mu-GDP bound refutes the claimed mu: it is
# a proof of violation, not a suggestion of one.
#
# What it covers.  The audit separates the public-X routes:
#   Route 1 (gen.syn.dp.ols.public)  audits at AUC 0.47, which is chance.  It
#       privatises sigma^2 through its own Gaussian channel
#       (Delta_sigma2 = 4 B_y^2 K_X/(n - p)), so the residual scale is noised.
#       It is the positive control, showing the audit is calibrated and does not
#       fire on everything.
#   Route 2 (gen.syn.dp.projected.subagg)  audits at 0.51.  synth.response draws
#       the response from the model's own family at (beta_DP, dispersion), so the
#       release depends on the data only through beta_DP, which is what the
#       subsample-and-aggregate proposition requires.  A release that also
#       depended on the data through the residuals would not meet that
#       hypothesis, and would show an AUC flat in epsilon, since such a leak
#       does not pass through the noise.
#
# The script also covers two routes with sensitive X:
#   V4 (gen.syn.dp.full)      sensitive X, Gaussian response.
#   V5 (gen.syn.dp.logistic)  sensitive X, binary response.
# These add three channels the public-X routes do not have (a DP histogram per
# continuous column, a DP categorical release, and a DP Gaussian copula), so the
# panel of attack statistics grows to watch each of them.  See Section B.
#
# Output: results/dp-audit-mugdp.csv, results/dp-audit-mugdp.log
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
source("lib/33-dp-projection-public-x.r")
source("lib/dp-public-bounds.r")                            # dp.public.suffstat.sens
source("lib/66-dp-suffstat.r")                              # Route 3
.dp.full.no.main <- TRUE
source("lib/35-dp-full.r")                                  # V4
# 53- re-sources the whole stack, and 00-shared.r clears globalenv on load, which
# wipes 33-'s generators.  Route 1 would then be undefined, every rep would error
# into the tryCatch, and the audit would have no positive control at all.
# .shared.loaded is 00-shared.r's own re-init guard (see its usage note).
.shared.loaded <- TRUE
.dp.logistic.no.main <- TRUE
source("lib/53-dp-logistic.r")                              # V5
# 53- rm()s .sobol.main and .dp.full.no.main on its way out, so only clear what
# is actually still here.
rm(list = intersect(c(".sobol.main", ".shared.loaded", ".dp.public.no.main",
                      ".dp.full.no.main", ".dp.logistic.no.main"), ls()))

SMOKE <- nzchar(Sys.getenv("AUDIT_SMOKE"))

s.n      <- 1000L        # background sample size
s.delta  <- 1e-6
s.m      <- NULL         # Route 2 block count, NULL applies floor(n / (5 p))
PROBE.R <- 20           # probe row residual, in residual SDs (an admissible row)
v.eps    <- c(1, 2, 5, 10)
s.B      <- 200L         # replications per world per (route, eps)

if(SMOKE) { v.eps <- c(1, 10); s.B <- 30L
            cat("[SMOKE] eps = {1, 10}, B = 30\n") }

# The design is public: fixed covariates, a public coefficient box, a public
# response bound.  Nothing here is a statistic of the sensitive sample.
v.beta   <- c(2, 1.5, -0.8)
s.sigma  <- 2
v.clip.lo <- c(`(Intercept)` = -8, x1 = -8, x2 = -8)
v.clip.hi <- c(`(Intercept)` =  8, x1 =  8, x2 =  8)
s.B.y     <- 120         # public bound on |y|; the probe row at 20 SDs stays inside

# +------------------------------------------------------------------------+
# |  The neighbouring pair                                                 |
# +------------------------------------------------------------------------+
# Background rows are fixed across the whole audit (they are the adversary's
# known context).  Only the probe row row differs between the two worlds, so the
# two datasets are neighbours under one-row replacement.
set.seed(6489L)
v.x1.bg <- rnorm(s.n - 1L)
v.x2.bg <- rnorm(s.n - 1L)
v.y.bg  <- v.beta[[1L]] + v.beta[[2L]] * v.x1.bg + v.beta[[3L]] * v.x2.bg +
           rnorm(s.n - 1L, 0, s.sigma)

# The probe row's covariates are the same in both worlds.  Only its residual moves.
s.cx1 <- 0; s.cx2 <- 0
s.cmu <- v.beta[[1L]] + v.beta[[2L]] * s.cx1 + v.beta[[3L]] * s.cx2

world <- function(w) {
    s.cy <- s.cmu + `if`(w == 1L, PROBE.R * s.sigma, 0)
    data.frame(response = c(v.y.bg, s.cy),
               x1 = c(v.x1.bg, s.cx1),
               x2 = c(v.x2.bg, s.cx2))
}
l.world <- list(`0` = world(0L), `1` = world(1L))

# +------------------------------------------------------------------------+
# |  The adversary's statistics                                            |
# +------------------------------------------------------------------------+
# A real auditor tries every channel and keeps the best, so the attack is a panel
# of statistics and the reported AUC is the maximum over them.  Reporting only
# the residual scale would understate the adversary once that channel is closed.
#
#   sigma   the residual scale of the release.  This is the channel that was wide
#           open: Route 2 reproduced the real sigma exactly.  With a public
#           dispersion it becomes a constant and carries nothing, which is the
#           point of the fix.
#   b0      the released intercept.  The probe row sits at x = (0, 0), so it moves
#           the intercept most.  This is the legitimate, DP-protected channel: it
#           must leak something (beta_DP does depend on the data) but no more than
#           the mu-GDP bound allows.  Closing sigma without checking b0 would be
#           auditing the wrong quantity.
#   b1, b2  the remaining coefficients, for completeness.
attack.stats <- function(d.syn) {
    mo <- tryCatch(lm(response ~ x1 + x2, d.syn), error = function(e) NULL)
    if(is.null(mo)) return(c(sigma = NA_real_, b0 = NA_real_,
                             b1 = NA_real_, b2 = NA_real_))
    v.b <- coef(mo)
    c(sigma = sigma(mo), b0 = v.b[[1L]], b1 = v.b[[2L]], b2 = v.b[[3L]])
}
V.STATS <- c("sigma", "b0", "b1", "b2")

# +------------------------------------------------------------------------+
# |  Section B: the sensitive-X routes, V4 and V5                          |
# +------------------------------------------------------------------------+
# V4 (gen.syn.dp.full) and V5 (gen.syn.dp.logistic) release the covariates too,
# through a DP histogram per column, a DP categorical release, and a DP Gaussian
# copula on the continuous block.  Those are three channels that Routes 1 and 2
# do not have, and an audit that only watched the response would not look at any
# of them.  Section A's panel is therefore not sufficient here.
#
# V4 and V5 appear sound on the axis where Route 2 failed. Both synthesise the
# covariates and draw a fresh response (V4 from N(0, sigma2_DP), V5 from
# Bernoulli(expit(x' beta_DP))), so neither carries the empirical residual law.
# Code that looks sound can still leak, so measure it.
#
# The background has A deliberate gap.  x1's bin 3 ([2, 3)) is empty in the
# background, and the probe row moves into it.  That is the neighbouring pair an
# adversary would choose against the histogram channel.  A left-packed count
# vector would slide every count above an empty interior bin down one slot,
# making the channel's true sensitivity O(n) rather than sqrt(2).  Keeping the
# gap here tests that end to end in the real generator rather
# than in the isolated mechanism (which is what verify-dp-soundness Task 5 does).
S.NB     <- 1000L
V.BR1    <- seq(0, 5, length.out = 6L)      # 5 unit bins; bin 3 is the gap
V.BR2    <- seq(-4, 4, length.out = 9L)     # 8 bins
L.XSPEC  <- list(
    x1 = list(type = "continuous",  breaks = V.BR1),
    x2 = list(type = "continuous",  breaks = V.BR2),
    x3 = list(type = "categorical", levels = c("0", "1")))
V.BETA.B <- c(2, 1.5, -0.8, 0.6)            # intercept, x1, x2, x3
S.SIG.B  <- 2
S.CLIP   <- 8                               # public scalar box, +/- 8 on every coord
V.DISP   <- c(lo = 0, hi = 40)              # public box on sigma^2 (true value 4)

set.seed(9317L)
# x1 avoids [2, 3) entirely, so bin 3 is empty until the probe row arrives.
v.b.x1 <- c(runif(ceiling((S.NB - 1L) / 2), 0, 2),
            runif(floor((S.NB - 1L) / 2),   3, 5))
# x2 is correlated with x1, so the copula channel has structure to protect.
v.b.x2 <- pmin(pmax(0.9 * (v.b.x1 - 2.5) + rnorm(S.NB - 1L, 0, 1), -3.9), 3.9)
v.b.x3 <- factor(rbinom(S.NB - 1L, 1L, 0.4), levels = c(0L, 1L),
                 labels = c("0", "1"))
v.b.eta <- V.BETA.B[[1L]] + V.BETA.B[[2L]] * v.b.x1 + V.BETA.B[[3L]] * v.b.x2 +
           V.BETA.B[[4L]] * (v.b.x3 == "1")
v.b.y   <- v.b.eta + rnorm(S.NB - 1L, 0, S.SIG.B)
v.b.yb  <- rbinom(S.NB - 1L, 1L, plogis(v.b.eta - 2))     # binary response for V5

# The probe row is adversarial in every channel at once, which is the full-row
# replacement V4/V5 claim to protect against:
#   world 0   x1 = 0.5 (a populated bin), x2 = 0, x3 = "0", y on the line
#   world 1   x1 = 2.5 (the empty bin),  x2 = 3.9 (the corner, against a
#             positively correlated background), x3 = "1", y 20 residual SDs off
world.B <- function(w, s.y.nm) {
    s.x1 <- `if`(w == 1L, 2.5, 0.5)
    s.x2 <- `if`(w == 1L, 3.9, 0.0)
    s.x3 <- `if`(w == 1L, "1", "0")
    d <- data.frame(
        x1 = c(v.b.x1, s.x1), x2 = c(v.b.x2, s.x2),
        x3 = factor(c(as.character(v.b.x3), s.x3), levels = c("0", "1")))
    s.eta <- V.BETA.B[[1L]] + V.BETA.B[[2L]] * s.x1 + V.BETA.B[[3L]] * s.x2 +
             V.BETA.B[[4L]] * (s.x3 == "1")
    d[[s.y.nm]] <- `if`(s.y.nm == "response",
        c(v.b.y,  s.eta + `if`(w == 1L, PROBE.R * S.SIG.B, 0)),
        c(v.b.yb, `if`(w == 1L, 1L, 0L)))
    d
}
l.world.B  <- list(`0` = world.B(0L, "response"), `1` = world.B(1L, "response"))
l.world.B5 <- list(`0` = world.B(0L, "ybin"),     `1` = world.B(1L, "ybin"))

# The panel now watches the covariate channels too.  Each statistic is aimed at
# one of them, so a violation names the channel it came through.
#   gapfrac  share of released x1 in [2, 3), the probe row's bin.  The histogram
#            Channel, and the statistic that fires on a left-packed vector.
#   cor12    correlation of the two continuous covariates in the release.  The
#            Copula channel.
#   p.x3     share of x3 == "1".  The categorical channel.
#   x2.max   an upper-tail statistic on x2's marginal, where the probe row sits.
#   sigma    residual scale (V4 only).  Route 2's old leak.  Must be closed here.
#   b0, b1   released coefficients.  The legitimate, DP-protected channel: it
#            must leak something, but no more than the mu-GDP bound allows.
x.stats <- function(d.syn) c(
    gapfrac = mean(d.syn[["x1"]] >= 2 & d.syn[["x1"]] < 3),
    cor12   = suppressWarnings(stats::cor(d.syn[["x1"]], d.syn[["x2"]])),
    p.x3    = mean(as.character(d.syn[["x3"]]) == "1"),
    x2.max  = max(d.syn[["x2"]]))

attack.stats.B <- function(d.syn) {
    mo <- tryCatch(lm(response ~ x1 + x2 + x3, d.syn), error = function(e) NULL)
    v.b <- `if`(is.null(mo), c(NA_real_, NA_real_),
                c(coef(mo)[[1L]], coef(mo)[[2L]]))
    c(sigma = `if`(is.null(mo), NA_real_, sigma(mo)),
      b0 = v.b[[1L]], b1 = v.b[[2L]], x.stats(d.syn))
}
V.STATS.B <- c("sigma", "b0", "b1", "gapfrac", "cor12", "p.x3", "x2.max")

attack.stats.B5 <- function(d.syn) {
    mo <- tryCatch(suppressWarnings(
              glm(ybin ~ x1 + x2 + x3, d.syn, family = binomial())),
          error = function(e) NULL)
    v.b <- `if`(is.null(mo), c(NA_real_, NA_real_),
                c(coef(mo)[[1L]], coef(mo)[[2L]]))
    c(ybar = mean(d.syn[["ybin"]]), b0 = v.b[[1L]], b1 = v.b[[2L]],
      x.stats(d.syn))
}
V.STATS.B5 <- c("ybar", "b0", "b1", "gapfrac", "cor12", "p.x3", "x2.max")

# The dispersion is declared public: the DGP draws N(0, 2), so sigma^2 = 4 is a
# constant of the design.  This is the mode that costs no budget.
disp.pub <- dp.dispersion.public(4,
    source = "DGP constant: the audit's own design draws residuals N(0, sd = 2)")

gen.route2 <- function(d, eps, seed)
    gen.syn.dp.projected.subagg(d, y.nm = "response", family = "gaussian",
        link = "identity", x.nms = c("x1", "x2"), eps = eps, delta = s.delta,
        m = s.m, clip.lo = v.clip.lo, clip.hi = v.clip.hi,
        dispersion = disp.pub, seed = seed)

gen.route1 <- function(d, eps, seed)
    gen.syn.dp.ols.public(d, y.nm = "response", x.nms = c("x1", "x2"),
        eps = eps, delta = s.delta, B.y = s.B.y, seed = seed)

# Route 3: the whitened sufficient-statistic release, audited under the
# same Section A worlds as Routes 1 and 2 (the d.syn emission publishes X, so
# the y-only probe row is the neighbouring pair the emission's premise covers).
# The audit design's public constants are declared inline, exactly as disp.pub
# is: the population second moment of (1, x1, x2) under the design's own law is
# the identity, so the whitening is the identity map, B_x is the box norm, and
# the response bound is the audit's own s.B.y.  The lambda rule is
# released-eigenmin, the deployment mode for a design with no published
# error curve, which also puts the fourth (eig) channel under audit.  The
# dispersion stays private (the y'y channel), so all channels are exercised.
l.w3 <- list(
    W = `dimnames<-`(diag(3), list(c("(Intercept)", "x1", "x2"),
                                   c("(Intercept)", "x1", "x2"))),
    B.x = sqrt(1 + 2 * 8^2),
    source = paste0("identity whitening: the audit design draws x1, x2 from ",
                    "N(0, 1), whose population second moment is the identity; ",
                    "B_x = max ||x|| over the public box |x_j| <= 8, a design ",
                    "constant"))
l.y3 <- list(centre = 0, scale = 1, B.y = s.B.y, y.clip = c(-s.B.y, s.B.y),
             source = paste0("the audit's own public response bound |y| <= ",
                             s.B.y, ", identity standardisation"))
l.s3 <- dp.public.suffstat.sens(l.w3, l.y3)

gen.suffstat <- function(d, eps, seed)
    gen.syn.dp.suffstats(d, y.nm = "response", sc = NA_integer_, eps = eps,
        delta = s.delta, whitening = l.w3, y.std = l.y3, sens = l.s3,
        lambda.method = "released-eigenmin", x.public = TRUE, seed = seed)

# V4: sensitive X, Gaussian response.  Every channel is a Gaussian mechanism and
# they compose in one mu-GDP quadrature (mu.split, default 5/5/90 on marg/corr/coef).
gen.v4 <- function(d, eps, seed)
    gen.syn.dp.full(d, y.nm = "response", x.spec = L.XSPEC,
        eps.total = eps, delta = s.delta, m = s.m,
        clip.lo = -S.CLIP, clip.hi = S.CLIP,
        disp.lo = V.DISP[["lo"]], disp.hi = V.DISP[["hi"]],
        use.copula = TRUE, score.method = "adaptive", seed = seed)

# V5: sensitive X, binary response.  The Bernoulli draw from expit(x' beta_DP) is
# post-processing (supplement S6, third extension).  No dispersion, so the
# aggregate stays p-dimensional.
gen.v5 <- function(d, eps, seed)
    gen.syn.dp.logistic(d, y.nm = "ybin", x.spec = L.XSPEC,
        x.public = FALSE, eps = eps, delta = s.delta, m = s.m,
        clip.lo = -S.CLIP, clip.hi = S.CLIP,
        use.copula = TRUE, score.method = "adaptive", seed = seed)

# Mann-Whitney AUC: P(T_1 > T_0), the optimal distinguisher's success rate.
#
# The index is seq_along(v.t1), not seq_along(v.t0).  The latter picks out the
# right elements only when the two worlds have the same number of usable reps,
# which they need not, since each world is filtered by is.finite independently.
# One failed generator rep would shift the window, and the AUC would then be
# computed from a mix of the two worlds' ranks.  V5's per-block glm does fail
# occasionally, so this matters.
mw.auc <- function(v.t0, v.t1) {
    s.r <- rank(c(v.t0, v.t1))
    (sum(s.r[seq_along(v.t1) + length(v.t0)]) -
         length(v.t1) * (length(v.t1) + 1) / 2) / (length(v.t0) * length(v.t1))
}

S.BOOT  <- 2000L        # bootstrap resamples for the advantage's lower bound
S.ALPHA <- 0.05         # one-sided level, Bonferroni-split across the channels

# Returns the point AUC and a lower confidence bound on the advantage |AUC - 0.5|.
#
# Why A confidence bound, and not the point estimate.  The headline advantage is a
# Max over the K channels the adversary tried, and the max of K noisy estimates is
# biased upward by construction: under a perfectly sound mechanism every channel
# sits at 0.5 +/- SE, and the largest of seven such deviations lands well above 0.5
# on noise alone.  At B = 30 per world (SE ~ 0.09) that artefact alone produces a
# max AUC near 0.62, which the eps = 1 bound of 0.566 would report as a violation.
# An audit that reports spurious violations is worth no more than one that passes
# everything.  A violation here is a proof that the mechanism does not deliver
# its stated mu, so it has to be a claim that withstands review.
#
# The verdict is therefore taken on the lower bound, Bonferroni-corrected across the K
# channels.  This is the standard construction in the DP-auditing literature
# (Jagielski et al. 2020, Nasr et al. 2021).  Power against a real leak is
# retained, since an AUC of 1.000 has a lower bound of 1.000.
#
# The bootstrap resamples each world independently and is seeded from the label, so
# the whole audit stays reproducible even though the generators reset the global RNG
# on every call.
#
# Note the hard stop.  An earlier version returned NA when the generator failed, and
# the verdict line read isTRUE(NA > bound) == FALSE, i.e. "within bound".  A run that
# produced no data reported as a pass.  That is the same silent-failure pattern that
# let the V4/V5 drivers drop every replication and write empty results, and an audit
# that passes when it measured nothing is worse than no audit.  Fail loudly instead.
auc.of <- function(v.t0, v.t1, s.label, s.k = 1L) {
    v.t0 <- v.t0[is.finite(v.t0)]; v.t1 <- v.t1[is.finite(v.t1)]
    if(length(v.t0) < 5L || length(v.t1) < 5L)
        stop(sprintf("audit produced no usable reps for %s (%d / %d finite). ",
                     s.label, length(v.t0), length(v.t1)),
             "The generator is erroring; fix it rather than reporting a pass.")
    # A degenerate (constant) statistic is a legitimate result, not a failure: it
    # means the channel carries nothing.  Chance is exactly 0.5, advantage 0.
    if(length(unique(c(v.t0, v.t1))) == 1L) return(c(auc = 0.5, adv.lo = 0))
    s.auc <- mw.auc(v.t0, v.t1)
    set.seed(sum(utf8ToInt(s.label)))
    v.bs <- vapply(seq_len(S.BOOT), function(b)
        mw.auc(sample(v.t0, replace = TRUE), sample(v.t1, replace = TRUE)),
        numeric(1))
    c(auc = s.auc,
      adv.lo = unname(quantile(abs(v.bs - 0.5), S.ALPHA / s.k, names = FALSE)))
}

# +------------------------------------------------------------------------+
# |  Run                                                                   |
# +------------------------------------------------------------------------+
cat("EMPIRICAL mu-GDP AUDIT, probe row distinguishing game on a neighbouring pair\n")
cat(strrep("=", 82), "\n")
cat(sprintf("n = %d, probe row residual = %g SDs, B = %d reps/world, delta = %g\n\n",
            s.n, PROBE.R, s.B, s.delta))
cat(sprintf("%-8s %5s %7s | %8s %8s %8s | %9s | %s\n",
            "route", "eps", "mu", "AUC(max)", "AUC lo", "via", "AUC bound",
            "verdict"))
cat(strrep("-", 92), "\n")

# One engine, five routes.  A route is (worlds, generator, statistic panel).  The
# distinguishing game is identical for all of them, and that is the point: the
# mu-GDP bound does not care which mechanism produced the release.
audit.route <- function(s.route, gen, l.w, v.stats, stats.fn) {
    Map(function(s.eps) {
        s.mu <- gdp.mu.from.eps.delta(s.eps, s.delta)

        # Independent mechanism randomness in each world, as the adversary faces.
        #
        # The offset must depend on eps.  A first version reused seeds 10000+ and
        # 90000+ at every eps, so the same standard normals were recycled across
        # the whole grid and the AUC was identical to four decimals at every
        # eps and for both routes -- an artefact masquerading as a result.  The
        # epsilon comparison is only meaningful if the mechanism's randomness is
        # resampled independently at each point.
        s.off <- 10000L * as.integer(match(s.eps, v.eps))
        draw <- function(w, off) Map(function(b) {
            res <- tryCatch(gen(l.w[[w]], s.eps, off + b), error = function(e) NULL)
            `if`(is.null(res), `names<-`(rep(NA_real_, length(v.stats)), v.stats),
                 stats.fn(res[["d.syn"]])[v.stats])
        }, seq_len(s.B)) |> Reduce(rbind, x = _)
        m.t0 <- draw("0", s.off)
        m.t1 <- draw("1", s.off + 500000L)

        # The auditor tries every channel and keeps the best, so the level is
        # Bonferroni-split across them.
        s.K <- length(v.stats)
        m.a <- vapply(v.stats, function(k)
            auc.of(m.t0[, k], m.t1[, k],
                   sprintf("%s/eps=%g/%s", s.route, s.eps, k), s.K),
            numeric(2))
        v.auc <- m.a["auc", ]; v.lo <- m.a["adv.lo", ]

        # A distinguisher may bet either way, so the advantage is |AUC - 0.5|.
        # The point estimate names the channel and is reported.  The verdict is
        # taken on the lower confidence bound (see auc.of).
        s.best   <- v.stats[[which.max(abs(v.auc - 0.5))]]
        s.auc    <- 0.5 + max(abs(v.auc - 0.5))
        s.auc.lo <- 0.5 + max(v.lo)

        s.auc.bd <- pnorm(s.mu / sqrt(2))                    # exact mu-GDP bound
        b.viol   <- isTRUE(s.auc.lo > s.auc.bd + 1e-9)

        cat(sprintf("%-8s %5g %7.3f | %8.3f %8.3f %8s | %9.3f | %s\n",
                    s.route, s.eps, s.mu, s.auc, s.auc.lo, s.best, s.auc.bd,
                    `if`(b.viol, "*** VIOLATES mu-GDP ***", "within bound")))
        # Every channel's AUC is kept, not only the winner: a route can be within
        # bound overall and still be the wrong shape underneath.
        data.frame(route = s.route, eps = s.eps, mu = s.mu, B = s.B,
                   auc = s.auc, auc.lo = s.auc.lo, best.stat = s.best,
                   auc.bound = s.auc.bd, violates = b.viol,
                   stat = names(v.auc), auc.stat = as.numeric(v.auc),
                   adv.lo.stat = as.numeric(v.lo),
                   row.names = NULL)
    }, v.eps)
}

l.rows <- c(
    # Section A: public X.  Route 1 is the positive control (it noises sigma^2
    # through its own channel), Route 2 the subsample-and-aggregate mechanism,
    # and Route 3 the suffstat release with all its channels private.
    audit.route("public",   gen.route1, l.world,    V.STATS,    attack.stats),
    audit.route("subagg",   gen.route2, l.world,    V.STATS,    attack.stats),
    audit.route("suffstat", gen.suffstat, l.world,  V.STATS,    attack.stats),
    # Section B: sensitive X.  Three extra channels, and the background carries an
    # empty interior bin that the probe row moves into (see above).
    audit.route("v4-full",  gen.v4,     l.world.B,  V.STATS.B,  attack.stats.B),
    audit.route("v5-logit", gen.v5,     l.world.B5, V.STATS.B5, attack.stats.B5))

d.out <- Reduce(rbind, l.rows)

cat("\n")
cat("AUC bound = Phi(mu / sqrt(2)), the exact optimal-distinguisher bound for a\n")
cat("mu-GDP mechanism.  An AUC above it is a PROOF that the release does not\n")
cat("satisfy the mu it claims.  The accuracy bound is the looser S6 form, shown\n")
cat("for comparison: at eps = 10 it is 0.99995 and cannot detect anything.\n")

# Per-channel breakdown for the sensitive-X routes.  The headline AUC is a max
# over channels, so it says whether the release is within bound and not through
# What.  V4/V5 add the histogram, categorical and copula channels, and each is
# reported here so a future violation names its own channel.
cat("\nPer-channel AUC, sensitive-X routes (v4-full, v5-logit):\n")
d.ch <- d.out[d.out[["route"]] %in% c("v4-full", "v5-logit"), ]
v.ch <- unique(d.ch[["stat"]])
cat(sprintf("%-9s %5s | %s\n", "route", "eps",
            paste(sprintf("%8s", v.ch), collapse = " ")))
invisible(Map(function(s.r) Map(function(s.e) {
    d.i <- d.ch[d.ch[["route"]] == s.r & d.ch[["eps"]] == s.e, ]
    cat(sprintf("%-9s %5g | %s\n", s.r, s.e,
        paste(sprintf("%8.3f", d.i[["auc.stat"]][match(v.ch, d.i[["stat"]])]),
              collapse = " ")))
}, v.eps), c("v4-full", "v5-logit")))
cat("\ngapfrac is the share of released x1 in the probe row's bin [2, 3), which is\n")
cat("EMPTY in the background.  It is the histogram channel, and the statistic that\n")
cat("fires on a left-packed count vector, where an empty interior bin would\n")
cat("shift every count above it.\n")

s.viol <- d.out[!duplicated(d.out[c("route", "eps")]) & d.out[["violates"]], ]
if(nrow(s.viol) > 0L)
    stop(sprintf("mu-GDP VIOLATED on %d (route, eps) cells: %s",
                 nrow(s.viol),
                 paste(sprintf("%s@eps=%g", s.viol[["route"]], s.viol[["eps"]]),
                       collapse = ", ")))

s.out <- `if`(SMOKE, "results/dp-audit-mugdp-smoke.csv",
                     "results/dp-audit-mugdp.csv")
write.csv(d.out, s.out, row.names = FALSE)
cat(sprintf("\nWritten to %s (%d rows)\n", sub("^\\.\\./", "", s.out), nrow(d.out)))
