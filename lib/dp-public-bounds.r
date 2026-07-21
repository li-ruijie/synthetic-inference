#!/usr/bin/env r
# +------------------------------------------------------------------------+
# |  dp-public-bounds.r -- the only permitted source of a DP bound         |
# +------------------------------------------------------------------------+
#
# Why this file exists.  Eight times now this project has calibrated DP noise to a
# sensitivity that was not the sensitivity of the function it released, and four of
# those were the same line, written independently in different files:
#
#     B.y <- 1.5 * max(abs(d[["response"]]))
#
# It appeared in twenty-four places across the code tree, including inside three
# Generators (28-, 33-, 34-), which silently substituted it whenever a caller
# omitted B.y.  A driver author could not see the substitution happening.  The
# bound looks like a modelling choice and is in fact a statistic of the very data
# being protected: it sets Delta_beta = 2 B_y max_i ||(X'X)^-1 x_i||, hence the
# noise scale, hence there is no valid sensitivity bound, hence no theorem.
#
# The fix is not to correct the twenty-four lines.  It is to make the incorrect line
# impossible to write by accident:
#
#   1. The three generators refuse a NULL B.y rather than defaulting it.
#   2. Every bound a driver uses comes from here, and every one carries a
#      Required `source` string that is printed in the run's log.
#
# The point of the `source` field is that it is self-indicting.  The unsound bound
# still exists as a construction, but to use it you now have to type
#
#     list(B.y = 1.5 * max(abs(y)), source = "the sample maximum")
#
# which no one will do, and which any reader of the log will catch.  A design that
# makes the unsound call unwriteable without saying so is worth more than any
# one fix.
#
# Everything here is A functional of the published DGP, computed by
# 64-dp-design-constants.r from an independent draw of N = 1e6 at the public seed
# 20260713.  None of it touches the n = 2000, seed-6489 analysis sample.
#
# Do not source A box from supplement S1.  Its column is headed "Original"
# under the note "Seed 6489, n = 2000": it is coef(lm(...)) on the analysis sample.
# Building a box from it would repeat the same error, carried through the
# paper's own supplement, where an auditor tracing the provenance reads
# "Supplement S1" and stops.  The population coefficient differs from S1's by up to
# 11.6%, so the two are not interchangeable.
#
# Requires: dp.dispersion.public() from 28-dp-projection.r (source that first).
# Reads:    results/dp-design-constants.rds (run 64-dp-design-constants.r first).
#
# +------------------------------------------------------------------------+

.DP.CONST.PATH <- "results/dp-design-constants.rds"

.dp.const <- function(sc) {
    if(!file.exists(.DP.CONST.PATH))
        stop("dp-public-bounds: ", .DP.CONST.PATH, " is missing. ",
             "Run 64-dp-design-constants.r first: it computes the population ",
             "constants from the published DGP, and every DP box rests on them.")
    l.k <- readRDS(.DP.CONST.PATH)
    s.nm <- paste0("sc", sc)
    if(is.null(l.k[[s.nm]]))
        stop(sprintf("dp-public-bounds: no design constants for scenario %s", sc))
    l.k[[s.nm]]
}

# +------------------------------------------------------------------------+
# |  The public covariate metadata                                         |
# +------------------------------------------------------------------------+
# The x.spec breaks and levels are public by assumption: the DP histogram and
# categorical mechanisms both require public bin edges, and their sqrt(2)
# sensitivity is a statement about a fixed, public binning.  Defining them once
# stops a driver from quietly widening a break to fit its own sample.
dp.public.xspec <- function()
    list(age    = list(type = "continuous",  breaks = seq(20, 90, length = 11)),
         female = list(type = "categorical", levels = c("0", "1")),
         race   = list(type = "categorical", levels = c("1", "2", "3")),
         income = list(type = "continuous",  breaks = seq(0, 4000, length = 11)))

# +------------------------------------------------------------------------+
# |  The public coefficient box                                            |
# +------------------------------------------------------------------------+
# +/- kappa * |beta_pop|, centred on zero (the approved construction, and the
# one 62-mc-dp-gdp.r already uses).  Centring on zero rather than on beta_pop
# costs nothing, since the width sets the noise (2 kappa |beta_pop_j|, identical
# either way), and it assumes strictly less: it claims a magnitude scale and not
# the sign or the location.  Take the conservatism at no cost.
#
# v.names lets the caller align the box to the coefficient names its own model
# matrix produces, which is what subsample.aggregate.beta.loglink expects.
dp.public.box <- function(sc, kappa, v.names = NULL) {
    l.k <- .dp.const(sc)
    v.b <- l.k[["beta"]]
    if(!is.null(v.names)) {
        if(!all(v.names %in% names(v.b)))
            stop(sprintf(paste0("dp-public-bounds: scenario %s has no population ",
                                "coefficient for %s. The box must cover every ",
                                "coefficient the model fits."),
                         sc, paste(setdiff(v.names, names(v.b)), collapse = ", ")))
        v.b <- v.b[v.names]
    }
    list(lo = -kappa * abs(v.b), hi = kappa * abs(v.b),
         source = sprintf("+/- %g * |beta_pop|, centred on zero; %s",
                          kappa, l.k[["source"]]))
}

# The sigma^2 box for the subagg dispersion coordinate (V4, disp = TRUE).
# sigma^2 is non-negative, so the box is anchored at zero rather than centred:
# [0, kappa * sigma2_pop].  Its width sets s_{p+1} = kappa * sigma2_pop / m.
dp.public.disp.box <- function(sc, kappa) {
    l.k <- .dp.const(sc)
    s.s2 <- l.k[["sigma2"]]
    if(is.na(s.s2))
        stop(sprintf(paste0("dp-public-bounds: scenario %s (%s) has NO free ",
                            "dispersion -- its variance is fixed by the mean -- ",
                            "so it needs no sigma^2 box and cannot leak one."),
                     sc, l.k[["family"]]))
    list(lo = 0, hi = kappa * s.s2,
         source = sprintf("[0, %g * sigma2_pop] = [0, %.4g]; %s",
                          kappa, kappa * s.s2, l.k[["source"]]))
}

# The dispersion as a public value, for the Route 2 generators, which take a
# dispersion object rather than a box, dispersion being a mode.  This is the
# free mode: the population residual variance is a functional of the published
# DGP, so re-using it spends no budget.
dp.public.dispersion <- function(sc) {
    l.k <- .dp.const(sc)
    s.s2 <- l.k[["sigma2"]]
    if(is.na(s.s2))
        stop(sprintf(paste0("dp-public-bounds: scenario %s (%s) has no free ",
                            "dispersion; do not pass one."), sc, l.k[["family"]]))
    dp.dispersion.public(s.s2, source = sprintf("population dispersion; %s",
                                                l.k[["source"]]))
}

# The logistic zero-model's box, for the ZI generator's second channel.  Only the
# two-part (quasi-gamma) scenarios have one.
dp.public.zero.box <- function(sc, kappa, v.names = NULL) {
    l.k <- .dp.const(sc)
    v.g <- l.k[["gamma"]]
    if(is.null(v.g))
        stop(sprintf(paste0("dp-public-bounds: scenario %s (%s) is not ",
                            "zero-inflated and has no zero-model box."),
                     sc, l.k[["family"]]))
    if(!is.null(v.names)) v.g <- v.g[v.names]
    list(lo = -kappa * abs(v.g), hi = kappa * abs(v.g),
         source = sprintf("+/- %g * |gamma_pop| (logistic zero-model), zero-centred; %s",
                          kappa, l.k[["source"]]))
}

# +------------------------------------------------------------------------+
# |  The public block-spread scale (two-stage box floor)                     |
# +------------------------------------------------------------------------+
# tau_j(m) = sqrt(v_j * m / n), with v_j the per-fitted-observation coefficient
# variance published by 64- (diag(vcov) x nobs of the population fit).  A block
# fit on b rows has sampling variance close to v_j / b, so tau_j(m) is the
# public scale of the block-estimate spread that a re-centred coefficient box
# must cover.  It is a functional of the published DGP and of the public design
# choices (n, m), not a statistic of the analysis sample.
#
# channel "beta" is the family fit.  On the two-part (quasi-gamma) scenarios
# that fit uses the positive rows only, so the per-block fitted count is
# (n / m) * pos.frac with the public positive fraction, and v_j is already
# per-positive-observation.  channel "gamma" is the ZI logistic zero-model,
# fitted on every row of a block.
dp.public.block.tau <- function(sc, m, n, v.names = NULL,
                                channel = c("beta", "gamma")) {
    channel <- match.arg(channel)
    l.k <- .dp.const(sc)
    v.v1 <- l.k[[`if`(channel == "beta", "beta.var1", "gamma.var1")]]
    if(is.null(v.v1))
        stop(sprintf(paste0("dp-public-bounds: scenario %s (%s) has no %s ",
                            "block-spread constants.  Re-run ",
                            "64-dp-design-constants.r, or the scenario is not ",
                            "zero-inflated and has no zero-model."),
                     sc, l.k[["family"]], channel))
    if(!is.null(v.names)) {
        if(!all(v.names %in% names(v.v1)))
            stop(sprintf(paste0("dp-public-bounds: scenario %s has no ",
                                "block-spread constant for %s."),
                         sc, paste(setdiff(v.names, names(v.v1)),
                                   collapse = ", ")))
        v.v1 <- v.v1[v.names]
    }
    s.frac <- `if`(channel == "beta" && !is.null(l.k[["pos.frac"]]),
                   l.k[["pos.frac"]], 1)
    s.b <- (n / m) * s.frac
    list(tau = sqrt(v.v1 / s.b),
         source = sprintf(paste0("tau_j(m) = sqrt(v_j / b), b = %.4g fitted ",
                                 "rows per block (n/m = %.4g, fitted fraction ",
                                 "%.4g), v_j = diag(vcov) x nobs of the ",
                                 "population %s fit; %s"),
                          s.b, n / m, s.frac, channel, l.k[["source"]]))
}

# +------------------------------------------------------------------------+
# |  The published zero rate (block-count rule)                             |
# +------------------------------------------------------------------------+
# The ZI positive channel fits the positive rows only, so its block-count rule
# scales n by the published positive fraction, n_eff = n * pos.frac
# = n * (1 - pi0).  pos.frac is a functional of the public draw (64-), and the
# rate travels as a dp.rate provenance object so a realised rate cannot slip
# in.
dp.public.zero.rate <- function(sc) {
    l.k <- .dp.const(sc)
    if(is.null(l.k[["pos.frac"]]))
        stop(sprintf(paste0("dp-public-bounds: scenario %s (%s) is not ",
                            "zero-inflated and has no published zero rate"),
                     sc, l.k[["family"]]))
    dp.rate.public(1 - l.k[["pos.frac"]],
                   source = sprintf("published zero rate 1 - pos.frac; %s",
                                    l.k[["source"]]))
}

# +------------------------------------------------------------------------+
# |  The public response bound, and the clipping step it implies (Route 1)  |
# +------------------------------------------------------------------------+
# Route 1 only.  Subsample-and-aggregate needs no bound on |y|: its sensitivity
# comes from clipping the per-block coefficients into a public box, so it holds
# however large or unbounded y is.  That is the reason the mechanism exists.
# Route 1 is different, and it is a theorem: Delta_beta = 2 B_y max_i
# ||(X'X)^-1 x_i|| is defined through B_y, because OLS beta_hat has unbounded
# global sensitivity when y is unbounded (one row at |y| -> infinity moves
# beta_hat arbitrarily).  No public B_y, no Route 1.
#
# Scenarios 1, 2 and 10 draw y = qnorm(u, 50, 15), which is genuinely unbounded,
# so no finite B_y exists for them and 1.5 * max|y| made the bounded-universe
# assumption true for the observed sample by construction.  That is the whole
# circularity.  A sound Route 1 there needs y clipped to a public range as an
# explicit pre-processing step, which moves the estimand to the clipped-OLS one
# and must be declared as such in the public-covariate exact release
# proposition.
dp.public.B.y <- function(sc) {
    l.k <- .dp.const(sc)
    if(is.na(l.k[["B.y"]]))
        stop(sprintf(paste0("dp-public-bounds: scenario %s has no public bound on ",
                            "|y|, and none can be invented from the sample. Route 1 ",
                            "is UNAVAILABLE here. Use Route 2 (subagg), which needs ",
                            "no bound on |y|."), sc))
    list(B.y = l.k[["B.y"]], y.clip = l.k[["y.clip"]], source = l.k[["y.source"]])
}

# Apply the public clipping step.  A no-op where the DGP is structurally bounded
# (scenario 9, y in [0, 100]).  At mean +/- 6 sd the clip binds with probability
# ~1e-9 per row, so on any realised sample it is the identity map, but the
# Estimand it defines is clipped-OLS and the paper must say so.
dp.public.clip.y <- function(d, y.nm, sc) {
    v.c <- dp.public.B.y(sc)[["y.clip"]]
    if(is.null(v.c)) return(d)
    `[[<-`(d, y.nm, value = pmin(pmax(d[[y.nm]], v.c[[1L]]), v.c[[2L]]))
}

# +------------------------------------------------------------------------+
# |  The public covariate bound                                            |
# +------------------------------------------------------------------------+
# The Gaussian-mechanism-on-X'y baseline has sensitivity Delta = 2 B_x B_y under
# one-row replacement, where B_x = max_i ||x_i||_2.  Taking either factor from
# the analysis sample, as in
#
#     B.x <- sqrt(max(rowSums(m.X^2)))          # the sample's max row norm
#     B.y <- 1.5 * max(abs(d[["response"]]))    # the sample maximum
#
# would make the "global" sensitivity a function of the very data it protects,
# voiding the epsilon reported.  The paper's central negative result is measured
# against this baseline, so its bounds must come from here.
#
# B_x is a functional of the public support box that dp.public.xspec() already
# declares, so no statistic of the sample is needed.  Over that box:
#
#   intercept    1                 -> 1
#   age          [20, 90]          -> 90^2
#   female       one 0/1 dummy     -> 1
#   race         two 0/1 dummies, mutually exclusive under treatment coding -> 1
#   income       [0, 4000]         -> 4000^2
#
# The premise this makes explicit: the bound holds only for rows inside the public
# box, so the baseline presumes X is publicly clipped to it, exactly as Route 1's
# B_y holds only alongside dp.public.clip.y.  Both are the same bounded-universe
# assumption, declared rather than assumed.
dp.public.B.x <- function(intercept = TRUE) {
    l.spec <- dp.public.xspec()
    v.sq <- vapply(l.spec, function(s) {
        if(identical(s[["type"]], "continuous")) {
            v.r <- range(s[["breaks"]])
            max(abs(v.r))^2
        } else {
            # k - 1 treatment dummies, each in {0, 1}, at most one of them 1.
            1
        }
    }, numeric(1))
    s.sq <- sum(v.sq) + if(isTRUE(intercept)) 1 else 0
    list(B.x = sqrt(s.sq),
         source = paste0("max ||x||_2 over the PUBLIC support box of ",
                         "dp.public.xspec(): age in [20, 90], income in [0, 4000], ",
                         "dummies in {0, 1}. A functional of the published design, ",
                         "not a statistic of the analysis sample."))
}

# The PHI application's public response bound.  Waist circumference has a
# physiological range, which is exactly what a public B_y is supposed to look
# like: external domain knowledge, fixed before the data are seen, and no more a
# statistic of the sample than the speed of light is.
dp.public.B.y.phi <- function()
    list(B.y = 250, y.clip = c(0, 250),
         source = paste0("adult waist circumference, physiological range ",
                         "[0, 250] cm. External domain knowledge, fixed before ",
                         "the data are seen, not a statistic of the sample."))

# +------------------------------------------------------------------------+
# |  B_x on the normalised scale, in closed form                           |
# +------------------------------------------------------------------------+
# Map every continuous covariate to [0, 1] by its public range and leave the
# dummies in {0, 1}.  Under treatment coding a factor activates at most one
# dummy, so any row of the public universe satisfies
#
#     ||x||^2 <= 1 (intercept) + 1 per continuous covariate + 1 per factor
#
# and therefore B_x = sqrt(1 + n.cont + n.factor), exactly, with no reference to
# any sample.  This is the smallest sound B_x there is, so it gives the Gaussian
# mechanism the most favourable Delta_2 it can soundly be given.
#
# 20-dp-precheck.r was reaching for it when it min-max scaled X and y by their
# Sample ranges and then declared "y in [0, 1], so B_y = 1".  That made the
# bounded-universe assumption true for the observed sample by construction, which
# is the circularity this whole file exists to prevent.  Two neighbouring
# datasets get two different scalings, hence two different noise scales, hence no
# sensitivity bound holds and no theorem can be instantiated.  Normalising X is
# harmless (the mechanism already treats X as public -- it inverts X'X in the
# clear).  Normalising Y by its own sample range is unsound.
dp.public.B.x.normalised <- function(n.cont, n.factor, intercept = TRUE) {
    s.sq <- n.cont + n.factor + if(isTRUE(intercept)) 1 else 0
    list(B.x = sqrt(s.sq),
         source = sprintf(paste0("sqrt(1 + %d continuous + %d factors) on the ",
                                 "unit-normalised scale: treatment coding ",
                                 "activates at most one dummy per factor, so this ",
                                 "bounds ||x||_2 over the whole public universe. ",
                                 "A functional of the DESIGN, not of the sample."),
                          n.cont, n.factor))
}

# B_x from a design matrix the mechanism already treats as public.
#
# Read this before calling it.  `sqrt(max(rowSums(m.X^2)))` is, syntactically,
# the exact line found in the baseline driver and removed.  It is sound here
# and nowhere else, and the reason is a premise, not a technicality:
#
#   The Gaussian-on-X'y baseline inverts X'X in the clear.  It cannot run
#   otherwise, so it assumes X is public, which is the same assumption the
#   public-covariate exact release proposition makes ("Assume X is public"),
#   and which is what makes the baseline a like-for-like comparison rather
#   than a baseline handed an advantage the sound routes are denied.  Under
#   that premise every functional of X -- its row norms, its column ranges,
#   X'X, its eigenvalues -- is public.  The only sensitive object is y.
#
# A bound read off X is acceptable, and a bound read off y is unsound.  That is the
# whole distinction, and it is why B.y has no analogue of this function and not
# will.  If your X is sensitive, this bound does not exist and neither does the
# mechanism: Delta_beta is a function of X, so the sensitivity constant would
# itself leak the design (the author decision), and the global-sensitivity
# route is unavailable even with a perfectly known B_y.  Use Route 2.
dp.public.B.x.design <- function(m.X, source) {
    if(missing(source) || !is.character(source) || !nzchar(source))
        stop("dp-public-bounds: B_x from a design matrix REQUIRES a `source` ",
             "saying WHY that design is public. If you cannot write the ",
             "sentence, the bound is not public and you do not have a mechanism.")
    list(B.x = sqrt(max(rowSums(m.X^2))), source = source)
}

# Clip X into the public support box.  The companion of dp.public.clip.y, and the
# step that makes dp.public.B.x a bound rather than a hope.  On the realised
# samples it is the identity map (income tops out near 3778 against the box's
# 4000), but the universe it defines is what the sensitivity is computed over.
dp.public.clip.x <- function(d) {
    l.spec <- dp.public.xspec()
    for(nm in names(l.spec)) {
        s <- l.spec[[nm]]
        if(!identical(s[["type"]], "continuous")) next
        if(is.null(d[[nm]])) next
        v.r <- range(s[["breaks"]])
        d[[nm]] <- pmin(pmax(d[[nm]], v.r[[1L]]), v.r[[2L]])
    }
    d
}

# +------------------------------------------------------------------------+
# |  The public whitening map and the statistic sensitivities              |
# +------------------------------------------------------------------------+
# The whitened sufficient-statistic release (66-dp-suffstat.r) releases
# vech(X'X), X'y, and y'y on the scale w = W x with W = Sigma^{-1/2} and
# Sigma = E[x x'] under the published DGP (64-dp-design-constants.r, N = 1e6
# at the public seed).  W is a fixed public affine map, so applying it costs
# no budget, and it conditions X'X near the identity: the survey measured
# condition 171.3 on the unit-box scale against 1.106 whitened, and the
# coefficient error moved from 484.5% to 29.8% at eps = 10 with X private
# (65-dp-suffstat-survey.r).
#
# B_x on the whitened scale is the maximum of ||W x|| over the public support
# box.  ||W x||^2 is a convex quadratic in the continuous coordinates for
# each fixed dummy pattern, so its box maximum sits at a vertex, and the
# categoricals contribute finitely many patterns.  The enumeration below is
# therefore exact, and it reads only xspec and Sigma, not the sample.
.dp.box.vertices <- function() {
    l.spec <- dp.public.xspec()
    v.age    <- range(l.spec[["age"]][["breaks"]])
    v.income <- range(l.spec[["income"]][["breaks"]])
    l.g <- expand.grid(age = v.age, female = c(0, 1),
                       race = c("1", "2", "3"), income = v.income)
    m.V <- cbind(1, l.g[["age"]], l.g[["female"]],
                 as.integer(l.g[["race"]] == "2"),
                 as.integer(l.g[["race"]] == "3"),
                 l.g[["income"]])
    colnames(m.V) <- c("(Intercept)", "age", "female1", "race2", "race3",
                       "income")
    m.V
}

dp.public.whitening <- function(sc) {
    l.k <- .dp.const(sc)
    if(is.null(l.k[["Sigma"]]))
        stop("dp-public-bounds: no Sigma for scenario ", sc, ". Re-run ",
             "64-dp-design-constants.r, which publishes it.")
    m.Sigma <- l.k[["Sigma"]]
    l.e <- eigen(m.Sigma, symmetric = TRUE)
    m.W <- l.e[["vectors"]] %*% diag(1 / sqrt(l.e[["values"]])) %*%
        t(l.e[["vectors"]])
    dimnames(m.W) <- dimnames(m.Sigma)
    s.Bx <- sqrt(max(rowSums((.dp.box.vertices() %*% t(m.W))^2)))
    list(W = m.W, B.x = s.Bx,
         source = sprintf(paste0("W = Sigma^{-1/2} with Sigma = E[xx'] over ",
                                 "the clipped public box; B_x = max ||W x|| ",
                                 "over the box vertices (exact: the quadratic ",
                                 "is convex per dummy pattern); %s"),
                          l.k[["source"]]))
}

# The public y standardisation the release uses alongside the whitening.  Any
# public affine map of y is valid (the refit maps back exactly).  The DGP's own
# published marginal parameters are the natural one.  B_y is the standardised
# public clip, so it exists exactly where dp.public.B.y does and the refusal
# on the unbounded scenarios is inherited rather than re-implemented.
dp.public.y.std <- function(sc) {
    l.y <- dp.public.B.y(sc)                 # refuses scenarios without a bound
    l.par <- if(sc %in% c(1L, 2L, 10L)) {
        list(centre = 50, scale = 15, s = "DGP marginal N(50, 15)")
    } else if(sc == 9L) {
        list(centre = 50, scale = 20,
             s = "DGP latent marginal N(50, 20), structural range [0, 100]")
    } else {
        stop("dp-public-bounds: no published y standardisation for scenario ",
             sc)
    }
    v.r <- `if`(is.null(l.y[["y.clip"]]), c(0, 100), l.y[["y.clip"]])
    s.By <- max(abs((v.r - l.par[["centre"]]) / l.par[["scale"]]))
    list(centre = l.par[["centre"]], scale = l.par[["scale"]], B.y = s.By,
         y.clip = v.r,
         source = sprintf("%s; standardised public range [%g, %g] gives B_y = %g",
                          l.par[["s"]],
                          (v.r[[1L]] - l.par[["centre"]]) / l.par[["scale"]],
                          (v.r[[2L]] - l.par[["centre"]]) / l.par[["scale"]],
                          s.By))
}

# The three statistic sensitivities, one-row replacement, on the whitened
# scale.  Every constant is a functional of (W, the box, the y clip), all
# public, so the trio exists even when X is sensitive, which is what the
# whole route is for (Delta_beta-style bounds read the realised X and are
# unavailable there, per dp.public.B.x.design's warning).
#
#   vech(X'X)  Delta = sqrt(2) B_x^2.  Proven upper bound: for the move
#              A = u u' - v v', ||vech(A)||^2 = (||A||_F^2 + sum_i A_ii^2)/2,
#              ||A||_F^2 = |u|^4 + |v|^4 - 2(u'v)^2 <= 2 B_x^4, and
#              sum_i A_ii^2 <= |u|^4 + |v|^4 <= 2 B_x^4 since
#              (a - b)^2 <= a^2 + b^2 for a, b >= 0 and sum u_i^4 <= |u|^4.
#              Equality needs u, v each one-hot at norm B_x with disjoint
#              support, and the whitened box contains no such pair, so the
#              bound is strict here (the A3 pattern: sound, over-noising in
#              the safe direction.  The measured worst pair is the verify
#              suite's diagnostic, not the calibration).
#   X'y        Delta = 2 B_x B_y, exactly as 20-dp-precheck.r.
#   y'y        Delta = B_y^2, attained at y_i = 0 -> y_i' = B_y.
dp.public.suffstat.sens <- function(l.white, l.ystd) {
    s.Bx <- l.white[["B.x"]]; s.By <- l.ystd[["B.y"]]
    list(Delta.vech = sqrt(2) * s.Bx^2,
         Delta.xty  = 2 * s.Bx * s.By,
         Delta.yy   = s.By^2,
         source = sprintf(paste0("sqrt(2) B_x^2 (strict upper bound on the ",
                                 "whitened box), 2 B_x B_y, and B_y^2 ",
                                 "(attained), at B_x = %.6g, B_y = %.6g"),
                          s.Bx, s.By))
}

# +------------------------------------------------------------------------+
# |  Provenance, printed                                                   |
# +------------------------------------------------------------------------+
# Call this once at the top of every driver.  A bound whose provenance is not
# printed is a bound nobody audits.
dp.bounds.announce <- function(...) {
    l.b <- list(...)
    cat("PUBLIC BOUNDS (none of these is a statistic of the analysis sample):\n")
    invisible(Map(function(s.nm, l.x)
        cat(sprintf("  %-10s %s\n", s.nm, l.x[["source"]])),
        names(l.b), l.b))
    cat("\n")
}
