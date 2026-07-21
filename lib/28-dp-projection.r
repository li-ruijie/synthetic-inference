#!/usr/bin/env r
# +------------------------------------------------------------+
# |  DP-Projection synthetic data generator                    |
# +------------------------------------------------------------+
#
# Method: output perturbation of beta_hat + sufficient statistics
# projection.  Two stages:
#
#   Stage 1:  beta_DP = beta_hat + eta,  eta ~ N(0, sigma^2 I_p),
#             sigma calibrated to a sensitivity of beta_hat.
#   Stage 2 (post-processing):  generate synthetic data y* such that
#             OLS(X_syn, y*) = beta_DP exactly, via projection.
#
# Caveat: the sensitivity here is local (DFBETA-based, see below).  Local
# sensitivity is data-dependent and does not provide (eps, delta)-DP, so
# gen.syn.dp.projected carries no formal DP claim.  It is retained as the projection engine and for the
# non-private inferential-fidelity experiments.  The sound mu-GDP routes
# are gen.syn.dp.ols.public (public X, global sensitivity, in the
# companion 33-dp-projection-public-x.r) and gen.syn.dp.projected.subagg
# (subsample-and-aggregate, worst-case per-coordinate sensitivity).  If
# beta_DP came from one of those sound mechanisms, the projected data
# would inherit its GDP by post-processing.
#
# Virtual-data trick (canonical-link GLMs):
# -----------------------------------------
# The canonical score equation is X'(y - mu(beta)) = 0.  At the
# MLE, X'y = X'mu(beta_hat).  Define
#     y_virtual = y + mu(beta_DP) - mu(beta_hat)
# Then X'y_virtual = X'mu(beta_DP), so the MLE of a model fit to
# (X, y_virtual) is beta_DP.  Passing (X, y_virtual) through
# gen.syn.projected() produces synthetic data that reproduces
# beta_DP exactly (up to Poisson rounding error).
#
# For identity link (Gaussian OLS): mu(beta) = X beta, so
#   y_virtual = y + X (beta_DP - beta_hat) = y + X eta
# and residuals e = y - X beta_hat are unchanged, preserving
# the dispersion sigma^2_OLS.
#
# Sensitivity:
# ------------
# Local sensitivity of beta_hat at dataset D:
#     LS_beta(D) = max_i ||(X'X)^{-1} x_i ( y_i - x_i' beta_hat )||_2
#                  (leave-one-out DFBETA formula, Cook 1977)
# This local sensitivity is data-dependent, so it is a disclosure metric
# only and does not yield (eps, delta)-DP (a public bound on |y| does not
# fix this).  The sound routes replace it with a global sensitivity
# (gen.syn.dp.ols.public) or a worst-case per-coordinate sensitivity over
# data-independent blocks (subsample-and-aggregate).
#
# Dependencies: 11-projection-generator.r (project.ols, gen.syn.projected)
#               gen-syn-family.r (gen.syn.family, the copula generator)

# +------------------------------------------------------------+
# |  Setup                                                     |
# +------------------------------------------------------------+
# standalone: source only if not already loaded
if(!exists("gen.syn.projected"))
    source("lib/11-projection-generator.r")
if(!exists("gen.syn.family"))
    source("lib/gen-syn-family.r")

# +------------------------------------------------------------+
# |  beta_hat sensitivity                                      |
# +------------------------------------------------------------+
# L2 local sensitivity of beta_hat under one-observation
# replacement, estimated via leave-one-out influence (DFBETA).
#
# For observation i replaced: change in beta_hat ~ DFBETA_i (removal)
# + DFBETA_{new} (addition).  Upper bound by 2 * max_i ||DFBETA_i||_2.
#
# LS(D) = 2 * max_i || beta_hat(D) - beta_hat(D_{-i}) ||_2
#
# This is local sensitivity.  For formal DP use smooth sensitivity
# (Nissim-Raskhodnikova-Smith 2007) which propagates LS through
# neighboring datasets.  For this first implementation, we return
# LS with the caveat that the DP guarantee is empirical rather
# than worst-case formal: valid when the underlying data generating
# process produces datasets with bounded influence similar to D.
#
# Comparison to worst-case global sensitivity with bounded-y:
#   Delta_global = 2 * y.max * max_i ||(X'X)^{-1} x_i||_2
# The worst-case bound is typically 10-100x larger than LS because
# it allows any y in [-y.max, y.max], whereas LS uses actual residuals
# e_i = y_i - x_i' beta_hat which are O(residual_sd) << y.max.
#
# Arguments:
#   mo       fitted lm or glm object.
# Returns: list with
#   Delta         local sensitivity (DFBETA, a disclosure metric only,
#                 not a DP-calibration sensitivity).
#   Delta.global  worst-case global sensitivity (for comparison).
sens.beta.ols <- function(mo, y.range = NULL) {
    m.dfbeta <- dfbeta(mo)                              # n x p
    v.dfb.norms <- sqrt(rowSums(m.dfbeta^2))
    Delta.local <- 2 * max(v.dfb.norms)

    Delta.global <- NA_real_
    if(!is.null(y.range)) {
        m.X <- model.matrix(mo)
        m.XtXinv <- tryCatch(solve(crossprod(m.X)),
                              error = function(e) NULL)
        if(!is.null(m.XtXinv)) {
            m.influence <- m.X %*% m.XtXinv
            v.row.norms <- sqrt(rowSums(m.influence^2))
            Delta.global <- 2 * max(abs(y.range)) * max(v.row.norms)
        }
    }
    list(Delta = Delta.local,
         Delta.local = Delta.local,
         Delta.global = Delta.global)
}

# +------------------------------------------------------------+
# |  Smooth sensitivity (Nissim-Raskhodnikova-Smith 2007)      |
# +------------------------------------------------------------+
# For a function f with local sensitivity LS(D), the beta-smooth
# sensitivity is S^beta(D) = sup_{D'} exp(-beta * d(D, D')) * LS(D')
# where d(,) is Hamming distance on datasets.  For Gaussian mechanism
# with (eps, delta)-DP: beta = eps / (2 * ln(2/delta)), and sigma = S^beta
# * sqrt(2 * ln(2/delta)) / eps.
#
# For OLS with bounded ||x_i||_2 <= B_x and |y_i| <= B_y publicly, the
# local sensitivity at k-neighbours is bounded by
#   LS(D'_k) <= LS(D) + k * dstep
# where dstep is a per-step LS increase from adding/replacing one
# extreme observation (x_ext, y_ext) with ||x_ext|| = B_x, |y_ext|= B_y:
#   dstep = 2 * B_x * B_y / lambda_min(X'X - x_ext x_ext^T)
# For practical use, approximate dstep via current design:
#   dstep ~ 2 * B_x * B_y / lambda_min(X'X)
# Fragility (documented, not hidden): this is an approximation, not an
# identity, and it errs in the unsafe direction.  Subtracting the
# rank-one PSD term x_ext x_ext^T cannot raise the smallest eigenvalue,
# so lambda_min(X'X - x_ext x_ext^T) <= lambda_min(X'X) and the
# approximated dstep sits at or below the true one.  Where k* > 0 the
# smooth sensitivity is increasing in dstep, since its derivative
# carries the factor (1 - beta * LS / dstep) and k* > 0 is exactly
# beta * LS / dstep < 1, so the sigma this returns sits at or below
# what the bound requires.  The gap is small for well-conditioned X
# with n >> p, which is the regime this route is used in, but it is a
# gap and not an equality.  Prefer the subsample-and-aggregate route
# (D2) for any formal claim.
#
# Closed-form maximisation of
#   g(k) = (LS(D) + k * dstep) * exp(-beta * k)
# yields k* = (1/(beta) - LS/dstep).  If k* > 0,
#   S^beta = (dstep / (beta * e)) * exp(beta * LS / dstep)
# Otherwise S^beta = LS(D) (the max is at k=0).
#
# Arguments:
#   mo        fitted lm/glm
#   eps       privacy budget epsilon > 0
#   delta     privacy budget delta in (0, 1)
#   B.x       public bound on ||x_i||_2
#   B.y       public bound on |y_i|
# Returns: list with S.smooth, S.local, beta, ratio (S.smooth/S.local)
smooth.sens.beta.ols <- function(mo, eps, delta, B.x = NULL, B.y = NULL) {
    m.X <- model.matrix(mo)
    y <- model.response(model.frame(mo))
    # local sensitivity (DFBETA-based upper bound)
    m.dfbeta <- dfbeta(mo)
    LS <- 2 * max(sqrt(rowSums(m.dfbeta^2)))
    # public bounds: default to observed (not formally DP but
    # caller can override with true public bounds).
    if(is.null(B.x)) B.x <- sqrt(max(rowSums(m.X^2)))
    if(is.null(B.y)) B.y <- max(abs(y))
    # per-step increase in LS from adding one extreme observation
    # dstep = 2 * B_x * B_y / lambda_min(X'X)
    lambda.min <- tryCatch(
        min(eigen(crossprod(m.X), symmetric = TRUE,
                  only.values = TRUE)[["values"]]),
        error = function(e) NA_real_)
    if(is.na(lambda.min) || lambda.min <= 0)
        return(list(S.smooth = NA_real_, S.local = LS, beta = NA_real_,
                    ratio = NA_real_, dstep = NA_real_))
    dstep <- 2 * B.x * B.y / lambda.min
    # beta parameter for Gaussian mechanism smooth sens
    beta.p <- eps / (2 * log(2 / delta))
    # closed-form maximiser
    k.star <- 1/beta.p - LS/dstep
    S.smooth <- `if`(k.star > 0,
        (dstep / (beta.p * exp(1))) * exp(beta.p * LS / dstep),
        LS)
    list(S.smooth = S.smooth, S.local = LS, beta = beta.p,
         ratio = S.smooth / LS, dstep = dstep, k.star = max(0, k.star))
}

# +------------------------------------------------------------+
# |  Smooth sensitivity for log-link GLM MLE (DP-D1, "try")    |
# +------------------------------------------------------------+
# Analytic smooth sensitivity for a log-link GLM coefficient, the
# Fisher-information analog of smooth.sens.beta.ols.  This is the smaller
# but fragile route (design-doc D1): the smooth-sensitivity growth bound
# is derived by treating the GLM one-step DFBETA like the OLS influence,
# with the Fisher information X'WX in place of X'X and a public bound on
# the score-residual magnitude g_i = (dmu/deta)_i / V(mu_i) * (y_i - mu_i)
# in place of B_y.  The subsample-and-aggregate route (D2) is the primary,
# assumption-free worst-case method.  This is carried for comparison.
#
# Fragility (documented, not hidden): the IRLS weights W depend on beta,
# so X'WX and g_i shift under one-row changes in a way this bound does not
# fully capture.  The result is a heuristic smooth sensitivity, valid when
# W and lambda_min(X'WX) are stable across neighbours (large n, well-
# conditioned design), not a proof-grade worst-case bound.  Prefer D2 for
# the formal claim.  Use this to quantify how much smaller an analytic
# route could be when its assumptions hold.
#
# Arguments:
#   mo         fitted log-link glm.
#   eps, delta smooth-sensitivity DP budget.
#   B.x        public bound on ||x_i||_2 (default: observed max).
#   B.r        public bound on |g_i| (default: observed max).
# Returns: list(S.smooth, S.local, dstep, ratio, lambda.min, B.x, B.r).
smooth.sens.beta.loglink <- function(mo, eps, delta, B.x = NULL, B.r = NULL) {
    m.X   <- model.matrix(mo)
    fam   <- family(mo)
    v.mu  <- fitted(mo)
    v.eta <- fam[["linkfun"]](v.mu)
    v.dmu.deta <- fam[["mu.eta"]](v.eta)
    v.V   <- fam[["variance"]](v.mu)
    v.w   <- v.dmu.deta^2 / v.V                         # IRLS weights
    m.XtWX <- crossprod(sweep(m.X, 1, sqrt(v.w), "*"))  # Fisher information

    LS <- 2 * max(sqrt(rowSums(dfbeta(mo)^2)))          # local sensitivity

    if(is.null(B.x)) B.x <- sqrt(max(rowSums(m.X^2)))
    v.y  <- model.response(model.frame(mo))
    v.g  <- (v.dmu.deta / v.V) * (v.y - v.mu)           # score-residual
    if(is.null(B.r)) B.r <- max(abs(v.g))

    lambda.min <- tryCatch(
        min(eigen(m.XtWX, symmetric = TRUE, only.values = TRUE)[["values"]]),
        error = function(e) NA_real_)
    if(is.na(lambda.min) || lambda.min <= 0)
        return(list(S.smooth = NA_real_, S.local = LS, dstep = NA_real_,
                    ratio = NA_real_, lambda.min = lambda.min,
                    B.x = B.x, B.r = B.r))

    # per-neighbour LS growth from adding one extreme observation.
    dstep  <- 2 * B.x * B.r / lambda.min
    beta.p <- eps / (2 * log(2 / delta))
    k.star <- 1 / beta.p - LS / dstep
    S.smooth <- `if`(k.star > 0,
        (dstep / (beta.p * exp(1))) * exp(beta.p * LS / dstep), LS)
    list(S.smooth = S.smooth, S.local = LS, dstep = dstep,
         ratio = S.smooth / LS, lambda.min = lambda.min,
         B.x = B.x, B.r = B.r)
}

# +------------------------------------------------------------+
# |  D1 made viable: public standardisation of the design     |
# +------------------------------------------------------------+
# smooth.sens.beta.loglink is unusable on a raw design with a large-scale
# covariate: B.x ~ max||x_i|| is dominated by that covariate and
# lambda_min(X'WX) is near zero, so dstep = 2 B.x B.r / lambda_min inflates the
# sensitivity by four to five orders of magnitude (income in scenarios 3/6).
#
# This wrapper makes the analytic route viable by standardising the design
# with a public per-coefficient scale c_j: X_std = X diag(c)^{-1}.  With c_j
# the public magnitude of covariate j (intercept and factor dummies 1, a
# continuous covariate its public scale), X_std has unit-scale columns, so
# X_std' W X_std is well-conditioned (lambda_min = O(1)) and the growth bound
# is finite and attained.
#
# DP.  The reparametrisation beta_std = diag(c) beta with c public is a public
# linear map, so the smooth-sensitivity Gaussian mechanism on beta_std is
# (eps, delta)-DP and beta_DP = diag(c)^{-1} beta_std_DP is post-processing.
# c is the public input the analytic route requires, exactly as B.x, B.y are
# for the OLS bound and clip.lo/clip.hi are for subsample-and-aggregate.
# Standardisation removes the conditioning failure but not the underlying
# heuristic: the IRLS weights W still shift under one-row changes, so this is
# an analytic smooth-sensitivity route whose regularity (stable, well-
# conditioned Fisher information) is induced by the public standardisation,
# not the assumption-free worst-case guarantee of subsample-and-aggregate.
# Prefer D2 for the formal Proposition.  D1-std is the smaller analytic
# alternative when a public covariate scaling is available.
#
# Arguments:
#   mo         fitted log-link glm.
#   eps, delta smooth-sensitivity DP budget.
#   x.scale    public named per-coefficient scale c (aligned to coef names),
#              positive.  Required.
#   B.x, B.r   public bounds on ||x_std,i||_2 and |g_i| (default: observed
#              max on the standardised design / raw score-residual).
# Returns: list(S.smooth.std, S.local.std, dstep, ratio, lambda.min,
#               sigma.native, x.scale, B.x, B.r).  sigma.native is the
#               per-coefficient native-scale Gaussian-mechanism noise sd.
smooth.sens.beta.loglink.std <- function(mo, eps, delta, x.scale,
                                         B.x = NULL, B.r = NULL) {
    v.beta <- coef(mo)
    m.X    <- model.matrix(mo)
    v.c    <- x.scale[names(v.beta)]
    if(any(is.na(v.c)) || any(v.c <= 0))
        stop("x.scale must give a positive public scale for every coefficient")

    # standardise the design by the public scales: X_std = X diag(c)^{-1}.
    m.X.std <- sweep(m.X, 2, v.c, "/")

    fam   <- family(mo)
    v.mu  <- fitted(mo)
    v.eta <- fam[["linkfun"]](v.mu)
    v.dmu.deta <- fam[["mu.eta"]](v.eta)
    v.V   <- fam[["variance"]](v.mu)
    v.w   <- v.dmu.deta^2 / v.V                              # IRLS weights (fit-invariant)
    m.XtWX <- crossprod(sweep(m.X.std, 1, sqrt(v.w), "*"))   # standardised Fisher info

    # local sensitivity of the standardised coefficient beta_std = c * beta.
    m.dfbeta.std <- sweep(dfbeta(mo), 2, v.c, "*")
    S.local <- 2 * max(sqrt(rowSums(m.dfbeta.std^2)))

    if(is.null(B.x)) B.x <- sqrt(max(rowSums(m.X.std^2)))
    v.y  <- model.response(model.frame(mo))
    v.g  <- (v.dmu.deta / v.V) * (v.y - v.mu)                # score-residual (scale-free)
    if(is.null(B.r)) B.r <- max(abs(v.g))

    lambda.min <- tryCatch(
        min(eigen(m.XtWX, symmetric = TRUE, only.values = TRUE)[["values"]]),
        error = function(e) NA_real_)
    if(is.na(lambda.min) || lambda.min <= 0)
        return(list(S.smooth.std = NA_real_, S.local.std = S.local,
                    dstep = NA_real_, ratio = NA_real_, lambda.min = lambda.min,
                    sigma.native = NA_real_, x.scale = v.c, B.x = B.x, B.r = B.r))

    dstep  <- 2 * B.x * B.r / lambda.min
    beta.p <- eps / (2 * log(2 / delta))
    k.star <- 1 / beta.p - S.local / dstep
    S.smooth <- `if`(k.star > 0,
        (dstep / (beta.p * exp(1))) * exp(beta.p * S.local / dstep), S.local)

    # native per-coefficient noise sd: isotropic N(0, (S.smooth * g)^2) on
    # beta_std, un-standardised by c_j (post-processing).
    s.gm <- dp.gm.multiplier(eps, delta)
    v.sigma.native <- `names<-`((S.smooth * s.gm) / v.c, names(v.beta))

    list(S.smooth.std = S.smooth, S.local.std = S.local, dstep = dstep,
         ratio = S.smooth / S.local, lambda.min = lambda.min,
         sigma.native = v.sigma.native, x.scale = v.c, B.x = B.x, B.r = B.r)
}

# +------------------------------------------------------------+
# |  Gaussian-mechanism noise multiplier (per unit sensitivity)|
# +------------------------------------------------------------+
# gdp.delta() and dp.gm.multiplier() are defined in dp-gaussian-mechanism.r,
# which is their single definition.  They are not defined here.
source("lib/dp-gaussian-mechanism.r")

# +------------------------------------------------------------+
# |  mu-Gaussian DP composition and calibration                |
# +------------------------------------------------------------+
# gdp.delta(eps, mu) above is the (eps, delta) profile of a mu-GDP Gaussian
# mechanism.  These wrap it for composition and calibration, stated purely in
# mu-GDP for all eps > 0 (Dong, Roth & Su 2022).

# Compose mu_1, ..., mu_k GDP mechanisms: total is sqrt(sum mu_i^2)-GDP.
gdp.compose <- function(v.mu) sqrt(sum(v.mu^2))

# (eps, delta) profile at multiplier mu, with argument order (mu, eps).
gdp.to.delta <- function(mu, eps) gdp.delta(eps, mu)

# Largest mu a mechanism may have at target (eps, delta).  Exact for all eps.
gdp.mu.from.eps.delta <- function(eps, delta) {
    stopifnot(delta > 0, delta < 1, eps > 0)
    uniroot(function(m) gdp.delta(eps, m) - delta,
            interval = c(1e-8, 1e3), tol = 1e-12)[["root"]]
}

# Gaussian sd to reach (eps, delta) at L2 sensitivity Delta: sigma = Delta / mu*.
gdp.sigma <- function(s.Delta, eps, delta)
    s.Delta / gdp.mu.from.eps.delta(eps, delta)

# +------------------------------------------------------------+
# |  Route 1 global sensitivities (public X, |y_i| <= B.y)     |
# +------------------------------------------------------------+
# OLS beta_hat = (X'X)^{-1} X'y.  One-row y replacement (range 2 B.y) moves
# beta by (X'X)^{-1} x_i * dy, so Delta = 2 B.y max_i ||(X'X)^{-1} x_i||_2.
# X must have full column rank: without it (X'X)^{-1} does not exist and Delta
# is undefined, which is an undefined noise scale rather than a loose one, so
# solve() is left to error.  This is a hypothesis of the public-covariate exact
# release proposition (the paper's supplementary privacy section).
gs.beta.ols.public <- function(m.X, B.y) {
    m.inv <- solve(crossprod(m.X))
    2 * B.y * max(sqrt(rowSums((m.X %*% m.inv)^2)))
}

# sigma_hat^2 = ||y - X beta_hat||^2 / (n - p) = y'(I - H)y / (n - p) with
# H = X (X'X)^{-1} X'.  Replacing row i (y' = y + d e_i, |d| <= 2 B.y) gives
# the exact identity RSS(y') - RSS(y) = 2 d r_i + d^2 (1 - h_ii), r_i the
# i-th residual.  Over the bounded universe |y_j| <= B.y the worst case is
# |Delta RSS| = 4 B.y^2 max_i max{1 - h_ii, sum_{j != i} |h_ij|}, attained at
# an adversarial sign pattern, so the design constant K_X multiplies the
# naive 4 B.y^2 bound.  K_X depends only on the public X (no privacy cost).
# Omitting K_X (3.3 to 6.1 on the paper's scenario designs) under-noises the
# variance channel.
gs.sigma2 <- function(m.X, B.y) {
    n <- nrow(m.X); p <- ncol(m.X)
    m.H <- m.X %*% solve(crossprod(m.X)) %*% t(m.X)
    v.h <- diag(m.H)
    v.s <- rowSums(abs(m.H)) - abs(v.h)
    4 * B.y^2 * max(pmax(1 - v.h, v.s)) / (n - p)
}

# +------------------------------------------------------------+
# |  DP noise: Gaussian mechanism                              |
# +------------------------------------------------------------+
# sigma = Delta * dp.gm.multiplier(eps, delta): the classical
# sqrt(2 log(1.25/delta))/eps for eps < 1, and the exact analytic /
# Gaussian-DP calibration for eps >= 1 (see dp.gm.multiplier above).
dp.gaussian <- function(v.x, eps, delta, s.Delta, seed = NA) {
    if(!is.na(seed)) set.seed(seed)
    s.sigma <- s.Delta * dp.gm.multiplier(eps, delta)
    list(v.x.dp = v.x + rnorm(length(v.x), 0, s.sigma),
         sigma  = s.sigma,
         eta    = rnorm(0))                             # drained by prev call
}

# dp.gaussian.noise: return noise vector instead of already-added.
dp.gaussian.noise <- function(s.p, eps, delta, s.Delta, seed = NA) {
    if(!is.na(seed)) set.seed(seed)
    s.sigma <- s.Delta * dp.gm.multiplier(eps, delta)
    list(v.eta = rnorm(s.p, 0, s.sigma),
         sigma = s.sigma)
}

# +------------------------------------------------------------+
# |  dp.perturb.beta: shared DP coefficient perturbation       |
# +------------------------------------------------------------+
# Factored-out coefficient-DP step shared by the DP generators.
# Given a fitted lm/glm, returns beta_DP = un-standardise(beta_tilde
# + eta), where beta_tilde = s * beta_hat is the column-standardised
# coefficient and eta is Gaussian mechanism noise calibrated to the
# (standardised) DFBETA local sensitivity.  This is exactly the step
# inlined in gen.syn.dp.projected (Step 2-3 below) and duplicated in
# 53-dp-logistic.r.  The full DP-ZI zero-model (gen.syn.dp.projected.zi)
# reuses it here instead of re-inlining.
#
# The standardisation, DFBETA sensitivity, seed offset (seed + 1L), and
# norm clipping match gen.syn.dp.projected coefficient-for-coefficient,
# so the helper is a faithful extraction (verified by test-dp-perturb-
# beta.R, which checks beta_DP agreement to machine precision).  The
# existing generator is left untouched: this is additive, not a refactor.
#
# Arguments:
#   mo            fitted lm or glm.
#   eps, delta    DP budget for this coefficient release.
#   seed          RNG seed.  Noise drawn under set.seed(seed + 1L).
#   standardise   column-standardise X before adding noise.
#   clip.sigma    truncate ||eta||_2 at clip.sigma * sqrt(p) * sigma
#                 (Lipschitz post-processing, so DP is preserved).
#   sens.method   "local" (DFBETA, empirical DP) or "smooth" (NRS 2007,
#                 gaussian-identity lm only.  Other fits fall back to
#                 local with a caveat, matching the inline code).
#   B.x, B.y      public bounds for smooth sensitivity.
# Returns: list(beta.hat, beta.dp, sigma, Delta, Delta.smooth,
#               col.scale, eta, sens.method).
dp.perturb.beta <- function(mo, eps, delta, seed = 6489,
                            standardise = TRUE, clip.sigma = Inf,
                            sens.method = c("local", "smooth"),
                            B.x = NULL, B.y = NULL) {
    sens.method <- match.arg(sens.method)
    v.beta.hat <- coef(mo)
    m.X <- model.matrix(mo)

    # column-standardise: s_j = sd(X[,j]) for non-constant columns.
    v.col.sd <- apply(m.X, 2, stats::sd)
    v.s <- ifelse(v.col.sd > 1e-8, v.col.sd, 1)
    names(v.s) <- colnames(m.X)
    v.s <- v.s[names(v.beta.hat)]                       # align order
    if(!standardise) v.s <- rep(1, length(v.beta.hat))
    v.beta.tilde <- v.beta.hat * v.s

    # standardised DFBETA local sensitivity.
    m.dfbeta <- dfbeta(mo)
    m.dfbeta.tilde <- sweep(m.dfbeta, 2, v.s, "*")
    s.LS <- 2 * max(sqrt(rowSums(m.dfbeta.tilde^2)))

    # smooth sensitivity: gaussian-identity (lm) only, else local
    # (log-link glms need non-linear treatment, deferred to the
    # worst-case route.  Mirrors the inline fallback).
    s.S.smooth <- NA_real_
    if(sens.method == "smooth" && inherits(mo, "lm") &&
       !inherits(mo, "glm")) {
        m.X.std <- sweep(m.X, 2, v.s, "/")              # standardised
        v.y <- model.response(model.frame(mo))
        d.std.refit <- data.frame(m.X.std[, -1, drop = FALSE])
        d.std.refit[[".dp.y"]] <- v.y
        names(d.std.refit) <- c(colnames(m.X.std)[-1], ".dp.y")
        fml.std <- reformulate(colnames(m.X.std)[-1], response = ".dp.y")
        mo.std <- lm(fml.std, d.std.refit)
        l.smooth <- smooth.sens.beta.ols(mo.std, eps, delta,
                                          B.x = B.x, B.y = B.y)
        s.S.smooth <- 2 * l.smooth[["S.smooth"]]        # scale factor
        if(is.finite(s.S.smooth)) s.LS <- s.S.smooth
    }
    s.Delta <- s.LS

    # Gaussian mechanism noise on the standardised scale, optional clip.
    set.seed(seed + 1L)
    l.dp <- dp.gaussian.noise(length(v.beta.hat), eps, delta, s.Delta)
    v.eta.std <- l.dp[["v.eta"]]
    s.p <- length(v.eta.std)
    s.eta.norm <- sqrt(sum(v.eta.std^2))
    s.eta.max <- clip.sigma * sqrt(s.p) * l.dp[["sigma"]]
    if(is.finite(s.eta.max) && s.eta.norm > s.eta.max)
        v.eta.std <- v.eta.std * (s.eta.max / s.eta.norm)
    v.beta.tilde.dp <- v.beta.tilde + v.eta.std
    v.beta.dp <- v.beta.tilde.dp / v.s                  # un-standardise
    names(v.beta.dp) <- names(v.beta.hat)

    list(beta.hat     = v.beta.hat,
         beta.dp      = v.beta.dp,
         sigma        = l.dp[["sigma"]],
         Delta        = s.Delta,
         Delta.smooth = s.S.smooth,
         col.scale    = v.s,
         eta          = v.beta.dp - v.beta.hat,
         sens.method  = sens.method)
}

# +------------------------------------------------------------+
# |  Dispersion provenance: the origin of sigma^2              |
# +------------------------------------------------------------+
# A GLM's likelihood has two parameters: the coefficient and the dispersion.  DP
# must cover both.  That means saying explicitly where the dispersion came from,
# because the two defensible answers carry different privacy claims:
#
#   Public   The dispersion is already public.  Either it is external knowledge
#            (a design constant, prior literature), or the data owner has already
#            Published it -- journal articles routinely report R^2 or a residual
#            SD alongside the coefficient table.  Post-processing permits any
#            function of (DP output, public information), so re-using it costs no
#            budget and discloses nothing new.  Where it was already published by
#            the owner, the claim is mu-GDP conditional on that publication: the
#            synthetic data adds nothing the table did not already give away.
#            This mode spends no budget and leaves beta_DP as noisy as before.
#
#   Private  Nothing has been published, so the dispersion is spent from the
#            budget through the subagg channel (disp = TRUE).  This costs: the L2
#            sensitivity becomes sqrt(p + 1) rather than sqrt(p), so beta_DP is
#            about 8% noisier at p = 6, and the standard errors inherit noise.
#
# `source` is required and is printed.  That is the point of the constructor.  A
# dispersion silently lifted from the sample -- which is what shipped -- now has
# to be written as dp.dispersion.public(var(resid(fit)), source = "the sample"),
# which is self-indicting on the page.
dp.dispersion.public <- function(sigma2, source) {
    if(missing(source) || !is.character(source) || !nzchar(source))
        stop("source is REQUIRED: name where this PUBLIC dispersion comes from ",
             "(external constant, prior literature, or the owner's own published ",
             "table).  If it came from the sample, it is not public.")
    if(!is.numeric(sigma2) || length(sigma2) != 1L || !is.finite(sigma2) || sigma2 <= 0)
        stop("sigma2 must be a single positive finite number")
    structure(list(mode = "public", sigma2 = sigma2, source = source),
              class = "dp.dispersion")
}

dp.dispersion.private <- function(lo, hi, source) {
    if(missing(source) || !is.character(source) || !nzchar(source))
        stop("source is REQUIRED: name where this PUBLIC BOX on the dispersion ",
             "comes from.  The box width sets the noise scale, so a box read off ",
             "the sample makes the noise scale a function of the private data.")
    if(!is.numeric(lo) || !is.numeric(hi) || length(lo) != 1L || length(hi) != 1L ||
       !is.finite(lo) || !is.finite(hi) || hi <= lo)
        stop("lo and hi must be single finite numbers with hi > lo")
    structure(list(mode = "private", lo = lo, hi = hi, source = source),
              class = "dp.dispersion")
}

print.dp.dispersion <- function(x, ...) {
    cat(sprintf("<dp.dispersion: %s>\n", x[["mode"]]))
    if(x[["mode"]] == "public")
        cat(sprintf("  sigma^2 = %g (no budget spent)\n", x[["sigma2"]]))
    else
        cat(sprintf("  public box [%g, %g] (spends budget: sqrt(p+1) sensitivity)\n",
                    x[["lo"]], x[["hi"]]))
    cat(sprintf("  source: %s\n", x[["source"]]))
    invisible(x)
}

# Families with a free dispersion parameter.  Poisson and binomial have none --
# their variance is determined by the mean -- so they need no dispersion channel
# and cannot leak one.
dp.has.dispersion <- function(family)
    family %in% c("gaussian", "gamma", "quasi-gamma")

# +------------------------------------------------------------+
# |  dp.rate.public / dp.subagg.blocks: the block-count rule   |
# +------------------------------------------------------------+
# m = floor(n_eff / (5 p)), the largest block count leaving five rows per
# released coefficient per block.  n is public under the replace-one
# neighbouring convention, p is the number of released coefficients, and rate
# (a dp.rate object) scales n to the expected positive count n * (1 - pi0) on
# a two-part positive channel, with pi0 the published zero rate.  Passing a
# realised rate (mean(y == 0)) would calibrate the mechanism to the data being
# protected, which is unsound, so the rate travels as a
# provenance object on the dp.dispersion.public pattern and
# dp.public.zero.rate() (dp-public-bounds.r) serves it for the scenario
# designs.
dp.rate.public <- function(pi0, source) {
    if(missing(source) || !is.character(source) || length(source) != 1L ||
       !nzchar(source))
        stop("source is REQUIRED: name where this PUBLIC zero rate comes from ",
             "(a published design constant, not the sample).  If it came ",
             "from the sample, it is not public.")
    if(!is.numeric(pi0) || length(pi0) != 1L || is.na(pi0) || pi0 < 0 || pi0 >= 1)
        stop("pi0 must be a single value in [0, 1)")
    structure(list(pi0 = pi0, source = source), class = "dp.rate")
}

print.dp.rate <- function(x, ...) {
    cat(sprintf("<dp.rate: pi0 = %g>\n  source: %s\n",
                x[["pi0"]], x[["source"]]))
    invisible(x)
}

dp.subagg.blocks <- function(n, p, rate = NULL) {
    if(!is.null(rate) && !inherits(rate, "dp.rate"))
        stop("rate must be built by dp.rate.public(), so that its provenance ",
             "is on the record")
    s.eff <- `if`(is.null(rate), n, n * (1 - rate[["pi0"]]))
    s.m <- floor(s.eff / (5 * p))
    if(s.m < 2L)
        stop(sprintf(paste0("dp.subagg.blocks: floor(n_eff / (5 p)) = %d < 2 ",
                            "(n_eff = %.1f, p = %d), too few rows to ",
                            "subsample-and-aggregate"), s.m, s.eff, p))
    as.integer(s.m)
}

# +------------------------------------------------------------+
# |  synth.response: build the release from the DP output only |
# +------------------------------------------------------------+
# This replaces inject.virtual.y in the DP generators, and the difference is the
# governs whether the post-processing clause applies at all.
#
# inject.virtual.y sets  y_virt = y + mu(beta_DP) - mu(beta_hat).  It retains the
# real y.  The projection then reproduces the real data's residuals in the
# release, in scale and in shape.  Measured, sigma_syn equals sigma_real to four
# significant figures.  Heavy-tailed input gave heavy-tailed output (kurtosis
# 7.0 -> 8.6), and skewed input gave skewed output (skew 1.8 -> 1.5), so the release
# did not depend on the data "only through beta_DP", which is the
# hypothesis the subsample-and-aggregate proposition needs for its
# post-processing clause.  The empirical audit (63-dp-audit-mugdp.r) scored the
# resulting release at AUC 1.000 against a mu-GDP bound of 0.566 at eps = 1: a
# perfect distinguisher, and flat in epsilon, because the leak did not go
# through the noise.
#
# synth.response instead draws the response from the model's own family, at the
# mean implied by beta_DP and at the declared dispersion.  The result is a
# function of (beta_DP, dispersion, X, seed) and nothing else, so post-processing
# genuinely applies.  The real y appears nowhere.
#
# Note this costs the release its distributional fidelity to the empirical
# residuals, and keeps its inferential fidelity, which is the trade the paper
# argues for throughout.
#
# Arguments:
#   v.beta.dp   the DP coefficient (the only data-dependent input, and it is DP).
#   s.sigma2    the dispersion, from a dp.dispersion object.  NA for poisson and
#               binomial, which have none.
# Returns: d with the response replaced.
synth.response <- function(d, y.nm, x.nms, family, link,
                            v.beta.dp, s.sigma2, seed) {
    m.X <- model.matrix(reformulate(x.nms), d)
    v.b <- `names<-`(rep(0, ncol(m.X)), colnames(m.X))
    v.common <- intersect(names(v.beta.dp), colnames(m.X))
    v.b[v.common] <- v.beta.dp[v.common]

    v.eta <- as.numeric(m.X %*% v.b)

    # Guard the log link against overflow.  At a small epsilon the public box is
    # wide and beta_DP is very noisy, so exp(X beta_DP) can reach Inf and the draw
    # below then yields NaN.  Clamping the synthetic linear predictor is a function
    # of beta_DP and X only, so it stays a post-processing and costs no budget.
    # It exists to keep the arithmetic finite and not to restore utility. A rep whose
    # coefficients are this extreme still fails the analyst's refit, and the drivers
    # report that as a dropped rep (their n.ok column).  Do not read this as
    # papering over the noise -- the noise is real and must be reported.
    v.eta <- pmin(pmax(v.eta, -500), 500)

    v.mu  <- switch(link, identity = v.eta, log = exp(v.eta),
                    logit = plogis(v.eta), stop("unsupported link: ", link))
    s.n <- nrow(m.X); s.p <- ncol(m.X)
    set.seed(seed + 11L)

    v.y <- switch(family,
        # Gaussian: draw, orthogonalise against the design, rescale to the target
        # dispersion.  The refit then recovers beta_DP and sigma2 exactly, so the
        # exact-fidelity property of the projection theorem survives the change
        # intact.
        gaussian = {
            v.e <- lm.fit(m.X, rnorm(s.n))[["residuals"]]
            s.ss <- sum(v.e^2)
            if(s.ss <= 0) stop("degenerate residual basis in synth.response")
            v.mu + v.e * sqrt(s.sigma2 * (s.n - s.p) / s.ss)
        },
        # Poisson and binomial carry no free dispersion: the family fixes the
        # variance given the mean, so the draw is complete on beta_DP alone.
        poisson  = rpois(s.n, pmax(v.mu, .Machine[["double.eps"]])),
        binomial = rbinom(s.n, 1L, pmin(pmax(v.mu, 1e-9), 1 - 1e-9)),
        # Gamma with dispersion phi: shape = 1/phi, rate = 1/(phi mu).
        gamma = ,
        `quasi-gamma` = {
            s.phi <- max(s.sigma2, 1e-8)
            rgamma(s.n, shape = 1 / s.phi,
                   rate = 1 / (s.phi * pmax(v.mu, .Machine[["double.eps"]])))
        },
        stop("unsupported family in synth.response: ", family))

    `[[<-`(d, y.nm, value = v.y)
}

# +------------------------------------------------------------+
# |  Virtual dataset trick: inject beta_DP                     |
# +------------------------------------------------------------+
# Given original d and target beta_DP, return d.virtual where
#   MLE(y ~ x, family = ..., d.virtual) == beta_DP.
# For canonical-link GLMs (identity-gaussian, log-poisson, log-gamma,
# logit-binomial), this uses the identity
#     X'y_virtual = X'mu(beta_DP)
# by setting y_virtual = y + mu(beta_DP) - mu(beta_hat).
#
# Warning, and it is the reason synth.response exists above.  y_virtual retains
# the real y.  Anything built from it therefore depends on the real data by more
# than beta_DP, and cannot claim the subsample-and-aggregate proposition's
# post-processing clause.  This function is appropriate for the analysis-scope
# generator (gen.syn.dp.projected), which makes no formal DP claim, but it must
# Not be used to build a release that does.
#
# Arguments:
#   d            original data frame.
#   y.nm         response column name.
#   v.beta.hat   original MLE coefficient vector.
#   v.beta.dp    DP-perturbed coefficient vector (same length + names).
#   family       one of "gaussian", "gamma", "poisson".
#   link         one of "identity", "log" (canonical only).
# Returns: modified data frame.
inject.virtual.y <- function(d, y.nm, v.beta.hat, v.beta.dp,
                              family, link) {
    v.x.nms <- setdiff(names(d), y.nm)
    fml <- reformulate(v.x.nms, response = y.nm)
    m.X <- model.matrix(fml, d)

    # align beta to columns of m.X (factor levels may add columns).
    v.beta.h.aligned <- rep(0, ncol(m.X));
    names(v.beta.h.aligned) <- colnames(m.X)
    v.beta.h.aligned[intersect(names(v.beta.hat), colnames(m.X))] <-
        v.beta.hat[intersect(names(v.beta.hat), colnames(m.X))]
    v.beta.d.aligned <- rep(0, ncol(m.X));
    names(v.beta.d.aligned) <- colnames(m.X)
    v.beta.d.aligned[intersect(names(v.beta.dp), colnames(m.X))] <-
        v.beta.dp[intersect(names(v.beta.dp), colnames(m.X))]

    v.eta.h  <- as.numeric(m.X %*% v.beta.h.aligned)
    v.eta.dp <- as.numeric(m.X %*% v.beta.d.aligned)
    v.mu.h   <- switch(link, identity = v.eta.h,  log = exp(v.eta.h))
    v.mu.dp  <- switch(link, identity = v.eta.dp, log = exp(v.eta.dp))

    d.virt <- d
    d.virt[[y.nm]] <- d[[y.nm]] + v.mu.dp - v.mu.h

    # Integer responses (poisson): round.  This is a projection onto
    # the count lattice.  Introduces O(1) rounding error per row but
    # preserves overall X'y_virtual ~= X'mu(beta_DP) to O(sqrt(n)).
    if(family == "poisson")
        d.virt[[y.nm]] <- pmax(0, round(d.virt[[y.nm]]))

    # Positive-support responses (gamma, quasi-gamma): clip to small
    # positive to avoid log(0).  Post-processing, so DP is preserved.
    if(family %in% c("gamma", "quasi-gamma"))
        d.virt[[y.nm]] <- pmax(1e-8, d.virt[[y.nm]])

    d.virt
}

# +------------------------------------------------------------+
# |  gen.syn.dp.projected: DP-projection generator             |
# +------------------------------------------------------------+
# End-to-end DP-projection generator.
#
# Arguments:
#   d         original data frame.
#   y.nm      response column name.
#   family    one of "gaussian", "gamma", "poisson".
#   link      one of "identity", "log".  Canonical links only.
#   supports  named vector specifying support of each variable
#             (see 11-projection-generator.r for options).
#   eps       noise budget epsilon > 0.  Not a privacy parameter here: see the
#             Caveat below.  It scales the noise and nothing more.
#   delta     as eps.
#   seed      RNG seed.
#
# The y.range argument is deliberately absent and must not be reintroduced, since
# a bound on |y| leaves a local-sensitivity calibration exactly as it was, and
# accepting one here would misrepresent the guarantee this route provides.
# Returns: list with
#   d.syn      synthetic dataset with OLS(y.nm ~ ...) = beta_DP
#   beta.hat   original MLE
#   beta.dp    DP-perturbed MLE (the target of projection)
#   sigma      Gaussian mechanism scale
#   Delta      sensitivity bound used
#   eta        noise vector (for diagnostics)
#
# Caveat: the sensitivity used here is local (DFBETA-based).  Local sensitivity
# does not provide (eps, delta)-DP, so this analysis-level generator carries no
# formal DP claim.  It is retained for the
# non-private inferential-fidelity experiments and as the projection engine.
# The sound DP routes are gen.syn.dp.ols.public (public X, global sensitivity)
# and gen.syn.dp.projected.subagg (subsample-and-aggregate), both mu-GDP.
gen.syn.dp.projected <- function(d, y.nm = "response",
                                  family = "gaussian",
                                  link = "identity",
                                  supports = NULL,
                                  eps = 1, delta = 1e-6,
                                  seed = 6489,
                                  standardise = TRUE,
                                  clip.sigma = Inf,
                                  sens.method = c("local", "smooth"),
                                  B.x = NULL, B.y = NULL) {
    # clip.sigma: if finite, the L2 norm of the DP noise eta is
    # truncated at clip.sigma * sigma (per-coordinate standardised).
    # This is Lipschitz post-processing, so DP is preserved.  For
    # log-link GLMs, clip.sigma=3 prevents overflow of mu(beta_DP)
    # at low eps.  For identity link, Inf (no clipping) is appropriate.
    # sens.method: "local" = DFBETA (empirical DP), "smooth" = NRS 2007
    # smooth sensitivity (formal (eps, delta)-DP when B.x, B.y are true
    # public bounds).
    # B.x, B.y: public L2 and L-infinity bounds for smooth sensitivity.
    sens.method <- match.arg(sens.method)
    if(sens.method == "local")
        warning("gen.syn.dp.projected (V1) is deprecated. It calibrates to the ",
                "DFBETA local sensitivity and carries no formal DP at any eps. ",
                "Its output is a disclosure metric, not a private release. ",
                "The sound routes are gen.syn.dp.projected.subagg (Route 2) and ",
                "gen.syn.dp.ols.public (Route 1).", call. = FALSE)
    if(!(family %in% c("gaussian", "gamma", "poisson", "quasi-gamma")))
        stop("family must be one of gaussian, gamma, poisson, quasi-gamma")
    if(!(link %in% c("identity", "log")))
        stop("link must be one of identity, log (canonical only)")

    v.x.nms <- setdiff(names(d), y.nm)
    fml <- reformulate(v.x.nms, response = y.nm)

    # Zero-inflation handling: for log-link with zeros, fit on
    # positive subset only (zeros carried through by copula).  DP is
    # for the positive-component coefficient.  Zero indicator is from
    # non-private copula (full-DP zero-model deferred to v2).
    b.zi <- (link == "log") && any(d[[y.nm]] == 0)
    d.fit <- `if`(b.zi, d[d[[y.nm]] > 0, , drop = FALSE], d)
    m.X <- model.matrix(fml, d.fit)

    # Step 1: fit the original model to obtain beta_hat.
    mo.orig <- tryCatch(switch(family,
        gaussian    = lm(fml, d.fit),
        gamma       = glm(fml, d.fit, family = Gamma(link = link)),
        poisson     = glm(fml, d.fit, family = poisson(link = link)),
        `quasi-gamma` = {
            s.np <- ncol(m.X)
            glm(fml, d.fit,
                family = quasi(link = "log", variance = "mu^2"),
                start = c(1, rep(0, s.np - 1)),
                control = glm.control(maxit = 200))
        }), error = function(e) NULL)
    if(is.null(mo.orig))
        stop("Failed to fit original model; check data and family/link.")
    v.beta.hat <- coef(mo.orig)

    # Step 2: column-standardise X for balanced sensitivity.
    # Set s_j = sd(X[,j]) for non-constant columns (1 otherwise).
    # Standardised coefficients: beta_tilde[j] = s_j * beta_hat[j].
    # After adding DP noise on the standardised scale and un-scaling,
    # the noise on beta_hat[j] scales as sigma / s_j.  This matches
    # per-coefficient uncertainty: covariates with high variance (high
    # information content) can absorb more noise than those with low
    # variance -- automatic scale matching.
    v.col.sd <- apply(m.X, 2, stats::sd)
    v.s <- ifelse(v.col.sd > 1e-8, v.col.sd, 1)
    names(v.s) <- colnames(m.X)
    v.s <- v.s[names(v.beta.hat)]                       # align order
    if(!standardise) v.s <- rep(1, length(v.beta.hat))
    v.beta.tilde <- v.beta.hat * v.s                    # standardised beta

    # Sensitivity of beta_tilde = sensitivity of (s * beta) via DFBETA
    # scaled column-wise.  Since dfbeta(mo) returns n x p with entry
    # (i, j) = beta_hat[j] - beta_hat[j, -i], we get:
    #   dfbeta_tilde[i, j] = s[j] * dfbeta[i, j]
    # Local sensitivity = 2 * max_i sqrt(sum_j dfbeta_tilde[i, j]^2).
    m.dfbeta <- dfbeta(mo.orig)
    m.dfbeta.tilde <- sweep(m.dfbeta, 2, v.s, "*")
    s.LS <- 2 * max(sqrt(rowSums(m.dfbeta.tilde^2)))
    s.S.smooth <- NA_real_
    if(sens.method == "smooth") {
        # Compute smooth sensitivity on the standardised scale.
        # Requires refitting lm on standardised X.
        d.std <- d
        v.x.cols <- setdiff(colnames(m.X), "(Intercept)")
        m.X.std <- sweep(m.X, 2, v.s, "/")              # standardised
        if(family == "gaussian" && link == "identity") {
            mo.std <- lm.fit(m.X.std, d[[y.nm]])
            mo.std.obj <- mo.orig                       # lm object
            # Re-fit via lm for dfbeta compatibility
            d.std.refit <- data.frame(m.X.std[, -1, drop = FALSE])
            d.std.refit[[y.nm]] <- d[[y.nm]]
            names(d.std.refit) <- c(colnames(m.X.std)[-1], y.nm)
            fml.std <- reformulate(colnames(m.X.std)[-1], response = y.nm)
            mo.std.obj <- lm(fml.std, d.std.refit)
            l.smooth <- smooth.sens.beta.ols(mo.std.obj, eps, delta,
                                              B.x = B.x, B.y = B.y)
            s.S.smooth <- 2 * l.smooth[["S.smooth"]]    # scale factor
            # Use smooth sens in place of LS if available and finite.
            if(is.finite(s.S.smooth)) s.LS <- s.S.smooth
        }
        # For log-link GLMs, smooth sens requires more care (non-linear
        # beta).  Fall back to local sensitivity with caveat.
    }
    s.Delta <- s.LS
    l.sens <- list(Delta = s.Delta,
                   Delta.local = s.LS,
                   Delta.smooth = s.S.smooth,
                   Delta.global = NA_real_,
                   sens.method = sens.method)

    # Step 3: add isotropic DP noise to beta_tilde, then un-scale.
    # Optional norm clipping (DP-preserving post-processing): if
    # ||eta||_2 > clip.sigma * sqrt(p) * sigma, rescale to the
    # radius.  With clip.sigma = 3, only ~ Pr(chi^2_p > 9p) ~ 1e-6
    # of samples get truncated in practice.  This caps runaway
    # perturbations that would overflow mu(beta_DP) for log links.
    set.seed(seed + 1L)
    l.dp <- dp.gaussian.noise(length(v.beta.hat), eps, delta, s.Delta)
    v.eta.standardised <- l.dp[["v.eta"]]
    s.p <- length(v.eta.standardised)
    s.eta.norm <- sqrt(sum(v.eta.standardised^2))
    s.eta.max <- clip.sigma * sqrt(s.p) * l.dp[["sigma"]]
    if(is.finite(s.eta.max) && s.eta.norm > s.eta.max)
        v.eta.standardised <- v.eta.standardised *
                              (s.eta.max / s.eta.norm)
    v.beta.tilde.dp <- v.beta.tilde + v.eta.standardised
    v.beta.dp <- v.beta.tilde.dp / v.s                  # un-standardise
    v.eta <- v.beta.dp - v.beta.hat                     # effective noise
    names(v.beta.dp) <- names(v.beta.hat)

    # Step 4: create virtual dataset where MLE = beta_DP.
    # For ZI: only modify positive y's.  Leave zeros alone so the
    # copula's categorical bootstrap of Z = I(Y>0) is preserved.
    if(b.zi) {
        # Inject virtual y only on positive subset, recombine with zeros.
        d.pos.virt <- inject.virtual.y(d.fit, y.nm, v.beta.hat,
                                        v.beta.dp, family, link)
        d.virt <- d
        v.pos.idx <- which(d[[y.nm]] > 0)
        d.virt[v.pos.idx, y.nm] <- d.pos.virt[[y.nm]]
    } else {
        d.virt <- inject.virtual.y(d, y.nm, v.beta.hat, v.beta.dp,
                                    family, link)
    }

    # Step 5: standard projection on the virtual dataset.
    # gen.syn.projected fits MLE on d.virt (= beta_DP), samples
    # synthetic X from copula, and projects synthetic y so that
    # MLE(synthetic) = beta_DP.  Output is DP by post-processing.
    d.syn <- gen.syn.projected(d.virt, y.nm = y.nm, family = family,
                                link = link, supports = supports,
                                seed = seed)

    list(d.syn        = d.syn,
         beta.hat     = v.beta.hat,
         beta.dp      = v.beta.dp,
         sigma        = l.dp[["sigma"]],
         Delta        = s.Delta,
         Delta.global = l.sens[["Delta.global"]],
         eta          = v.eta,
         col.scale    = v.s,
         eps          = eps,
         delta        = delta)
}

# +------------------------------------------------------------+
# |  gen.syn.dp.projected.zi: full DP-ZI with a private zero-model |
# +------------------------------------------------------------+
# Full differentially-private generator for zero-inflated log-link
# scenarios.  Where gen.syn.dp.projected's b.zi branch DP-protects only
# the positive-component coefficient and carries the zero indicator
# through the non-private copula bootstrap (Supp S6 residual-limitation
# 1), this closes the gap by privately releasing the zero-model too.
#
# Mechanism (fidelity-preserving order).  A DP-sampled zero indicator
# selects a different positive set than the copula did, so the positive
# values are projected after the indicator is fixed, not before.
#
#   1. Zero-model.  Fit logistic Z ~ X on the full data, Z = I(y > 0),
#      giving gamma_hat, and privatise it to gamma_DP by subsample-and-
#      aggregate (mu_Z-GDP, public clip box).
#   2. Positive component.  Fit the log-link GLM on d[y > 0], giving
#      beta_hat, and privatise it to beta_DP by subsample-and-aggregate
#      (mu_P-GDP, public clip box).
#   3. Virtual data.  Inject beta_DP into the positive rows so a log-link
#      refit on the release's positive subset returns beta_DP.  Zeros stay
#      0 (inject.virtual.y on the positive subset, recombined).
#   4. Synthetic covariates.  X_syn from the copula on the virtual full
#      data (gen.syn.family).  The copula response is kept only as a
#      projection warm-start.
#   5. Private zero indicator.  Z_syn_i ~ Bernoulli(expit(x_syn,i'
#      gamma_DP)).
#   6. Positive projection on the selected set.  On the Z_syn = 1 rows
#      only, warm-start the response (copula positive where available,
#      resampled from the copula positive pool otherwise) and project
#      (project.glm), reproducing beta_DP on the release's positive
#      subset, and y = 0 where Z_syn = 0.
#
# Privacy accounting (Remedy B, mu-GDP).  gamma_DP is mu_Z-GDP and beta_DP
# is mu_P-GDP, each by the subsample-and-aggregate Gaussian mechanism on a
# Public per-coordinate clip box (worst-case sensitivity, data-independent).
# By gdp.compose the pair is sqrt(mu_Z^2 + mu_P^2)-GDP, and the downstream
# synthetic assembly is post-processing, so the release inherits it.  The
# total mu_total = gdp.mu.from.eps.delta(eps, delta) is split in quadrature,
# mu_Z = mu_total sqrt(w_Z) and mu_P = mu_total sqrt(w_P) with w from
# mu.split, so mu_Z^2 + mu_P^2 = mu_total^2 charges exactly (eps, delta).
#
# Arguments:
#   x.nms       predictor names, default setdiff(names(d), y.nm).  Pass
#               explicitly to match a downstream analyst's model.
#   mu.split    named c(Z, P) fractions of the mu^2 (GDP) budget for the
#               zero-model and the positive component.  Must sum to 1.
#   m           subsample-and-aggregate block count.  Default NULL applies the
#               block-count rule per channel, floor(n / (5 p)) for the zero-model and
#               floor(n (1 - pi0) / (5 p)) for the positive fit, and then
#               Requires pos.rate (a dp.rate.public object carrying the
#               Published zero rate pi0, and dp.public.zero.rate(sc) serves it).
#               An explicit scalar m covers both channels as before.
#   clip.lo, clip.hi      public per-coordinate bounds on the positive-component
#               coefficients.  Required: they set s_j = (hi_j - lo_j)/m and
#               therefore the noise scale, so they must not be derived from d.
#   clip.z.lo, clip.z.hi  the same, for the logistic zero-model coefficients.
#
# The boxes are required arguments rather than an internal default because the
# noise scale is a function of the box alone.  A box built from the full-data
# fit would make the noise scale a function of the private data, a one-row
# change moving it by up to 5.7%.  That violates the
# requirement that one box be bound and the same box applied to both
# neighbours, so the mechanism had no DP guarantee at all.
# It also made the box centre equal to beta_hat, and failed blocks are imputed
# with the centre, so the non-private full-data estimate entered the release
# directly whenever a block fit failed.  Guarded by verify-gdp-mechanisms.R
# Task 11.
#
# Returns: list with d.syn, beta.hat, beta.dp, gamma.hat, gamma.dp, Z.syn,
#   the two noise scales, the public boxes, and the mu budgets (mu.total,
#   mu.Z, mu.P).
gen.syn.dp.projected.zi <- function(d, y.nm = "response",
                                    family = "quasi-gamma",
                                    link = "log",
                                    x.nms = NULL,
                                    eps = 1, delta = 1e-6,
                                    mu.split = c(Z = 0.5, P = 0.5),
                                    seed = 6489, m = NULL, pos.rate = NULL,
                                    clip.lo = NULL, clip.hi = NULL,
                                    clip.z.lo = NULL, clip.z.hi = NULL,
                                    dispersion = NULL) {
    if(is.null(clip.lo) || is.null(clip.hi))
        stop("clip.lo and clip.hi (public bounds on the positive-component ",
             "coefficients) are required")
    if(is.null(clip.z.lo) || is.null(clip.z.hi))
        stop("clip.z.lo and clip.z.hi (public bounds on the zero-model ",
             "coefficients) are required")
    if(dp.has.dispersion(family) && is.null(dispersion))
        stop("dispersion is required for the ", family, " family.  Supply ",
             "dp.dispersion.public(sigma2, source): the positive component has a ",
             "free dispersion parameter and DP must cover it.")
    if(!is.null(dispersion) && !inherits(dispersion, "dp.dispersion"))
        stop("dispersion must be built by dp.dispersion.public() or ",
             "dp.dispersion.private(), so that its provenance is on the record")
    if(!is.null(dispersion) && dispersion[["mode"]] == "private")
        stop("dp.dispersion.private() is not implemented for the zero-inflated ",
             "generator (its positive component is quasi-gamma, and the disp ",
             "channel of subsample.aggregate.beta.loglink is Gaussian-only).  ",
             "Supply dp.dispersion.public().")
    s.sigma2 <- `if`(is.null(dispersion), NA_real_, dispersion[["sigma2"]])
    if(link != "log")
        stop("gen.syn.dp.projected.zi handles the log link only.")
    if(!(family %in% c("gamma", "poisson", "quasi-gamma")))
        stop("family must be one of gamma, poisson, quasi-gamma")
    if(abs(sum(mu.split) - 1) > 1e-8)
        stop("mu.split must sum to 1")
    if(!any(d[[y.nm]] == 0))
        stop("no zeros in response; use gen.syn.dp.projected.subagg for non-ZI data.")

    if(is.null(x.nms)) x.nms <- setdiff(names(d), y.nm)

    # The block-count rule, per channel.  The zero-model fits every row, the
    # positive fit sees only the expected positive share, and both partitions
    # cover the fixed n rows, so both read the public nrow(d).  The positive
    # share must be the published rate: reading it off d would calibrate the
    # mechanism to the data it protects.
    if(is.null(m)) {
        if(is.null(pos.rate))
            stop("gen.syn.dp.projected.zi: with m = NULL (the block-count ",
                 "rule) pos.rate is required.  Supply dp.rate.public(pi0, ",
                 "source) with the PUBLISHED zero rate; dp.public.zero.rate(sc) ",
                 "serves it for the scenario designs.")
        if(!inherits(pos.rate, "dp.rate"))
            stop("pos.rate must be built by dp.rate.public(), so that its ",
                 "provenance is on the record")
        if(length(clip.lo) < 2L || length(clip.z.lo) < 2L)
            stop("gen.syn.dp.projected.zi: with m = NULL the clip boxes must ",
                 "be per-coefficient vectors (dp.public.box and ",
                 "dp.public.zero.box serve them): the rule reads p from the ",
                 "box length.")
        m.z <- dp.subagg.blocks(nrow(d), length(clip.z.lo))
        m.p <- dp.subagg.blocks(nrow(d), length(clip.lo), rate = pos.rate)
    } else {
        m.z <- m
        m.p <- m
    }

    # Split the total mu-GDP budget in quadrature: mu_Z^2 + mu_P^2 = mu_total^2.
    s.mu.total <- gdp.mu.from.eps.delta(eps, delta)
    s.mu.Z <- s.mu.total * sqrt(mu.split[["Z"]])
    s.mu.P <- s.mu.total * sqrt(mu.split[["P"]])

    v.pos <- d[[y.nm]] > 0
    d.pos <- d[v.pos, , drop = FALSE]

    # ---- Step 1: DP zero-model gamma_DP via subsample-and-aggregate ---
    # Logistic Z = I(y > 0) ~ X privatised to mu_Z-GDP by the assumption-free
    # subsample-and-aggregate mechanism, whose worst-case per-coordinate
    # sensitivity comes from the caller'S public clip box.  No full-data fit is
    # taken here: the box must not be a function of d, and subagg does its own
    # reference fit.
    d.z <- d[x.nms]
    d.z[[".Z"]] <- as.integer(v.pos)
    l.z <- subsample.aggregate.beta.loglink(d.z, ".Z", x.nms, "binomial",
                "logit", eps, delta, m = m.z, clip.lo = clip.z.lo,
                clip.hi = clip.z.hi, seed = seed + 100L, mu = s.mu.Z)
    v.gamma.dp <- l.z[["beta.dp"]]

    # ---- Step 2: DP positive-component beta_DP via subsample-and-aggregate -
    # subsample.aggregate.beta.loglink runs the same family/link/start logic for
    # its reference fit and errors if it fails, so no separate full-data fit is
    # needed here.
    #
    # The subsample-and-aggregate partition must be over a row set that is fixed
    # under one-row replacement.  Partitioning d.pos directly makes the partition
    # a function of |P| = sum(v.pos), which is sensitive (the zero-model channel
    # protects it): a replacement that toggles a row's positivity changes |P|,
    # reshuffles the whole partition, and moves the released aggregate by more
    # than the (hi_j - lo_j)/m the noise is calibrated to (attains about 3x the
    # sqrt(p) bound.  Verify-dp-soundness.R Task 6).  Instead partition the fixed
    # full N rows data-independently and restrict that assignment to the positive
    # rows, so a toggle changes exactly one block's positive subset and the
    # per-coordinate sensitivity stays (hi_j - lo_j)/m.
    set.seed(seed + 3L)
    v.block.full <- sample(rep_len(seq_len(m.p), nrow(d)))
    l.p <- subsample.aggregate.beta.loglink(d.pos, y.nm, x.nms, family,
                "log", eps, delta, m = m.p, clip.lo = clip.lo,
                clip.hi = clip.hi, seed = seed, mu = s.mu.P,
                block = v.block.full[v.pos])
    v.beta.hat <- l.p[["beta.hat"]]
    v.beta.dp  <- l.p[["beta.dp"]]

    # ---- Step 3: private zero indicator Z_syn -------------------------
    # Z_syn ~ Bernoulli(expit(X gamma_DP)): a function of gamma_DP and the (public)
    # covariates alone.  This channel was always sound -- a textbook
    # post-processing of the DP output.  It is the pattern the response channel
    # below should have followed from the start, and did not.
    #
    # X is released as it stands, and like Route 2 this generator assumes public
    # covariates.
    m.X.syn <- model.matrix(reformulate(x.nms), d)
    v.gamma.aligned <- `names<-`(rep(0, ncol(m.X.syn)), colnames(m.X.syn))
    v.common.z <- intersect(names(v.gamma.dp), colnames(m.X.syn))
    v.gamma.aligned[v.common.z] <- v.gamma.dp[v.common.z]
    v.p.syn <- plogis(as.numeric(m.X.syn %*% v.gamma.aligned))
    set.seed(seed + 2L)
    v.Z.syn <- rbinom(nrow(m.X.syn), 1L, v.p.syn)

    # ---- Step 4: draw the positive response from the family -----------
    # An earlier design warm-started the positive rows from the copula (fed the real
    # y) and then projected against d.virt, whose sufficient statistics carried
    # the real RSS, so the release reproduced the real residual distribution and
    # broke the subsample-and-aggregate proposition's post-processing
    # hypothesis.  Draw instead, from the model's own family at (beta_DP,
    # dispersion), and zero the rows the private zero-model switched off.  The
    # real y is not an input below.
    d.syn <- synth.response(d, y.nm, x.nms, family, "log",
                            v.beta.dp, s.sigma2, seed)
    d.syn <- `[[<-`(d.syn, y.nm,
                    value = `[<-`(d.syn[[y.nm]], v.Z.syn == 0L, 0))

    # beta.hat, gamma.hat, beta.agg and the n.failed counts are non-private
    # intermediates, returned for Monte Carlo evaluation only.  They are not
    # part of the release and must not be reported as such.
    list(d.syn        = d.syn,
         beta.hat     = v.beta.hat,
         beta.dp      = v.beta.dp,
         gamma.hat    = l.z[["beta.hat"]],
         gamma.dp     = v.gamma.dp,
         Z.syn        = v.Z.syn,
         sigma2       = s.sigma2,
         dispersion   = dispersion,
         sigma.beta   = mean(l.p[["sigma"]]),
         sigma.gamma  = mean(l.z[["sigma"]]),
         sigma.beta.j  = l.p[["sigma"]],
         sigma.gamma.j = l.z[["sigma"]],
         s.j.beta     = l.p[["s.j"]],
         s.j.gamma    = l.z[["s.j"]],
         clip.lo      = l.p[["clip.lo"]],
         clip.hi      = l.p[["clip.hi"]],
         clip.z.lo    = l.z[["clip.lo"]],
         clip.z.hi    = l.z[["clip.hi"]],
         mu.total     = gdp.compose(c(s.mu.Z, s.mu.P)),
         mu.Z         = s.mu.Z,
         mu.P         = s.mu.P,
         m            = c(z = m.z, p = m.p),
         n.failed.Z   = l.z[["n.failed"]],
         n.failed.P   = l.p[["n.failed"]],
         eps          = eps,
         delta        = delta)
}

# +------------------------------------------------------------+
# |  subsample.aggregate.beta.loglink: worst-case DP (DP-D2)   |
# +------------------------------------------------------------+
# Worst-case (eps, delta)-DP coefficient for a log-link GLM via
# subsample-and-aggregate (Nissim-Raskhodnikova-Smith 2007).  Where
# dp.perturb.beta calibrates Gaussian noise to the DFBETA local
# sensitivity (empirical DP, Supp S6 residual-limitation 2, because the
# log-link MLE has no closed smooth-sensitivity bound), this route gives
# a true worst-case guarantee that reuses the existing per-block log-link
# fits.
#
# Mechanism.
#   1. Partition the n rows into m disjoint blocks by a data-independent
#      random assignment of row positions (seeded).  Replacing one row's
#      values keeps it in its block, so it perturbs exactly one block.
#   2. Fit the log-link GLM on each block, giving beta^(1..m) (non-private
#      intermediate quantities, not released).
#   3. Widened Winsorized Mean per coordinate j: clip each beta^(b)_j to a
#      Public interval [lo_j, hi_j] (failed blocks imputed with the box
#      centre, also in the interval), then average the m clipped values.
#      Replacing one row moves one clipped value by at most (hi_j - lo_j),
#      so the mean a_j moves by at most (hi_j - lo_j)/m: a worst-case,
#      data-independent per-coordinate sensitivity.
#   4. Gaussian mechanism on the scaled aggregate.  With s_j =
#      (hi_j - lo_j)/m the scaled query Q_j = a_j/s_j has L2 sensitivity
#      <= sqrt(p) (one block change moves each scaled coordinate by <= 1),
#      so isotropic noise N(0, sigma_Q^2 I) with sigma_Q from the exact
#      GDP calibration (gdp.sigma) is (sqrt(p)/sigma_Q)-GDP, converted to
#      (eps, delta) by the exact dual.  Un-scaling by s_j (post-processing)
#      gives beta_DP_j = a_j + N(0, (s_j sigma_Q)^2).
#
# The m block fits are non-private, so the whole DP charge is the single
# Gaussian mechanism.  The public clip box is the only public input the
# guarantee requires (analogous to the OLS smooth-sensitivity B.x, B.y).
# No data-dependent sensitivity anywhere.
#
# The noise multiplier is the exact mu-GDP value from gdp.sigma (the
# Balle-Wang / Dong-Roth-Su calibration), not the classical
# sqrt(2 ln(1.25/delta))/eps constant.  Only the sensitivity is specific to
# this route, and that is the worst-case per-coordinate contribution.
#
# Arguments:
#   d, y.nm, x.nms   data, response name, predictor names.
#   family, link     GLM family and link.  gaussian/identity, gamma/log,
#                    poisson/log, quasi-gamma/log, and binomial/logit (the
#                    zero-inflated zero-model) are supported.
#   eps, delta       worst-case DP budget (ignored when mu is supplied).
#   m                block count.  Default NULL applies the block-count rule
#                    dp.subagg.blocks(n, p) = floor(n / (5 p)).  Required
#                    explicitly when `block` is supplied (see the guard below).
#   clip.lo, clip.hi public per-coefficient bounds (named vectors aligned
#                    to coef names, or scalars recycled).  Required.
#   mu               optional GDP budget.  When given, sigma.Q = sqrt(p)/mu
#                    directly, so two releases compose in quadrature (used by
#                    the zero-inflated generator).  Else from (eps, delta).
#   seed             RNG seed.  Partition under set.seed(seed), noise under
#                    set.seed(seed + 1L).
#   block            optional length-nrow(d) block assignment (values in 1..m).
#                    When supplied it replaces the internal seeded partition, so
#                    the caller can partition a fixed superset of the rows: the
#                    zero-inflated generator partitions the full N and restricts
#                    that assignment to the positive rows, which keeps the
#                    partition data-independent even when a one-row replacement
#                    toggles a row's positivity and so changes |P|.
# Returns: list(beta.hat, beta.dp, beta.agg, sigma, sigma.Q, s.j,
#               mu.total, m, n.failed, clip.lo, clip.hi, block).
# disp = TRUE appends the OLS dispersion sigma^2 as a (p+1)-th aggregated
# coordinate, with its own public box (disp.lo, disp.hi).  This makes the
# release sound when X is sensitive: the residual-norm sigma^2 channel of
# gen.syn.dp.full assumed response-row replacement holding X fixed, which is
# the wrong neighbouring relation once a neighbour may change the whole row.
# Here every coordinate, sigma^2 included, is clipped into a public box over a
# data-independent partition, so a one-row change moves the mean by at most
# s_j = (hi_j - lo_j)/m whether it altered x, y, or both.  No design constant,
# nothing data-dependent.  The scaled aggregate then has L2 sensitivity
# sqrt(p + 1) and the release is (sqrt(p + 1) / sigma.Q)-GDP.
#
# Not private, evaluation only: beta.hat is the full-data fit, and beta.agg
# and n.failed are data-dependent intermediates.  The DP release is beta.dp
# (and sigma2.dp).  Do not compose the others into a downstream release.
subsample.aggregate.beta.loglink <- function(d, y.nm, x.nms, family, link,
                                              eps, delta, m = NULL,
                                              clip.lo = NULL, clip.hi = NULL,
                                              seed = 6489, mu = NULL,
                                              block = NULL,
                                              disp = FALSE,
                                              disp.lo = NULL, disp.hi = NULL) {
    if(is.null(clip.lo) || is.null(clip.hi))
        stop("clip.lo and clip.hi (public per-coefficient bounds) are required")
    if(disp && (is.null(disp.lo) || is.null(disp.hi)))
        stop("disp.lo and disp.hi (public bounds on sigma^2) are required when disp = TRUE")
    if(disp && family != "gaussian")
        stop("disp = TRUE is defined only for the gaussian family (OLS dispersion)")
    fml <- reformulate(x.nms, response = y.nm)
    n <- nrow(d)

    fit.block <- function(dd) tryCatch(switch(family,
        gaussian      = lm(fml, dd),
        binomial      = glm(fml, dd, family = binomial(link = link)),
        gamma         = glm(fml, dd, family = Gamma(link = link)),
        poisson       = glm(fml, dd, family = poisson(link = link)),
        `quasi-gamma` = {
            s.np <- ncol(model.matrix(fml, dd))
            glm(fml, dd, family = quasi(link = "log", variance = "mu^2"),
                start = c(1, rep(0, s.np - 1)),
                control = glm.control(maxit = 200))
        }), error = function(e) NULL)

    mo.ref <- fit.block(d)
    if(is.null(mo.ref)) stop("reference fit failed; check data and family/link.")
    v.beta.hat <- coef(mo.ref)
    v.names <- names(v.beta.hat)
    s.p <- length(v.names)
    v.qnames <- if(disp) c(v.names, ".sigma2") else v.names
    s.q <- length(v.qnames)

    # align the public clip bounds to coefficient names (recycle scalars).
    align.clip <- function(cl) {
        if(is.null(names(cl))) return(`names<-`(rep_len(cl, s.p), v.names))
        v <- cl[v.names]
        if(any(is.na(v))) stop("named clip bounds must cover every coefficient")
        v
    }
    v.lo <- align.clip(clip.lo); v.hi <- align.clip(clip.hi)
    if(disp) {
        v.lo <- `names<-`(c(v.lo, disp.lo), v.qnames)
        v.hi <- `names<-`(c(v.hi, disp.hi), v.qnames)
    }
    if(any(v.hi <= v.lo))
        stop("clip.hi must exceed clip.lo for every coordinate")
    v.center <- (v.lo + v.hi) / 2

    # The block-count rule.  Self-default only on a self-partitioned call: a
    # supplied block assignment means the caller partitioned a fixed superset,
    # so nrow(d) here is a restricted, sensitive row count the rule must not
    # read.  On that path m stays a required argument, computed upstream from
    # published constants.
    if(is.null(m)) {
        if(!is.null(block))
            stop("subsample.aggregate.beta.loglink: m is required when block ",
                 "is supplied.  The block-count rule reads nrow(d), which on ",
                 "a restricted row set is sensitive; compute m upstream from ",
                 "published constants (dp.subagg.blocks) and pass it.")
        m <- dp.subagg.blocks(n, s.p)
    }

    # Partition the rows into m blocks.  When the caller supplies `block` (the
    # zero-inflated generator computes it over the fixed full N and restricts it
    # to the positive rows), use it verbatim, so a one-row replacement changes
    # exactly one block even when it toggles a row's positivity and so changes
    # |P|.  Otherwise partition by row position here, which is data-independent
    # only when the row set is fixed (Route 2 on the full n).
    if(is.null(block)) {
        set.seed(seed)
        v.block <- sample(rep_len(seq_len(m), n))
    } else {
        if(length(block) != n)
            stop("block must have length nrow(d) = ", n)
        v.block <- block
    }

    # per-block fit -> (beta^(b), sigma2^(b)) (non-private intermediates).
    l.beta <- Map(function(b) {
        mo.b <- fit.block(d[v.block == b, , drop = FALSE])
        v.out <- `names<-`(rep(NA_real_, s.q), v.qnames)
        if(!is.null(mo.b)) {
            cf <- coef(mo.b)
            v.c <- intersect(v.names, names(cf))
            v.out[v.c] <- cf[v.c]
            # residual df can vanish in a thin block.  NA then falls through to
            # the public box centre below, exactly as a failed fit does.
            if(disp) v.out[[".sigma2"]] <- tryCatch(stats::sigma(mo.b)^2,
                                                    error = function(e) NA_real_)
        }
        v.out
    }, seq_len(m))
    m.beta <- do.call(rbind, l.beta)                    # m x (p or p+1)

    # Widened Winsorized Mean: clip to the public box, impute failed
    # blocks with the (public) box centre, average coordinate-wise.
    m.clip <- mapply(function(j) {
        v <- pmin(pmax(m.beta[, j], v.lo[j]), v.hi[j])
        v[is.na(v) | !is.finite(v)] <- v.center[j]
        v
    }, seq_len(s.q))
    if(is.null(dim(m.clip))) m.clip <- matrix(m.clip, nrow = m)
    v.agg <- `names<-`(colMeans(m.clip), v.qnames)

    # worst-case Gaussian mechanism on the scaled aggregate, calibrated in
    # mu-GDP: the scaled aggregate has L2 sensitivity sqrt(q), so the release
    # is (sqrt(q) / sigma.Q)-GDP, with q = p (+ 1 when the dispersion rides
    # along).  Supplying mu directly lets this compose in quadrature.
    v.s <- (v.hi - v.lo) / m
    s.sigma.Q <- if(!is.null(mu)) sqrt(s.q) / mu else gdp.sigma(sqrt(s.q), eps, delta)
    set.seed(seed + 1L)
    v.rel <- `names<-`(v.agg + v.s * s.sigma.Q * rnorm(s.q), v.qnames)

    list(beta.hat = v.beta.hat,
         beta.dp = v.rel[v.names],
         sigma2.dp = if(disp) max(v.rel[[".sigma2"]], 1e-6) else NA_real_,
         beta.agg = v.agg[v.names],
         sigma = v.s * s.sigma.Q, sigma.Q = s.sigma.Q, s.j = v.s,
         mu.total = sqrt(s.q) / s.sigma.Q,
         m = m, n.failed = sum(is.na(m.beta[, 1])),
         clip.lo = v.lo, clip.hi = v.hi, block = v.block)
}

# +------------------------------------------------------------+
# |  gen.syn.dp.projected.subagg: worst-case DP generator (D2) |
# +------------------------------------------------------------+
# End-to-end worst-case DP-projection generator for canonical-link GLMs.  The
# coefficient comes from subsample-and-aggregate on a public clip box.  The
# response is then drawn from the model's family at (beta_DP, dispersion), so the
# release is a function of the DP output and public information alone.
#
# Scope: log-link GLMs (scenarios 3, 6) and Gaussian/identity (scenarios
# 1, 2, 9, 10), matching v.route2 in 62-mc-dp-gdp.r.
#
# Route 2 assumes X is public.  That is the route's standing assumption (with
# sensitive X you need V4/V5, which add a DP copula on the covariates), and the
# covariates are therefore released as they are.  Releasing X plainly makes the
# assumption visible in the output.
#
# Two design points:
#
#   1. There is no virtual-y projection tail.  inject.virtual.y retains the real
#      y, so a release built that way carries the real residuals in scale and in
#      shape, and cannot claim the subsample-and-aggregate proposition's
#      post-processing clause, whose hypothesis is a release depending on the
#      data "only through beta_DP".  An empirical audit scores such a release at
#      AUC 1.000 against a mu-GDP bound of 0.566.  See synth.response.
#
#   2. `dispersion` is now a required argument for the families that have one.
#      A GLM has two parameters and DP must cover both.  Say where sigma^2 comes
#      from: dp.dispersion.public() if it is external or already published (no budget,
#      no budget, beta_DP exactly as noisy as before), dp.dispersion.private() if
#      it must be paid for (costs ~8% more coefficient noise at p = 6).
#
# There is deliberately no y.range argument.  Subsample-and-aggregate takes its
# sensitivity from clipping the per-block coefficients into the public box, so it
# holds however large or unbounded y is.
# That is the reason the mechanism exists (Nissim, Raskhodnikova and Smith
# designed it for statistics with no bounded global sensitivity).  A bound on
# |y| is needed by Route 1, not here.
gen.syn.dp.projected.subagg <- function(d, y.nm = "response",
                                         family = "gamma", link = "log",
                                         x.nms = NULL,
                                         eps = 1, delta = 1e-6, m = NULL,
                                         clip.lo = NULL, clip.hi = NULL,
                                         dispersion = NULL, seed = 6489) {
    if(is.null(x.nms)) x.nms <- setdiff(names(d), y.nm)

    b.needs <- dp.has.dispersion(family)
    if(b.needs && is.null(dispersion))
        stop("dispersion is required for the ", family, " family: it has a free ",
             "dispersion parameter, and DP must cover it.  Supply ",
             "dp.dispersion.public(sigma2, source) if sigma^2 is external or ",
             "already published, or dp.dispersion.private(lo, hi, source) to ",
             "spend budget on it.")
    if(!is.null(dispersion) && !inherits(dispersion, "dp.dispersion"))
        stop("dispersion must be built by dp.dispersion.public() or ",
             "dp.dispersion.private(), so that its provenance is on the record")

    b.priv <- b.needs && dispersion[["mode"]] == "private"
    if(b.priv && family != "gaussian")
        stop("dp.dispersion.private() is implemented for the gaussian family ",
             "only (subsample.aggregate.beta.loglink's disp channel is the OLS ",
             "dispersion).  For ", family, " supply dp.dispersion.public().")

    l.sa <- subsample.aggregate.beta.loglink(d, y.nm, x.nms, family, link,
                eps, delta, m = m, clip.lo = clip.lo, clip.hi = clip.hi,
                seed = seed, disp = b.priv,
                disp.lo = `if`(b.priv, dispersion[["lo"]], NULL),
                disp.hi = `if`(b.priv, dispersion[["hi"]], NULL))

    s.sigma2 <- `if`(!b.needs, NA_real_,
                `if`(b.priv, l.sa[["sigma2.dp"]], dispersion[["sigma2"]]))

    # The release: X as it stands (public, by the route's assumption) and a
    # response drawn from the family at (beta_DP, sigma^2).  The real y is not an
    # input to this line, which is the point.
    d.syn <- synth.response(d, y.nm, x.nms, family, link,
                            l.sa[["beta.dp"]], s.sigma2, seed)

    # beta.hat, beta.agg and n.failed are non-private intermediates, returned for
    # Monte Carlo evaluation only.  They are not part of the release.
    list(d.syn     = d.syn,
         beta.hat  = l.sa[["beta.hat"]],
         beta.dp   = l.sa[["beta.dp"]],
         beta.agg  = l.sa[["beta.agg"]],
         sigma2    = s.sigma2,
         dispersion = dispersion,
         mu.total  = l.sa[["mu.total"]],
         sigma     = l.sa[["sigma"]],
         s.j       = l.sa[["s.j"]],
         m         = l.sa[["m"]],
         n.failed  = l.sa[["n.failed"]],
         route     = "subsample-aggregate",
         eps       = eps,
         delta     = delta,
         clip.lo   = l.sa[["clip.lo"]],
         clip.hi   = l.sa[["clip.hi"]])
}

# +------------------------------------------------------------+
# |  dp.project: unified wrapper                               |
# +------------------------------------------------------------+
# Auto-selects one of the four DP-Projection variants based on
# user inputs:
#
#   dp.scope    "analysis"  -> gen.syn.dp.projected (any X, leaks via
#                              copula, no formal DP guarantee, a
#                              local-sensitivity fidelity engine)
#               "full"      -> gen.syn.dp.ols.public (public X) or
#                              gen.syn.dp.full (sensitive X, via x.spec)
#   target      "ols"       -> default OLS estimator
#               "ridge"     -> gen.syn.dp.ridge.public (requires lambda)
#   x.public    TRUE/FALSE  -> whether X is publicly releasable
#   x.spec      list(...)   -> DP marginal specs (for sensitive-X full DP)
#
# The wrapper enforces the combinations that are supported:
#   - ridge target requires x.public = TRUE (public-X only)
#   - full DP sensitive X requires x.spec
#   - Gamma/Poisson only supported under "analysis" scope
#
# Delegates to the appropriate variant and returns its result.
# Requires the public-X variant source files (33, 34, 35) to be loaded
# when the corresponding scope/target is requested.
dp.project <- function(d, y.nm, family = "gaussian", link = "identity",
                        dp.scope = c("analysis", "full"),
                        target = c("ols", "ridge", "logistic"),
                        x.public = FALSE, x.spec = NULL,
                        lambda = 1, supports = NULL,
                        eps = 10, delta = 1e-6,
                        B.y = NULL, seed = 6489,
                        ...) {
    dp.scope <- match.arg(dp.scope)
    target <- match.arg(target)

    if(target == "logistic") {
        if(!x.public && is.null(x.spec))
            stop("logistic sensitive-X requires x.spec (or set x.public = TRUE)")
        if(!exists("gen.syn.dp.logistic"))
            stop("source('53-dp-logistic.r') first")
        return(gen.syn.dp.logistic(d, y.nm = y.nm, eps = eps,
                                    delta = delta, x.public = x.public,
                                    x.spec = x.spec, seed = seed, ...))
    }

    if(target == "ridge") {
        if(!x.public) stop("ridge target requires public X (use x.public = TRUE)")
        if(family != "gaussian" || link != "identity")
            stop("ridge target requires gaussian identity link")
        if(!exists("gen.syn.dp.ridge.public"))
            stop("source('34-dp-ridge-projection.r') first")
        if(is.null(B.y)) stop("B.y required")
        return(gen.syn.dp.ridge.public(d, y.nm = y.nm, lambda = lambda,
                                        eps = eps, delta = delta, B.y = B.y,
                                        seed = seed, ...))
    }

    # dp.scope = "analysis" dispatches to V1, which is deprecated and carries no
    # formal DP, calibrating to the DFBETA local sensitivity.  It is kept as
    # the no-guarantee baseline the paper contrasts the sound routes against,
    # and gen.syn.dp.projected warns on every call.  It takes no bound on |y|,
    # because no bound on |y| would make a local-sensitivity calibration sound.
    if(dp.scope == "analysis") {
        return(gen.syn.dp.projected(d, y.nm = y.nm, family = family,
                                     link = link, supports = supports,
                                     eps = eps, delta = delta,
                                     seed = seed, ...))
    }

    # full DP
    if(family != "gaussian" || link != "identity")
        stop("full DP currently requires gaussian identity link; ",
             "use dp.scope = 'analysis' for other families")
    if(is.null(B.y)) stop("B.y required for full DP")

    if(x.public) {
        if(!exists("gen.syn.dp.ols.public"))
            stop("source('33-dp-projection-public-x.r') first")
        return(gen.syn.dp.ols.public(d, y.nm = y.nm, eps = eps,
                                      delta = delta, B.y = B.y,
                                      seed = seed, ...))
    }

    # sensitive X full DP
    if(is.null(x.spec))
        stop("x.spec required for full DP with sensitive X; ",
             "see gen.syn.dp.full() docs")
    if(!exists("gen.syn.dp.full"))
        stop("source('35-dp-full.r') first")
    gen.syn.dp.full(d, y.nm = y.nm, x.spec = x.spec, eps.total = eps,
                     delta = delta, B.y = B.y, seed = seed, ...)
}

# +------------------------------------------------------------+
# |  Self-check when sourced directly                        |
# +------------------------------------------------------------+
if(!exists(".dp.suppress.demo", inherits = FALSE)) {
    cat("28-dp-projection.r loaded.\n")
    cat("  Main entry: dp.project(), which selects the variant.\n")
    cat("  Variants:   gen.syn.dp.projected, .ols.public, .ridge.public, .full\n")
}
