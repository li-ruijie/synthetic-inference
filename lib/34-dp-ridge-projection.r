#!/usr/bin/env r
# +------------------------------------------------------------+
# |  DP-Ridge-Projection (for ill-conditioned / sparse strata) |
# +------------------------------------------------------------+
#
# When the design matrix X is ill-conditioned (sparse categorical
# strata, n/p near 1, high leverage), OLS local sensitivity LS(beta)
# is dominated by extreme observations, giving DP noise that destroys
# coefficient estimates (cf. PHI raw: max|%| = 6115% at eps=10 with
# gen.syn.dp.projected).
#
# This variant uses ridge regression as the target estimator:
#   beta_ridge = (X'X + lambda * I)^{-1} X'y
# which has bounded sensitivity via the regularisation term.  Ridge
# DFBETA is bounded by 2*B_x*B_y / (lambda_min(X'X) + lambda), and the
# lambda term provides a floor: sensitivity <= 2*B_x*B_y / lambda.
#
# The target changes: analyst running ridge(lambda) on synthetic data
# recovers beta_ridge_DP (not beta_OLS_DP).  Suitable when:
#   - Pre-registered analysis uses ridge anyway (common in prediction,
#     not in inference).
#   - Analyst accepts regularisation bias for DP tractability.
#   - Sparse strata preclude stable OLS.
#
# Depends on: base R, plus dp.gm.multiplier() (28-dp-projection.r).
if(!exists("dp.gm.multiplier")) {
    .dp.suppress.demo <- TRUE; source("lib/28-dp-projection.r")
}

rlaplace <- function(n, location = 0, scale = 1) {
    v.u <- runif(n) - 0.5
    location - scale * sign(v.u) * log(1 - 2 * abs(v.u))
}

# +------------------------------------------------------------+
# |  Ridge coefficient sensitivity, public X (the sound one)    |
# +------------------------------------------------------------+
# Global sensitivity of beta_ridge under response-row replacement with a public
# design, which is Route 1's threat model (X public, |y_i| <= B.y).
#
# Derivation, and it is two lines.  With A = X'X + lambda I,
#
#     beta_ridge = A^{-1} X' y.
#
# Replacing row i's response by y_i' moves X'y by x_i (y_i' - y_i) and nothing
# else, so
#
#     beta' - beta = A^{-1} x_i (y_i' - y_i),
#     ||beta' - beta||_2 = |y_i' - y_i| ||A^{-1} x_i||_2 <= 2 B_y ||A^{-1} x_i||_2.
#
# Hence Delta_beta = 2 B_y max_i ||A^{-1} x_i||_2, which depends on X, lambda and
# B.y only.  X is public here by the route's premise, so this is a genuine global
# sensitivity.  It is exactly gs.beta.ols.public with A in place of X'X, and at
# lambda -> 0 it reduces to it (asserted in verify-gdp-mechanisms.R Task 6).
gs.beta.ridge.public <- function(m.X, lambda, B.y) {
    m.A.inv <- solve(crossprod(m.X) + lambda * diag(ncol(m.X)))
    2 * B.y * max(sqrt(rowSums((m.X %*% m.A.inv)^2)))
}

# +------------------------------------------------------------+
# |  Deprecated: ridge DFBETA local sensitivity                 |
# +------------------------------------------------------------+
# This is invalid as a DP calibration, and a sound path must not call it.
#
# It is the local sensitivity of beta_ridge, read off the realised residuals
# (v.e <- v.y - v.mu below).  A sensitivity that is a function of the data it is
# meant to protect gives neighbouring datasets different noise scales, so no
# valid bound exists and there is no theorem to instantiate at any epsilon.
#
# Kept only for reproducibility.  The sound alternative is gs.beta.ridge.public.
sens.beta.ridge <- function(m.X, v.y, lambda) {
    m.XtX <- crossprod(m.X)
    m.inv <- solve(m.XtX + lambda * diag(ncol(m.X)))
    v.beta <- m.inv %*% crossprod(m.X, v.y)
    v.mu <- as.numeric(m.X %*% v.beta)
    v.e <- v.y - v.mu                        # <-- the realised residuals. The bug.
    m.H <- m.X %*% m.inv %*% t(m.X)
    v.h <- diag(m.H)
    # leverage stays < 1 with lambda > 0
    m.dfb <- sweep(m.X %*% m.inv, 1, v.e / (1 - v.h), "*")
    list(beta = as.numeric(v.beta),
         LS   = 2 * max(sqrt(rowSums(m.dfb^2))))
}

# +------------------------------------------------------------+
# |  Ridge variance sensitivity (public X, |y_i| <= B.y)       |
# +------------------------------------------------------------+
# sigma_hat^2 = ||y - X beta_ridge||^2 / (n - p) = y'My / (n - p) with
# M = (I - H_l)'(I - H_l) and H_l = X (X'X + lambda I)^{-1} X'.  With
# A = X'X + lambda I this expands to M = I - X A^{-1} (X'X + 2 lambda I)
# A^{-1} X', an O(n^2 p) form that avoids the n x n product.  The one-row
# identity behind gs.sigma2 (28-dp-projection.r) holds verbatim with M in
# place of I - H: RSS(y') - RSS(y) = 2 d (My)_i + d^2 M_ii, so over the
# bounded universe |y_j| <= B.y the worst case is
# |Delta RSS| = 4 B.y^2 max_i max{M_ii, sum_{j != i} |M_ij|}.  At lambda = 0
# this reduces to gs.sigma2's hat-row constant K_X.  The naive 4 B.y^2 / (n - p)
# (row constant 1) would under-noise the variance channel.
gs.sigma2.ridge <- function(m.X, lambda, B.y) {
    n <- nrow(m.X); p <- ncol(m.X)
    m.XtX <- crossprod(m.X)
    m.A.inv <- solve(m.XtX + lambda * diag(p))
    m.M <- diag(n) - m.X %*% m.A.inv %*% (m.XtX + 2 * lambda * diag(p)) %*%
        m.A.inv %*% t(m.X)
    v.m <- diag(m.M)
    v.s <- rowSums(abs(m.M)) - abs(v.m)
    4 * B.y^2 * max(pmax(v.m, v.s)) / (n - p)
}

# +------------------------------------------------------------+
# |  DP-Ridge synthesiser (public X)                           |
# +------------------------------------------------------------+
# Arguments:
#   d         original data frame
#   y.nm      response column
#   x.nms     predictor names (if NULL, all non-y columns)
#   lambda    ridge parameter (positive)
#   eps       total privacy budget
#   delta     privacy budget delta
#   B.y       public bound on |y|
#   eps.split fraction of eps allocated to beta (vs sigma^2)
#   seed      RNG seed
# Returns: list with d.syn, beta.dp, sigma2.dp, etc.
gen.syn.dp.ridge.public <- function(d, y.nm = "response", x.nms = NULL,
                                      lambda = 1, eps = 1, delta = 1e-6,
                                      B.y = NULL, eps.split = 0.8,
                                      seed = 6489) {
    if(is.null(B.y)) stop("B.y required")
    if(lambda <= 0) stop("lambda must be > 0 for ridge")
    if(is.null(x.nms)) x.nms <- setdiff(names(d), y.nm)
    fml <- reformulate(x.nms, response = y.nm)
    m.X <- model.matrix(fml, d)
    n <- nrow(m.X); p <- ncol(m.X)
    v.y <- d[[y.nm]]

    # Fit ridge
    m.XtX <- crossprod(m.X)
    m.inv <- solve(m.XtX + lambda * diag(p))
    v.beta.ridge <- as.numeric(m.inv %*% crossprod(m.X, v.y))
    v.mu.r <- as.numeric(m.X %*% v.beta.ridge)
    v.e.r <- v.y - v.mu.r
    s.sigma2.r <- sum(v.e.r^2) / (n - p)                # approx

    # +--------------------------------------------------------------------+
    # |  Both channels use a public global sensitivity                     |
    # +--------------------------------------------------------------------+
    # A local sensitivity read off the realised residuals is a function of the
    # sensitive data, so the two neighbours get different noise scales and no DP
    # theorem can be instantiated.  The design here is public, so a global
    # sensitivity is always available (see gs.beta.ridge.public).
    #
    # No column standardisation is applied.  It would balance DFBETA noise across
    # columns of very different scale, but the global sensitivity already bounds
    # the raw-scale L2 movement of beta, so the noise is isotropic on the raw
    # scale, exactly as V2 does it.  No 3-sigma truncation is applied either,
    # since a Gaussian mechanism at a global sensitivity needs no such step.
    #
    # The accounting is mu-GDP, matching V2.  A Laplace channel does not enter a
    # mu quadrature at all, so both channels are Gaussian and compose in one
    # quadrature, mu^2 = mu_beta^2 + mu_v^2, converted to (eps, delta) by the
    # exact dual.
    s.mu.total <- gdp.mu.from.eps.delta(eps, delta)
    s.mu.beta  <- s.mu.total * sqrt(eps.split)
    s.mu.v     <- s.mu.total * sqrt(1 - eps.split)

    # Step 1: DP beta, Gaussian mechanism at the public-X global sensitivity
    # Delta_beta = 2 B.y max_i ||(X'X + lambda I)^{-1} x_i||.
    s.Delta.beta <- gs.beta.ridge.public(m.X, lambda, B.y)
    s.sigma.beta <- s.Delta.beta / s.mu.beta
    set.seed(seed + 1L)
    v.beta.dp <- `names<-`(v.beta.ridge + rnorm(p, 0, s.sigma.beta),
                           colnames(m.X))

    # Step 2: DP sigma^2, Gaussian mechanism at the exact global sensitivity
    # Delta_v = 4 B.y^2 K_M / (n - p), K_M the row constant of
    # M = (I - H_lambda)'(I - H_lambda) computed from the public X
    # (gs.sigma2.ridge).  The 1e-6 floor is the numerical form of the truncation
    # at zero in the public-covariate exact release proposition.
    s.Delta.sigma <- gs.sigma2.ridge(m.X, lambda, B.y)
    s.sigma.v     <- s.Delta.sigma / s.mu.v
    set.seed(seed + 2L)
    s.sigma2.dp <- max(s.sigma2.r + rnorm(1, 0, s.sigma.v), 1e-6)

    # Construct synthetic y: y_syn = X (I + lambda (X'X)^{-1}) beta_DP + e_null.
    # The lambda-corrected fitted part gives X'y_syn = (X'X + lambda I) beta_DP
    # (e_null is orthogonal to col(X)), so the analyst's ridge(lambda) refit
    # returns beta_DP exactly.
    m.XtX.inv <- solve(m.XtX)
    v.fitted <- as.numeric(m.X %*% (v.beta.dp + lambda * (m.XtX.inv %*% v.beta.dp)))
    m.H.ols <- m.X %*% m.XtX.inv %*% t(m.X)
    set.seed(seed + 3L)
    v.e.iid <- rnorm(n, 0, sqrt(s.sigma2.dp))
    v.e.proj <- as.numeric((diag(n) - m.H.ols) %*% v.e.iid)
    s.rss <- sum(v.e.proj^2)
    if(s.rss > 0)
        v.e.proj <- v.e.proj * sqrt((n - p) * s.sigma2.dp / s.rss)

    v.y.syn <- v.fitted + v.e.proj

    d.syn <- d
    d.syn[[y.nm]] <- v.y.syn

    # Not private, evaluation only: beta.ridge and sigma2.hat are the full-data
    # fits.  The DP release is d.syn.  LS, eps.beta and eps.sigma are gone with the
    # local-sensitivity channel that produced them.  The accounting is mu-GDP,
    # and mu.check re-composes the two channels in quadrature.
    list(d.syn       = d.syn,
         beta.ridge  = v.beta.ridge,
         beta.dp     = v.beta.dp,
         sigma2.hat  = s.sigma2.r,
         sigma2.dp   = s.sigma2.dp,
         Delta.beta  = s.Delta.beta,
         Delta.sigma = s.Delta.sigma,
         sigma.beta  = s.sigma.beta,
         sigma.v     = s.sigma.v,
         lambda      = lambda,
         eps         = eps,
         delta       = delta,
         mu.total    = s.mu.total,
         mu.beta     = s.mu.beta,
         mu.v        = s.mu.v,
         mu.check    = gdp.compose(c(s.mu.beta, s.mu.v)))
}

# Ridge regression fitter for verification.  With the lambda-corrected
# fitted component in gen.syn.dp.ridge.public, refitting ridge with the
# same lambda on the synthetic data returns beta_DP exactly (X public).
fit.ridge <- function(d, y.nm, x.nms, lambda) {
    fml <- reformulate(x.nms, response = y.nm)
    m.X <- model.matrix(fml, d)
    v.y <- d[[y.nm]]
    m.XtX <- crossprod(m.X)
    m.inv <- solve(m.XtX + lambda * diag(ncol(m.X)))
    as.numeric(m.inv %*% crossprod(m.X, v.y))
}

# +------------------------------------------------------------+
# |  Monte Carlo driver                                         |
# +------------------------------------------------------------+
# This file defines a mechanism.  The Monte Carlo demonstration lives in
# 36-mc-all-variants.r, the single DP-variant driver, which takes every bound
# from dp-public-bounds.r and refuses to run on an invented one.
