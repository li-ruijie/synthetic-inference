#!/usr/bin/env r
# +------------------------------------------------------------+
# |  Full-DP variant: public X + DP synthesised y (OLS only)   |
# +------------------------------------------------------------+
#
# When X is treated as public (common in health/regulatory contexts
# where covariates are demographic/non-sensitive), we can construct a
# Fully (eps, delta)-DP synthetic dataset for OLS with:
#
#   1. beta_DP   = beta_hat + Gaussian-mechanism noise, mu_beta-GDP at
#                  the global sensitivity Delta_beta.
#   2. sigma2_DP = sigma2   + Gaussian-mechanism noise, mu_v-GDP at the
#                  global sensitivity Delta_sigma2 = 4 B_y^2 K_X/(n-p),
#                  K_X = max_i max{1-h_ii, sum_{j!=i}|h_ij|} (public X).
#   3. e_syn = randomly drawn Gaussian residuals (no data used)
#   4. e_syn_proj = (I - H) e_syn, scaled to variance sigma2_DP
#   5. y_syn = X beta_DP + e_syn_proj
#   6. Release (X_orig, y_syn)
#
# Total privacy: mu_total-GDP with mu_total = sqrt(mu_beta^2 + mu_v^2)
# (mu-quadrature composition, Dong-Roth-Su 2022), equivalent to the input
# (eps, delta) by the exact GDP dual.  X is public (no DP spent).  The
# analyst fits OLS(X_orig, y_syn) and obtains beta_DP exactly by the
# projection construction, with SE from e_syn_proj scaled to sigma2_DP.
#
# This variant eliminates the copula-based leak of y through X that
# limits the gen.syn.dp.projected() scope.  It is limited to OLS
# (identity link) in this implementation.  Log-link extension follows
# the same pattern but uses the linearised score equation.
#
# Dependencies: base R, plus dp.gm.multiplier() (28-dp-projection.r).
if(!exists("dp.gm.multiplier")) {
    .dp.suppress.demo <- TRUE; source("lib/28-dp-projection.r")
}

# +------------------------------------------------------------+
# |  Full-DP OLS synthesiser (public X)                        |
# +------------------------------------------------------------+
# Arguments:
#   d         original data frame
#   y.nm      response column
#   x.nms     predictor names (if NULL, all non-y columns)
#   eps       total privacy budget epsilon
#   delta     privacy budget delta
#   B.y       public bound on |y|
#   eps.split numeric in (0, 1): fraction of the mu^2 (GDP) budget
#             allocated to beta (1 - eps.split to sigma^2)
#   seed      RNG seed
# Returns: list with d.syn, beta.dp, sigma2.dp, the GDP budgets mu.total,
#          mu.beta, mu.v, and the nominal eps shares eps.beta, eps.sigma.
gen.syn.dp.ols.public <- function(d, y.nm = "response", x.nms = NULL,
                                    eps = 1, delta = 1e-6,
                                    B.y = NULL,
                                    eps.split = 0.8, seed = 6489) {
    if(is.null(B.y)) stop("B.y (public y bound) required")
    if(is.null(x.nms)) x.nms <- setdiff(names(d), y.nm)
    fml <- reformulate(x.nms, response = y.nm)
    m.X <- model.matrix(fml, d)
    n <- nrow(m.X); p <- ncol(m.X)

    mo <- lm(fml, d)
    v.beta <- coef(mo)
    s.sigma2 <- summary(mo)[["sigma"]]^2

    # Allocate the mu-GDP budget: mu_total from (eps, delta), split the mu^2
    # budget between beta and sigma^2 by eps.split.
    s.mu.total <- gdp.mu.from.eps.delta(eps, delta)
    s.mu.beta  <- s.mu.total * sqrt(eps.split)
    s.mu.v     <- s.mu.total * sqrt(1 - eps.split)

    # Step 1: DP beta via the Gaussian mechanism at the public-X global
    # sensitivity Delta_beta = 2 B.y max_i ||(X'X)^{-1} x_i|| (was local DFBETA).
    s.Delta.beta <- gs.beta.ols.public(m.X, B.y)
    s.sigma.beta <- s.Delta.beta / s.mu.beta
    set.seed(seed + 1L)
    v.beta.dp <- v.beta + rnorm(p, 0, s.sigma.beta)
    names(v.beta.dp) <- names(v.beta)

    # Step 2: DP sigma^2 via the Gaussian mechanism (Gaussianised, was Laplace)
    # at the exact global sensitivity Delta_v = 4 B.y^2 K_X / (n - p), K_X the
    # hat-row design constant computed from the public X (see gs.sigma2).  The
    # 1e-6 floor is the numerical form of the truncation at zero in the
    # proposition.
    s.Delta.v  <- gs.sigma2(m.X, B.y)
    s.sigma.v  <- s.Delta.v / s.mu.v
    set.seed(seed + 2L)
    s.sigma2.dp <- max(s.sigma2 + rnorm(1, 0, s.sigma.v), 1e-6)
    s.sigma.dp  <- sqrt(s.sigma2.dp)

    # Step 3: draw iid Gaussian residuals (no data dependence)
    set.seed(seed + 3L)
    v.e.iid <- rnorm(n, 0, s.sigma.dp)

    # Step 4: project onto (I - H) so the analyst's OLS refit recovers beta.dp
    m.H <- m.X %*% solve(crossprod(m.X)) %*% t(m.X)
    v.e.proj <- as.numeric((diag(n) - m.H) %*% v.e.iid)
    s.rss <- sum(v.e.proj^2)
    if(s.rss > 0)
        v.e.proj <- v.e.proj * sqrt((n - p) * s.sigma2.dp / s.rss)

    # Step 5: y_syn = X beta_DP + e_syn_proj
    v.y.syn <- as.numeric(m.X %*% v.beta.dp) + v.e.proj

    d.syn <- d
    d.syn[[y.nm]] <- v.y.syn

    list(d.syn      = d.syn,
         beta.hat   = v.beta,
         beta.dp    = v.beta.dp,
         sigma2.hat = s.sigma2,
         sigma2.dp  = s.sigma2.dp,
         mu.total   = s.mu.total,
         mu.beta    = s.mu.beta,
         mu.v       = s.mu.v,
         eps.beta   = eps * eps.split,
         eps.sigma  = eps * (1 - eps.split),
         sigma.beta = s.sigma.beta,
         sigma.v    = s.sigma.v,
         B.y        = B.y)
}

# Laplace random generator (base R has no rlaplace).
rlaplace <- function(n, location = 0, scale = 1) {
    v.u <- runif(n) - 0.5
    location - scale * sign(v.u) * log(1 - 2 * abs(v.u))
}

# +------------------------------------------------------------+
# |  Monte Carlo driver                                         |
# +------------------------------------------------------------+
# This file defines a mechanism.  The Monte Carlo demonstration lives in
# 36-mc-all-variants.r, the single DP-variant driver, which takes every bound
# from dp-public-bounds.r and refuses to run on an invented one.
