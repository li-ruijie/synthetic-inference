#!/usr/bin/env r
# +------------------------------------------------------------------------+
# | 09-symmetric-generator.r -- Direction-agnostic covariance calibration  |
# +------------------------------------------------------------------------+
#
# Purpose:
# Implements the direction-agnostic (symmetric) synthetic data generator
# from Section 5 of the paper. The analyst knows the model family (e.g.
# GLM with log link) but does not know which variable will serve as Y.
# The generator must therefore preserve inferential fidelity for every
# regression direction Z_j ~ g^{-1}(Z_{-j}' beta).
#
# Scenario 1: known model family, unknown dependent variable.
# Scenario 2: known model + known Y (reduces to 07-calibrated).
#
# Theoretical basis -- Theorem 6 (direction-agnostic fidelity):
# If the synthetic distribution P* satisfies three conditions within
# every stratum s:
#   (i)   P*(S = s) = P(S = s)            (stratum prevalence match)
#   (ii)  E*[W|S=s] = E[W|S=s]            (conditional mean match)
#   (iii) Var*(W|S=s) = Cov(W|S=s)        (conditional covariance match)
# where W = h(Z_cnt) is the vector of continuous variables on the
# h-scale (Definition 3: h = id for real, log for positive, logit for
# bounded), then for every GLM of the form Z_j ~ g^{-1}(Z_{-j}' beta),
# and any choice of dependent variable j, the synthetic data estimator
# is consistent up to a bias of order nu_3 / 3!, and exactly so when
# h_j = id (the identity-transform case below).
#
# Key insight: for a Gaussian copula on the h-scale, the conditional
# moments E[h(Z_j)|Z_{-j}] and Var[h(Z_j)|Z_{-j}] are both
# determined by the covariance matrix of h(Z). Matching Cov(h(Z))
# within each stratum suffices for direction-agnostic inferential
# fidelity to order k=2, regardless of which variable is Y.
#
# Gaussian copula symmetry: the copula models the joint
# distribution symmetrically -- no variable is privileged as Y or X.
# This is why covariance calibration on h(Z) yields
# direction-agnostic fidelity, unlike 07-calibrated which targets
# a specific response.
#
# Theorem 6 at h_j = id: when h_j = id for all j (i.e. all
# variables have real support), the bias bound is exactly zero:
# matching Cov(Z) within strata preserves every OLS regression.
#
# This generalises the OLS result (h=identity -> match Cov(Z) ->
# all OLS regressions preserved) to non-canonical GLMs via the
# support transform h (Definition 3).
#
# Algorithm overview:
# 1. Generate a base synthetic dataset via gen.syn.family() (Cornish-
#    Fisher quantile inversion per stratum, preserving marginals).
# 2. Apply calibrate.joint(): within each stratum s, compute an affine
#    transform on the h-scale that simultaneously matches the target
#    mean vector mu_s and covariance matrix Sigma_s.
#    The transform is: W*_new = (W* - mu_s*) L_s*^{-1} L_s + mu_s
#    where L_s = chol(Sigma_s). This matches both mean and covariance
#    in a single operation (see Theorem 6 conditions ii-iii).
# 3. ZI variables are handled separately (positive-part only) because
#    the structural zeros break the Gaussian copula assumption.
#
# Limitations:
# Remark (non-elliptical copulas) names the cases where the
# elliptical assumption behind the direction-agnostic theorem
# fails:
#   (a) asymmetric tail dependence (lambda_U != lambda_L),
#   (b) sector-dependent correlation structure,
#   (c) discrete responses.
# The Gaussian copula used here has lambda_U = lambda_L = 0, so
# (a) does not bite.  The zero-inflated and count scenarios are
# case (c).  Every h in use (id, log, logit) is monotone, which
# this implementation requires and that remark does not mention.
# the non-elliptical error bound (supp S9): the residual bias
# delta_s is controlled by the upper/lower tail dependence
# coefficients lambda_U, lambda_L of the copula. For the Gaussian
# copula, lambda_U = lambda_L = 0, so the error vanishes in the
# elliptical case.
#
# Dependencies:
# - 00-shared.r: environment setup, functional primitives
# - gen-syn-family.r: base generator (CF quantile), support detection,
#   h-transforms (Definition 3), EB shrinkage
# - d.sim.scenarios.RData: the 11 simulation scenarios for testing
#
# +------------------------------------------------------------------------+
# | Setup                                                                  |
# +------------------------------------------------------------------------+
source("lib/00-shared.r")
source("lib/gen-syn-family.r")
# +------------------------------------------------------------------------+
# | Joint calibration: multivariate affine transform on h-scale           |
# +------------------------------------------------------------------------+
#
# calibrate.joint() is the core of the symmetric generator. It enforces
# Theorem 6 conditions (ii) and (iii) by performing a per-stratum
# multivariate affine transform on the h-scale:
#
#   W*_new = (W* - mu_s*) L_s*^{-1} L_s_tilde + mu_s_tilde
#
# where:
#   W* = h(Z*_cnt) is the synthetic continuous block on the h-scale,
#   mu_s*, L_s* = mean vector and Cholesky factor of Cov(W*|S=s),
#   mu_s_tilde, L_s_tilde = EB-shrunk target mean and Cholesky factor.
#
# This simultaneously matches both the mean vector and the full
# covariance matrix within each stratum -- a single operation that
# satisfies both conditions of Theorem 6 for all regression directions.
#
# For ZI variables (zero-inflated, count): the structural zeros violate
# the continuity assumption of the Gaussian copula, so these are
# calibrated per-variable on the positive part only (mean + sd matching).
#
calibrate.joint <- function(d.syn, d.orig, supports = NULL,
                            shrink = TRUE) {

    # Classify variables: categorical (cat/bin) vs continuous (cnt).
    # The categorical variables define the strata S.  The continuous
    # variables are the W = h(Z_cnt) block to be jointly calibrated.
    v.types <- mapply(id.cat, d.orig)
    v.cat.nm <- names(v.types[v.types != "cnt"])
    v.cnt.nm <- names(v.types[v.types == "cnt"])

    # Need both categorical (for strata) and continuous (for calibration)
    if(length(v.cat.nm) == 0 || length(v.cnt.nm) == 0) return(d.syn)

    # Determine variable supports and build h-transforms (Definition 3).
    # h = id for "real", log for "positive", logit for "bounded".
    # detect.support() auto-classifies from data.  Explicit supports override.
    if(is.null(supports))
        supports <- `names<-`(mapply(detect.support, d.orig[v.cnt.nm]),
                              v.cnt.nm)
    l.tr <- Map(function(nm) make.transform(supports[[nm]], d.orig[[nm]]),
                v.cnt.nm)

    # Build strata from the categorical cross-classification.
    # Theorem 6 condition (i) -- P*(S=s) = P(S=s) -- is already satisfied
    # by the base generator gen.syn.family(), which bootstraps strata.
    v.strata.orig <- interaction(d.orig[v.cat.nm], drop = TRUE)
    v.strata.syn  <- interaction(d.syn[v.cat.nm], drop = TRUE)
    l.orig.rows <- split(seq_len(nrow(d.orig)), v.strata.orig)
    l.syn.rows  <- split(seq_len(nrow(d.syn)), v.strata.syn)

    # Separate ZI/count variables (per-variable calibration) from the
    # non-ZI continuous block (joint multivariate calibration).
    # ZI variables have structural zeros that break the Gaussian copula
    # assumption, so they cannot participate in the joint affine transform.
    v.zi.nm    <- names(supports[supports %in% c("zeroinfl", "count")])
    v.joint.nm <- setdiff(v.cnt.nm, v.zi.nm)
    s.p        <- length(v.joint.nm)

    # +--------------------------------------------------------------------+
    # | Part 1: joint calibration for non-ZI continuous variables         |
    # +--------------------------------------------------------------------+
    # This block implements Theorem 6 conditions (ii)-(iii) for the
    # non-ZI continuous variables via a multivariate affine transform.
    if(s.p >= 1) {

        # Compute per-stratum h-scale target moments from the original data.
        # For each stratum s, we need:
        #   mu_s    = E[h(Z_cnt) | S=s]   -- the p-vector of means
        #   Sigma_s = Cov(h(Z_cnt) | S=s) -- the p x p covariance matrix
        # These are the targets that the affine transform will match.
        # Strata with fewer than p+1 complete cases are skipped (NULL).
        l.strata.moments <- Map(function(s.key) {
            v.oi <- l.orig.rows[[s.key]]
            # Apply h-transform to each variable within this stratum
            m.h <- mapply(function(nm) {
                v.h <- l.tr[[nm]][["h"]](d.orig[[nm]][v.oi])
                v.h[!is.finite(v.h)] <- NA
                v.h
            }, v.joint.nm)
            if(!is.matrix(m.h)) m.h <- matrix(m.h, ncol = 1)
            m.h <- `colnames<-`(m.h, v.joint.nm)
            b.cc <- complete.cases(m.h)
            s.ncc <- sum(b.cc)
            # Need at least p+1 observations for a non-degenerate covariance
            if(s.ncc < s.p + 1) return(NULL)
            m.h <- m.h[b.cc, , drop = FALSE]
            list(mu    = colMeans(m.h),
                 sigma = `if`(s.p > 1, cov(m.h), matrix(var(m.h), 1, 1)),
                 n     = s.ncc)
        }, names(l.orig.rows))

        # +----------------------------------------------------------------+
        # | Empirical Bayes shrinkage for stratum mean vectors and         |
        # | covariance matrices                                            |
        # +----------------------------------------------------------------+
        # Sparse strata (small n_s) yield noisy estimates of mu_s and
        # Sigma_s. EB shrinkage borrows strength across strata to
        # stabilise both targets. This is the multivariate extension of
        # the scalar shrinkage in 07-calibrated (shrink.moments).
        #
        # Mean shrinkage (per-component normal-normal EB):
        #   mu_s_tilde_j = w_sj * mu_s_hat_j + (1 - w_sj) * mu_grand_j
        #   w_sj = tau^2_j / (tau^2_j + sigma^2_jj / n_s)
        # where tau^2_j = between-stratum variance of mu_s_hat_j minus
        # the expected sampling noise.
        #
        # Covariance shrinkage (toward MANOVA-style pooled Sigma):
        #   Sigma_s_tilde = w_cov * Sigma_s_hat + (1-w_cov) * Sigma_pool
        #   w_cov = (n_s - 1) / (n_s + 1)
        # This ensures every stratum has a positive-definite target
        # covariance (Sigma_pool is PD), which is needed for the
        # Cholesky decomposition in the affine transform.
        v.valid <- !mapply(is.null, l.strata.moments)
        if(shrink && sum(v.valid) >= 3) {
            # Collect all stratum means into a matrix (K x p) and sizes
            m.mu.all <- Reduce(rbind, Map(function(l) l[["mu"]],
                                          l.strata.moments[v.valid]))
            v.ns.all <- mapply(function(l) l[["n"]],
                               l.strata.moments[v.valid])
            v.grand  <- colMeans(m.mu.all)

            # Pooled within-stratum covariance (MANOVA-style):
            # Sigma_pool = sum_s (n_s - 1) Sigma_s / sum_s (n_s - 1)
            s.df.total <- sum(v.ns.all) - sum(v.valid)
            m.sigma.pool <- Reduce(`+`, Map(function(l) {
                (l[["n"]] - 1) * l[["sigma"]]
            }, l.strata.moments[v.valid])) / s.df.total

            Map(function(s.key) {
                l <- l.strata.moments[[s.key]]
                if(is.null(l)) return(NULL)
                s.ns <- l[["n"]]

                # Per-component EB shrinkage for the mean vector.
                # tau^2_j = max(0, Var(mu_hat_j across strata) -
                #                  E[sigma^2_jj / n_s])
                # This is method-of-moments estimation of the between-
                # stratum variance, subtracting the expected sampling noise.
                v.var.diag <- diag(as.matrix(l[["sigma"]]))
                v.tau2 <- pmax(0, mapply(function(j) {
                    var(m.mu.all[, j]) -
                        mean(v.var.diag[j] / v.ns.all)
                }, seq_len(s.p)))
                # w_sj = tau^2_j / (tau^2_j + sigma^2_jj / n_s)
                v.w <- mapply(function(tau2, var.j) {
                    s.d <- tau2 + var.j / s.ns
                    `if`(s.d > 0, tau2 / s.d, 1)
                }, v.tau2, v.var.diag)
                l.strata.moments[[s.key]][["mu"]] <<-
                    v.w * l[["mu"]] + (1 - v.w) * v.grand

                # Shrink covariance toward pooled: w_cov = (n_s-1)/(n_s+1).
                # As n_s -> infinity, w_cov -> 1 (no shrinkage).
                # As n_s -> 1, w_cov -> 0 (fully pooled).
                s.w.cov <- (s.ns - 1) / (s.ns + 1)
                l.strata.moments[[s.key]][["sigma"]] <<-
                    s.w.cov * l[["sigma"]] +
                    (1 - s.w.cov) * m.sigma.pool
            }, names(l.strata.moments)[v.valid])
        }

        # +----------------------------------------------------------------+
        # | Multivariate affine transform per stratum                      |
        # +----------------------------------------------------------------+
        # This is the key operation that enforces Theorem 6. For each
        # stratum s, we transform the synthetic h-scale data W* to match
        # the target mean mu_s_tilde and covariance Sigma_s_tilde:
        #
        #   W*_new = (W* - mu_s*) L_s*^{-1} L_s_tilde + mu_s_tilde
        #
        # where L = chol(Sigma) is the upper Cholesky factor.
        #
        # Why this works:
        # - Centring by mu_s* and left-multiplying by L_s*^{-1} produces
        #   uncorrelated, unit-variance residuals (whitening).
        # - Left-multiplying by L_s_tilde re-introduces the target
        #   covariance structure.
        # - Adding mu_s_tilde shifts the mean to the target.
        #
        # Result: E[W*_new|S=s] = mu_s_tilde, Cov(W*_new|S=s) = Sigma_s_tilde.
        # Both conditions (ii) and (iii) of Theorem 6 are satisfied in
        # a single affine operation.
        #
        # Theorem 6 at h = id: when all supports are "real" (h=id),
        # this exactly preserves every OLS regression coefficient because
        # the affine transform is a bijection on the original scale.
        Map(function(s.key) {
            l.target <- l.strata.moments[[s.key]]
            if(is.null(l.target)) return(NULL)
            if(!s.key %in% names(l.syn.rows)) return(NULL)
            v.si <- l.syn.rows[[s.key]]
            if(length(v.si) < s.p + 1) return(NULL)

            # Compute current synthetic h-scale moments
            m.h.syn <- mapply(function(nm)
                l.tr[[nm]][["h"]](d.syn[[nm]][v.si]), v.joint.nm)
            if(!is.matrix(m.h.syn))
                m.h.syn <- matrix(m.h.syn, ncol = 1)

            v.mu.syn    <- colMeans(m.h.syn)
            m.sigma.syn <- `if`(s.p > 1, cov(m.h.syn),
                                matrix(var(m.h.syn), 1, 1))

            # Cholesky decomposition of both target and synthetic covariance.
            # If either is not positive definite (e.g. degenerate stratum),
            # skip this stratum -- the base CF generator's moments remain.
            b.ok <- tryCatch({
                m.L.target <- chol(l.target[["sigma"]])
                m.L.syn    <- chol(m.sigma.syn)
                TRUE
            }, error = function(e) FALSE)
            if(!b.ok) return(NULL)

            # The affine transform: whiten, then re-colour + shift.
            # sweep(..., 2, v.mu.syn) centres each column (subtract mean).
            # %*% solve(m.L.syn) whitens (decorrelate + unit variance).
            # %*% m.L.target imposes the target covariance structure.
            # sweep(..., 2, -mu_tilde) shifts to the target mean.
            m.h.new <- sweep(m.h.syn, 2, v.mu.syn) %*%
                       solve(m.L.syn) %*% m.L.target
            m.h.new <- sweep(m.h.new, 2, -l.target[["mu"]])

            # Apply h^{-1} to return to the original scale and write back.
            # h^{-1} is: id for real, exp for positive, plogis for bounded.
            Map(function(nm, j) {
                d.syn[[nm]][v.si] <<- l.tr[[nm]][["hinv"]](m.h.new[, j])
            }, v.joint.nm, seq_len(s.p))
        }, names(l.orig.rows))
    }

    # +--------------------------------------------------------------------+
    # | Part 2: per-variable calibration for ZI/count variables            |
    # +--------------------------------------------------------------------+
    # Zero-inflated and count variables have structural zeros (point mass
    # at 0) that prevent them from participating in the joint Gaussian
    # copula. Instead, they are calibrated univariately on the positive
    # Part only: the zero/non-zero indicator was already matched by the
    # categorical bootstrap (it appears as a binary stratum variable),
    # and the positive-part mean and SD are matched via location-scale
    # adjustment on the h-scale (h = log for ZI, log(x+0.5) for count).
    #
    # This is a weaker form of calibration than the joint transform --
    # it only matches the first two marginal moments, not the full
    # covariance -- but ZI variables inherently break the multivariate
    # normality assumption anyway (the non-elliptical error bound of supp S9 applies).
    if(length(v.zi.nm) > 0) {

        # Pre-compute EB-shrunk targets for each ZI variable.
        # Same scalar EB approach as 07-calibrated's calibrate():
        #   mu_s_tilde = w_s * mu_s_hat + (1-w_s) * mu_grand
        #   sd_s_tilde = sqrt(w_var * sd_s^2 + (1-w_var) * sd_pool^2)
        l.zi.targets <- Map(function(nm) {
            h    <- l.tr[[nm]][["h"]]
            # Initialise per-stratum raw moments (NA = missing)
            v.mu.raw <- v.sd.raw <- v.ns <- `names<-`(
                rep(NA_real_, length(l.orig.rows)), names(l.orig.rows))
            # Compute h-scale mean and SD for positive part only
            Map(function(s.key) {
                v.raw <- d.orig[[nm]][l.orig.rows[[s.key]]]
                v.raw <- v.raw[v.raw > 0]
                v.h <- h(v.raw)
                v.h <- v.h[is.finite(v.h)]
                v.ns[[s.key]]     <<- length(v.h)
                if(length(v.h) >= 1) v.mu.raw[[s.key]] <<- mean(v.h)
                if(length(v.h) >= 2) v.sd.raw[[s.key]] <<- sd(v.h)
            }, names(l.orig.rows))
            if(!shrink) return(list(mu = v.mu.raw, sd = v.sd.raw))

            # EB shrinkage for means (normal-normal model)
            b.ok <- !is.na(v.mu.raw) & !is.na(v.sd.raw) & v.ns >= 2
            v.mu.sh <- v.mu.raw
            if(sum(b.ok) >= 3) {
                s.grand  <- mean(v.mu.raw[b.ok])
                s.sigma2 <- mean(v.sd.raw[b.ok]^2)
                s.tau2   <- max(0, var(v.mu.raw[b.ok]) -
                             s.sigma2 * mean(1 / v.ns[b.ok]))
                v.w <- mapply(function(ns, ok) {
                    if(!ok) return(1)
                    s.d <- s.tau2 + s.sigma2 / ns
                    `if`(s.d > 0, s.tau2 / s.d, 1)
                }, v.ns, b.ok)
                v.mu.sh <- ifelse(b.ok,
                    v.w * v.mu.raw + (1 - v.w) * s.grand, v.mu.raw)
            }

            # EB shrinkage for SDs (toward pooled variance)
            b.ok.sd <- !is.na(v.sd.raw) & v.ns >= 2
            v.sd.sh <- v.sd.raw
            if(sum(b.ok.sd) >= 3) {
                s.var.pool <- mean(v.sd.raw[b.ok.sd]^2)
                v.w.var <- ifelse(b.ok.sd,
                                  (v.ns - 1) / (v.ns + 1), 0)
                v.sd.sh <- ifelse(b.ok.sd,
                    sqrt(pmax(0, v.w.var * v.sd.raw^2 +
                                 (1 - v.w.var) * s.var.pool)),
                    sqrt(s.var.pool))
                v.sd.sh[is.na(v.sd.raw)] <- NA_real_
            }
            list(mu = v.mu.sh, sd = v.sd.sh)
        }, v.zi.nm)
        l.zi.targets <- `names<-`(l.zi.targets, v.zi.nm)

        # Per-variable, per-stratum calibration (positive part only).
        # For each ZI variable in each stratum:
        #   1. Shift h(positive part) to match the target mean.
        #   2. Scale residuals to match the target SD.
        #   3. Apply h^{-1} to return to the original scale.
        # The structural zeros are left untouched.
        Map(function(nm) {
            h    <- l.tr[[nm]][["h"]]
            hinv <- l.tr[[nm]][["hinv"]]
            Map(function(s.key) {
                if(!s.key %in% names(l.syn.rows))  return(NULL)
                v.oi <- l.orig.rows[[s.key]]
                v.si <- l.syn.rows[[s.key]]
                if(length(v.si) < 2 || length(v.oi) < 2) return(NULL)

                v.syn.raw <- d.syn[[nm]][v.si]
                b.pos.s <- v.syn.raw > 0
                if(sum(b.pos.s) < 2) return(NULL)

                v.syn.h <- h(v.syn.raw[b.pos.s])
                v.syn.h <- v.syn.h[is.finite(v.syn.h)]
                if(length(v.syn.h) < 2) return(NULL)

                # Step 1: shift to match target mean
                s.mu.o <- l.zi.targets[[nm]][["mu"]][[s.key]]
                if(is.na(s.mu.o)) return(NULL)
                s.mu.s <- mean(v.syn.h)
                v.syn.h <- v.syn.h + (s.mu.o - s.mu.s)

                # Step 2: scale residuals to match target SD
                s.sd.o <- l.zi.targets[[nm]][["sd"]][[s.key]]
                if(!is.na(s.sd.o)) {
                    s.sd.s <- sd(v.syn.h)
                    if(s.sd.s > .Machine[["double.eps"]])
                        v.syn.h <- s.mu.o + (v.syn.h - mean(v.syn.h)) *
                                   (s.sd.o / s.sd.s)
                }

                # Step 3: return to original scale via h^{-1}
                v.syn.raw[b.pos.s] <- hinv(v.syn.h)
                d.syn[[nm]][v.si] <<- v.syn.raw
            }, names(l.orig.rows))
        }, v.zi.nm)
    }

    d.syn
}
# +------------------------------------------------------------------------+
# | Public API: gen.syn.symmetric()                                        |
# +------------------------------------------------------------------------+
# The main entry point. Two-stage pipeline:
#   Stage 1 (gen.syn.family): generate a base synthetic dataset using
#     stratified CF quantile inversion. This produces marginals that match
#     per-stratum mean and higher cumulants, but the joint covariance
#     structure is only approximately preserved (via pooled within-stratum
#     correlation, not per-stratum covariance).
#   Stage 2 (calibrate.joint): apply the multivariate affine transform
#     from Theorem 6 to exactly match per-stratum mean vectors and
#     covariance matrices on the h-scale.
#
# After both stages, the synthetic data satisfies all three conditions
# of Theorem 6, yielding direction-agnostic inferential fidelity for
# every GLM Z_j ~ g^{-1}(Z_{-j}' beta).
#
gen.syn.symmetric <- function(d, order.cf = 4, supports = NULL,
                              shrink = TRUE, seed = 6489) {
    # Stage 1: base CF generator (marginal fidelity)
    d.syn <- gen.syn.family(d, order = order.cf, supports = supports,
                            shrink = shrink, seed = seed)
    # Stage 2: joint covariance calibration (Theorem 6)
    calibrate.joint(d.syn, d, supports = supports, shrink = shrink)
}
# +------------------------------------------------------------------------+
# | Test: all-direction OLS comparison                                     |
# +------------------------------------------------------------------------+
# Guard: skip tests when this file is sourced by the Sobol sensitivity
# analysis (18-mc-projection.r sets .sobol.main before sourcing).
#
# These tests empirically verify the theoretical predictions:
#
# Test 1 -- Scenario 1 (Gaussian/identity, all real support):
#   Theorem 6 at h = id: h = id for all variables, so the affine
#   transform on h-scale is a transform on the original scale. The bias
#   bound is exactly zero. We verify by regressing each continuous
#   variable on the rest and comparing OLS coefficients.
#   Expected: "symm" column shows near-zero % deviation for all directions.
#
# Test 2 -- Scenario 4 (ZI quasi-gamma):
#   the non-elliptical error bound (supp S9): ZI variables violate
#   the elliptical assumption. The residual bias delta_s is controlled
#   by tail dependence (lambda_U, lambda_L), which is zero for the
#   Gaussian copula. We test in three directions: response (ZI), age
#   (real), income (positive). The symmetric generator should improve
#   on the uncalibrated base in all directions, not just one.
#
# Comparison:
#   - "orig"  = OLS on original data (ground truth)
#   - "uncal" = OLS on base CF generator (gen.syn.family, no joint cal.)
#   - "symm"  = OLS on symmetric generator (gen.syn.symmetric)
#   The % deviation shows (beta_syn - beta_orig) / |beta_orig| * 100.
#
if(!exists(".sobol.main")) {
load("results/d.sim.scenarios.RData")

# Support specifications for simulation scenarios 1 to 10 (Definition 3).
# Each named vector maps variable -> support type:
#   "real"     -> h = identity
#   "positive" -> h = log
#   "zeroinfl" -> h = log (positive part only, zero indicator categorical)
#   "count"    -> h = log(x + 0.5)
#   "bounded"  -> h = logit (scaled to [0,1])
l.supports <- list(
    c(response = "real",     age = "real", income = "real"),
    c(response = "real",     age = "real", income = "positive"),
    c(response = "positive", age = "real", income = "positive"),
    c(response = "zeroinfl", age = "real", income = "positive"),
    c(response = "zeroinfl", age = "real", income = "positive"),
    c(response = "count",    age = "real", income = "positive"),
    c(response = "positive", age = "real", income = "positive"),
    c(response = "zeroinfl", age = "real", income = "positive"),
    c(response = "bounded",  age = "real", income = "positive"),
    c(response = "real",     age = "real", income = "real")
)

# Helper: fit OLS with y.nm as response, all other standard variables
# as predictors. Returns NULL on failure (e.g. singular design).
fit.ols <- function(d, y.nm) {
    v.pred <- setdiff(intersect(
        c("response", "age", "female", "race", "income"), colnames(d)),
        y.nm)
    fml <- reformulate(v.pred, response = y.nm)
    tryCatch(lm(fml, d), error = function(e) NULL)
}

# -- Test 1: Scenario 1 (Gaussian/identity, all real) ---------------------
# Theorem 6 at h = id => bias bound = 0 => all OLS regressions preserved.
cat("=== Direction-agnostic test: Scenario 1 (Gaussian/identity) ===\n")
cat("Any continuous variable as Y, OLS in all directions\n\n")

d.orig <- l.scenarios[[1]][["data"]]
sup <- l.supports[[1]]
d.sym <- gen.syn.symmetric(d.orig, order.cf = 4, supports = sup)
d.cal <- gen.syn.family(d.orig, order = 4, supports = sup)

# Test all three regression directions: Y = response, age, income
v.cnt.test <- c("response", "age", "income")
Map(function(y.nm) {
    mo.o <- fit.ols(d.orig, y.nm)
    mo.s <- fit.ols(d.sym, y.nm)
    mo.u <- fit.ols(d.cal, y.nm)
    if(is.null(mo.o) || is.null(mo.s) || is.null(mo.u)) return(NULL)
    v.o <- coef(mo.o)
    v.s <- coef(mo.s)
    v.u <- coef(mo.u)
    v.pct.s <- round(100 * (v.s - v.o) / abs(v.o), 1)
    v.pct.u <- round(100 * (v.u - v.o) / abs(v.o), 1)
    cat(sprintf("--- Y = %s ---\n", y.nm))
    d.out <- data.frame(
        orig   = sprintf("%8.4f", v.o),
        uncal  = sprintf("%8.4f (%+.1f%%)", v.u, v.pct.u),
        symm   = sprintf("%8.4f (%+.1f%%)", v.s, v.pct.s))
    print(d.out)
    cat("\n")
}, v.cnt.test)

# -- Test 2: Scenario 4 (ZI quasi-gamma) ---------------------------------
# the non-elliptical error bound (supp S9): ZI breaks elliptical assumption, but Gaussian copula
# has lambda_U = lambda_L = 0, so error bound delta_s is small.
# Remark (non-elliptical copulas): this scenario exercises the
# Discrete-response case, the zero mass being the discreteness.
# The tail-dependence case does not apply at lambda_U = lambda_L = 0.
cat("=== Direction-agnostic test: Scenario 4 (ZI quasi-gamma) ===\n\n")
d.orig <- l.scenarios[[4]][["data"]]
sup <- l.supports[[4]]
d.sym <- gen.syn.symmetric(d.orig, order.cf = 4, supports = sup)
d.cal <- gen.syn.family(d.orig, order = 4, supports = sup)

# Test three directions: response (ZI), age (real), income (positive)
Map(function(y.nm) {
    mo.o <- fit.ols(d.orig, y.nm)
    mo.s <- fit.ols(d.sym, y.nm)
    mo.u <- fit.ols(d.cal, y.nm)
    if(is.null(mo.o) || is.null(mo.s) || is.null(mo.u)) return(NULL)
    v.o <- coef(mo.o)
    v.s <- coef(mo.s)
    v.u <- coef(mo.u)
    v.pct.s <- round(100 * (v.s - v.o) / abs(v.o), 1)
    v.pct.u <- round(100 * (v.u - v.o) / abs(v.o), 1)
    cat(sprintf("--- Y = %s ---\n", y.nm))
    d.out <- data.frame(
        orig   = sprintf("%8.4f", v.o),
        uncal  = sprintf("%8.4f (%+.1f%%)", v.u, v.pct.u),
        symm   = sprintf("%8.4f (%+.1f%%)", v.s, v.pct.s))
    print(d.out)
    cat("\n")
}, c("response", "age", "income"))
} # end if(!exists(".sobol.main"))
# +--------------------------------------------------------+
