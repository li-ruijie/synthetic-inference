#!/usr/bin/env r
# gen.syn.family: support-aware synthetic data generation
#
# Handles: real, positive, zero-inflated, count, bounded
# Uses: stratified categorical bootstrap + CF quantile on transformed scale
# Features: empirical Bayes shrinkage for sparse strata,
#           KNN conditioning for continuous predictors
#
# source() this file to get gen.syn.family()
#
# +------------------------------------------------------------------------+
# | Overview: connection to the paper                                      |
# +------------------------------------------------------------------------+
#
# This file implements the copula-based synthetic data generator that
# enters the paper in two roles (Section 1):
#   (1) constructive proof that the moment-matching conditions of
#       Theorem 2 (Partial moment matching) are achievable, and
#   (2) a generator that is competitive in practice and that,
#       combined with sufficient statistics projection (Theorem 7),
#       provides exact inferential fidelity.
#
# The algorithm proceeds in five stages:
#   A. Detect variable types and support.  Build h / h^{-1} transforms
#      (Definition 3, support transform).
#   B. Stratify by categorical predictors.  Compute per-stratum
#      conditional cumulants on the h-scale (Theorem 2, bias bound).
#   C. Apply empirical Bayes shrinkage to stratum means
#      (Theorem 5, sparse strata MSE and the optimal shrinkage rate, supplement S8).
#   D. Estimate the pooled within-stratum correlation matrix on the
#      h-scale, then generate correlated uniforms via a Gaussian
#      copula (Sklar's theorem [sklar-fonction-1959], and Nelsen (2007)
#      [nelsen-introduc-2007]).  This targets Cov(h(Z) | S=s) as
#      required by Theorem 6 (direction-agnostic fidelity).
#   E. Map uniforms -> CF quantile -> h^{-1} to obtain synthetic
#      observations on the original scale (Cornish & Fisher (1938)
#      [cornish-moments-1938].  Hall (1992) [hall-bootstra-1992]).

# +------------------------------------------------------------------------+
# | Helpers: variable-type detection                                       |
# +------------------------------------------------------------------------+
#
# id.cat classifies a variable as categorical ("cat"), binary ("bin"),
# or continuous ("cnt").  The distinction matters because categorical
# variables define the stratification (the "S" in Theorem 6,
# direction-agnostic fidelity), while continuous variables are
# synthesised via the CF quantile expansion on the h-scale.
id.cat <- function(v, threshold = 10) {
    b.cat <- "cnt"
    v.ulen <- length(unique(v))
    if(is.factor(v) | is.character(v) | is.logical(v))
        b.cat <- "cat"
    if(v.ulen <= threshold)
        b.cat <- "cat"
    if((b.cat == "cat") && (v.ulen == 2))
        b.cat <- "bin"
    b.cat
}

# detect.support infers the support of a continuous variable so that
# we can select the correct support transform h (Definition 3).
# The canonical choices mirror the paper's table:
#   R          -> identity
#   (0, inf)   -> log
#   zero-inflated (spike at 0 plus positive mass) -> log on the
#              positive part, with a separate binary indicator for
#              the zero/non-zero split.
detect.support <- function(v) {
    if(all(v > 0))                                         return("positive")
    if(any(v == 0) && all(v >= 0) && mean(v == 0) > 0.01) return("zeroinfl")
    "real"
}

# +------------------------------------------------------------------------+
# | Support transforms h and h^{-1} (Definition 3)                        |
# +------------------------------------------------------------------------+
#
# Definition 3 (support transform): h is a measurable injection from
# the support Y into R with a measurable left inverse h^{-1} defined
# on all of R (h^{-1}(h(y)) = y), chosen so that the transformed
# variable Z = h(Y) lives on the real line.  On the continuous
# supports h is a bijection onto R.  All moment matching
# (Theorem 2) and the CF expansion operate on the h-scale.
#
# The canonical choices are:
#   Y = R          : h = id,    h^{-1} = id
#   Y = (0, inf)   : h = log,   h^{-1} = exp
#   Y = (0,1)      : h = logit, h^{-1} = expit
#   Y = N_0 (count): h = log(y + 1/2), h^{-1} = round(exp(z) - 1/2)
#
# The "bounded" case generalises the (0,1) logit to an arbitrary
# [lo, hi] interval by first rescaling to (0,1) with boundary
# padding eps = 1/(2n).
make.transform <- function(support, v = NULL) {
    switch(support,
        real = list(
            h    = identity,
            hinv = identity),
        positive = list(
            h    = log,
            hinv = exp),
        zeroinfl = list(
            h    = log,
            hinv = exp),
        count = {
            # Shift by 0.5 before log, matching Definition 3's
            # count transform: h(y) = log(y + 1/2).  The inverse
            # rounds to the nearest non-negative integer.
            s.shift <- 0.5
            list(h    = function(x) log(x + s.shift),
                 hinv = function(x) pmax(round(exp(x) - s.shift), 0))
        },
        bounded = {
            # Generalised logit for bounded support [lo, hi].
            # Boundary padding eps prevents logit(0) = -Inf and
            # logit(1) = +Inf, keeping Z = h(Y) finite.
            s.lo  <- min(v, na.rm = TRUE)
            s.hi  <- max(v, na.rm = TRUE)
            s.eps <- 1 / (2 * length(v))
            list(h = function(x) {
                     x.sc <- (x - s.lo) / (s.hi - s.lo)
                     x.sc <- pmin(pmax(x.sc, s.eps), 1 - s.eps)
                     qlogis(x.sc)
                 },
                 hinv = function(x) {
                     s.lo + (s.hi - s.lo) * plogis(x)
                 })
        },
        stop(paste("Unknown support:", support)))
}

# +------------------------------------------------------------------------+
# | Cornish-Fisher quantile expansion (Theorem 2 mechanism)                |
# +------------------------------------------------------------------------+
#
# The CF expansion [cornish-moments-1938] inverts the Edgeworth
# expansion [hall-bootstra-1992] of the CDF to express quantiles in
# terms of cumulants.  This is the mechanism by which conditional
# moments on the h-scale are embedded into the synthetic marginals.
#
# Theorem 2 (partial moment matching) proves that matching k
# conditional moments of Z = h(Y) yields a score bias bounded by
# O(kappa_{k+1} / (k+1)!), where kappa_{k+1} is the unmatched
# (k+1)-th conditional cumulant.  The CF expansion to order k
# achieves exactly this: it reproduces the first k cumulants by
# construction and truncates higher-order terms.
#
# Proposition 1 (summary statistics) further shows that k = 2 on
# the h-scale suffices to control the full GLM output (coefficients,
# dispersion, standard errors, test statistics).

# std.cumulants: extract standardised cumulants gamma_3, ..., gamma_k
# from a sample.  These are the shape parameters that enter the CF
# expansion beyond the Gaussian baseline (mu, sigma).
#
# The conversion from raw central moments mu_r to cumulants gamma_r
# uses the standard identities:
#   gamma_1 = mu_3                           (skewness)
#   gamma_2 = mu_4 - 3                       (excess kurtosis)
#   gamma_3 = mu_5 - 10 * mu_3              (5th cumulant)
#   gamma_4 = mu_6 - 15 * mu_4 - 10*mu_3^2 + 30  (6th cumulant)
std.cumulants <- function(v, order) {
    s.n <- length(v)
    s.sigma <- sd(v)
    if(s.sigma < .Machine[["double.eps"]] || s.n < order + 1)
        return(rep(0, max(order - 2, 0)))
    v.c <- (v - mean(v)) / s.sigma
    v.mu <- mapply(function(r) mean(v.c^r), seq(3, order))
    v.gamma <- numeric(length(v.mu))
    if(order >= 3) v.gamma[1] <- v.mu[1]
    if(order >= 4) v.gamma[2] <- v.mu[2] - 3
    if(order >= 5) v.gamma[3] <- v.mu[3] - 10 * v.mu[1]
    if(order >= 6) v.gamma[4] <- v.mu[4] - 15 * v.mu[2] -
                                  10 * v.mu[1]^2 + 30
    v.gamma
}

# cf.quantile: Cornish-Fisher quantile function.
#
# Given a probability p, returns the quantile of a distribution with
# mean mu, standard deviation sigma, and standardised cumulants gamma.
# The expansion is:
#   Q(p) = mu + sigma * w(z)
# where z = Phi^{-1}(p) is the standard normal quantile and w(z) is
# the CF polynomial correction.  Each term in the polynomial
# corresponds to one order in the Edgeworth expansion:
#
#   Order 3: skewness correction   (g1/6)(z^2 - 1)
#   Order 4: kurtosis + skewness^2 cross-term
#   Order 5: 5th cumulant + mixed terms
#   Order 6: 6th cumulant + all cross-terms from Edgeworth inversion
#
# The expansion matches the first `order` cumulants by construction.
# By Theorem 2 (partial moment matching), this means the score bias
# is O(kappa_{order+1} / (order+1)!).
cf.quantile <- function(p, mu, sigma, gamma, order) {
    z <- qnorm(p)
    w <- z
    if(order >= 3) {
        g1 <- gamma[1]
        w <- w + (g1 / 6) * (z^2 - 1)
    }
    if(order >= 4) {
        g2 <- gamma[2]
        w <- w + (g2 / 24) * (z^3 - 3 * z) -
                 (g1^2 / 36) * (2 * z^3 - 5 * z)
    }
    if(order >= 5) {
        g3 <- gamma[3]
        w <- w + (g3 / 120) * (z^4 - 6 * z^2 + 3) -
                 (g1 * g2 / 24) * (z^4 - 5 * z^2 + 2) +
                 (g1^3 / 324) * (12 * z^4 - 53 * z^2 + 17)
    }
    if(order >= 6) {
        g4 <- gamma[4]
        # Leading term from the 6th cumulant plus all cross-terms
        # arising from the Edgeworth inversion [hall-bootstra-1992].
        w <- w + (g4 / 720) * (z^5 - 10 * z^3 + 15 * z) -
                 (g1 * g3 / 180) * (2 * z^5 - 17 * z^3 + 21 * z) -
                 (g2^2 / 192) * (z^5 - 7 * z^3 + 9 * z) +
                 (g1^2 * g2 / 216) * (4 * z^5 - 32 * z^3 + 41 * z) -
                 (g1^4 / 3888) * (14 * z^5 - 118 * z^3 + 155 * z)
    }
    mu + sigma * w
}

# cf.is.monotone: check that the CF quantile function is monotone
# over (0,1).  The CF expansion is a polynomial in z, so it can lose
# monotonicity when higher cumulants are large.  Non-monotonicity
# means the expansion is no longer a valid quantile function.
cf.is.monotone <- function(mu, sigma, gamma, order) {
    v.p <- seq(0.001, 0.999, length.out = 500)
    v.q <- cf.quantile(v.p, mu, sigma, gamma, order)
    all(diff(v.q) > 0)
}

# cf.safe.order: find the highest order <= `order` at which the CF
# expansion remains monotone.  If the full expansion is non-monotone,
# we reduce the order one step at a time until monotonicity holds.
# This is a safety mechanism: by Theorem 2, reducing the order from
# k to k-1 increases the bias bound by one factorial denominator
# (from kappa_{k+1}/(k+1)! to kappa_k/k!), but a non-monotone
# quantile function would be invalid.  In practice, order = 4
# (matching up to excess kurtosis) rarely loses monotonicity.
cf.safe.order <- function(mu, sigma, gamma, order) {
    if(sigma < .Machine[["double.eps"]]) return(2)
    while(order > 2) {
        b.mono <- tryCatch(cf.is.monotone(mu, sigma, gamma, order),
                           error = function(e) FALSE)
        if(isTRUE(b.mono)) break
        order <- order - 1
    }
    order
}

# +------------------------------------------------------------------------+
# | Empirical Bayes shrinkage (Theorem 5, sparse strata MSE)              |
# +------------------------------------------------------------------------+
#
# Theorem 5 (sparse strata MSE) shows that using raw per-stratum
# sample means to calibrate the generator yields MSE with a variance
# component proportional to sigma^2 / n_s for each stratum s.  When
# strata are small (n_s << n), this estimation noise dominates.
#
# Empirical Bayes shrinkage (normal-normal model) reduces this
# variance at the cost of a bias toward the grand mean:
#   tilde{theta}_s = w_s * hat{theta}_s + (1 - w_s) * bar{theta}
# where w_s = tau^2 / (tau^2 + sigma^2 / n_s).
#
# The optimal shrinkage rate (supplement S8) shows that this weight
# minimises the shrinkage contribution in the MSE bound.  When
# tau^2 >> sigma^2 * S / n (heterogeneous strata), w_s -> 1 (no
# shrinkage).  When tau^2 << sigma^2 * S / n (homogeneous), w_s -> 0
# (strong shrinkage).
#
# The shrinkage-calibration pipeline: the paper describes exactly
# this implementation: (i) compute per-stratum sample means and
# variances on the h-scale, (ii) estimate tau^2 by method of moments,
# (iii) compute w_s, (iv) replace stratum means with shrunk estimates.
shrink.moments <- function(l.strata.cf, l.strata.rows, v.cnt.nm) {
    # shrink per-stratum means toward the grand mean
    # using empirical Bayes (normal-normal model):
    #   theta_s ~ N(theta_grand, tau^2)
    #   theta_hat_s | theta_s ~ N(theta_s, sigma^2 / n_s)
    #   posterior mean = w_s * theta_hat_s + (1 - w_s) * theta_grand
    #   w_s = tau^2 / (tau^2 + sigma^2 / n_s)
    Map(function(nm) {
        # extract per-stratum mu and n_s
        v.mu <- mapply(function(l) l[[nm]][["mu"]], l.strata.cf)
        v.sigma <- mapply(function(l) l[[nm]][["sigma"]], l.strata.cf)
        v.ns <- mapply(length, l.strata.rows)
        b.valid <- !is.na(v.mu) & v.sigma > 0

        if(sum(b.valid) < 3) return(l.strata.cf)

        # grand mean and pooled within-stratum variance
        s.grand <- mean(v.mu[b.valid])
        s.sigma2.pool <- mean(v.sigma[b.valid]^2)

        # between-stratum variance (method of moments):
        #   Var(hat{theta}_s) = tau^2 + sigma^2 / n_s
        # so tau^2 = Var(hat{theta}_s) - sigma^2 * E[1/n_s].
        # Clamped at 0 (no negative variance).
        s.tau2 <- max(0, var(v.mu[b.valid]) -
                       s.sigma2.pool * mean(1 / v.ns[b.valid]))

        # shrinkage weights per stratum:
        #   w_s = tau^2 / (tau^2 + sigma^2 / n_s)
        # Large strata (large n_s) get w_s near 1 (little shrinkage).
        # Small strata (small n_s) get w_s near 0 (heavy shrinkage).
        v.w <- mapply(function(ns) {
            `if`(s.tau2 + s.sigma2.pool / ns > 0,
                 s.tau2 / (s.tau2 + s.sigma2.pool / ns),
                 1)
        }, v.ns)

        # apply shrinkage to means
        Map(function(s.key, w) {
            s.mu.old <- l.strata.cf[[s.key]][[nm]][["mu"]]
            if(!is.na(s.mu.old))
                l.strata.cf[[s.key]][[nm]][["mu"]] <<-
                    w * s.mu.old + (1 - w) * s.grand
        }, names(l.strata.cf), v.w)
    }, v.cnt.nm)
    l.strata.cf
}

# +------------------------------------------------------------------------+
# | KNN conditioning for continuous predictors                             |
# +------------------------------------------------------------------------+
#
# When there are no categorical variables to define strata, we create
# approximate strata via k-means clustering on the standardised
# continuous predictors.  Each cluster acts as a stratum S for the
# per-stratum moment computation (Theorem 2) and the EB shrinkage
# (Theorem 5).  The number of clusters is chosen as n / k, giving
# approximately k observations per stratum.
#
# This relates to the high-dimensional strata discussion: with continuous
# predictors, cross-classification is infeasible, so we discretise the
# predictor space into neighbourhoods that approximate the conditional
# distributions E[h(Y) | X in neighbourhood].
knn.strata <- function(d.cnt, k = 30) {
    # assign each observation to a KNN-based stratum
    # returns a factor of stratum labels (one per observation)
    s.n <- nrow(d.cnt)
    s.k <- min(k, s.n - 1)

    # standardise for distance computation
    m.sc <- scale(d.cnt)
    m.sc[is.na(m.sc)] <- 0

    # for each row, find its k nearest neighbours
    # return the sorted index set as a stratum key
    # (observations with the same neighbourhood share a stratum)
    # for efficiency, use the observation itself as the stratum centre
    # and assign strata via k-means clustering
    s.n.clusters <- max(1, round(s.n / s.k))
    l.km <- kmeans(m.sc, centers = s.n.clusters,
                   nstart = 5, iter.max = 50)
    factor(l.km[["cluster"]])
}

# +------------------------------------------------------------------------+
# | Main generator: gen.syn.family                                         |
# +------------------------------------------------------------------------+
#
# This is the core copula generator.  The algorithm:
#
# 1. Classify variables as categorical or continuous.  Detect supports.
# 2. Build support transforms h (Definition 3).
# 3. Stratify by categorical cross-classification (or KNN for
#    continuous-only data).  Strata define the conditioning sets
#    for the per-stratum moments E[h(Y) | S = s] and
#    Cov(h(Y) | S = s) required by Theorem 6 (direction-agnostic
#    fidelity).
# 4. For each stratum and each continuous variable, compute
#    conditional cumulants on the h-scale.  These parameterise
#    the CF quantile expansion (Theorem 2 mechanism).
# 5. Apply EB shrinkage to stratum means (Theorem 5).
# 6. Estimate the pooled within-stratum correlation on the h-scale.
#    By Sklar's theorem [sklar-fonction-1959], the joint distribution
#    of h(Y) can be decomposed into marginals (handled by CF) and a
#    copula (handled by correlated Gaussian draws).  We use the
#    Gaussian copula, which is the semiparametric estimator of
#    Genest et al. (1995) [genest-semipara-1995] applied to the
#    h-transformed variables.  At zero tail dependence
#    (lambda_U = lambda_L = 0) the tail-driven part of the
#    non-elliptical error bound vanishes (supplement S9), which is
#    the reason a Gaussian copula is used here.  That bound does not
#    control departures from ellipticity in the body, which a
#    tail-independent copula can carry, and the same section records
#    a vine [bedford-vinesne-2002] as the practical route to
#    non-zero tail dependence, improving distributional realism on
#    real trial data.
# 7. Generate correlated uniforms (pseudo-random or Sobol QMC),
#    apply the per-stratum CF quantile, then h^{-1} to get synthetic
#    observations on the original scale.
#
# Parameters:
#   d        : data.frame of original data
#   order    : CF expansion order (default 4 = match up to kurtosis).
#              By Proposition 1 (summary statistics), k = 2 suffices
#              for the full GLM summary.  Higher orders reduce the
#              O(kappa_{k+1}/(k+1)!) residual.
#   supports : named character vector of support types per variable
#              (auto-detected if NULL)
#   shrink   : whether to apply EB shrinkage (Theorem 5)
#   seed     : random seed for reproducibility
#   method   : "pseudo" (Mersenne Twister) or "sobol" (QMC sequence)
gen.syn.family <- function(d, order = 4, supports = NULL,
                           shrink = TRUE, seed = 6489,
                           method = c("pseudo", "sobol"),
                           strata.exact = FALSE) {
    method <- match.arg(method)
    set.seed(seed)
    s.n <- nrow(d)
    v.types <- mapply(id.cat, d)
    v.cat.nm <- names(v.types[v.types != "cnt"])
    v.cnt.nm <- names(v.types[v.types == "cnt"])

    # +--------------------------------------------------------------------+
    # | Stage 1: detect supports and build transforms (Definition 3)       |
    # +--------------------------------------------------------------------+

    if(is.null(supports)) {
        supports <- `names<-`(mapply(detect.support, d[v.cnt.nm]), v.cnt.nm)
    }
    v.missing <- v.cnt.nm[!v.cnt.nm %in% names(supports)]
    if(length(v.missing) > 0)
        supports <- c(supports,
                      `names<-`(mapply(detect.support, d[v.missing]), v.missing))

    # Build h and h^{-1} for each continuous variable.
    # All downstream operations (cumulants, CF quantile, correlation)
    # operate on h(Y), not on Y itself.
    l.transforms <- Map(function(nm) make.transform(supports[[nm]], d[[nm]]),
                        v.cnt.nm)

    # +--------------------------------------------------------------------+
    # | Stage 1b: zero-inflated variable decomposition                     |
    # +--------------------------------------------------------------------+
    # For zero-inflated variables, we decompose Y into:
    #   (a) a binary indicator I(Y > 0) -- joined to the categorical
    #       stratification, and
    #   (b) the positive part Y | Y > 0 -- synthesised via CF on the
    #       log scale.
    # This two-part model respects the support: log is defined only
    # on (0, inf), so the zero mass must be handled separately.  The
    # binary indicator enters the stratification, ensuring that the
    # conditional moments in Theorem 2 are computed over the correct
    # (positive) sub-population.
    v.zi.nm <- names(supports[supports == "zeroinfl"])
    d.aug <- d

    if(length(v.zi.nm) > 0) {
        Map(function(nm) {
            s.ind.nm <- paste0(nm, ".pos")
            d.aug[[s.ind.nm]] <<- factor(as.integer(d[[nm]] > 0))
        }, v.zi.nm)
        v.cat.nm <- c(v.cat.nm, paste0(v.zi.nm, ".pos"))
    }

    # +--------------------------------------------------------------------+
    # | Stage 2: stratification                                            |
    # +--------------------------------------------------------------------+
    # Strata define the conditioning variable S in the theoretical
    # results.  All conditional moment conditions (Theorem 2) and
    # direction-agnostic conditions (Theorem 6) are per-stratum.
    #
    # With categorical predictors, strata are the full cross-
    # classification of all factor levels.  Without categorical
    # predictors, we fall back to KNN-based clustering.
    #
    # The categorical bootstrap (sample with replacement from strata
    # labels) preserves P(S in .) -- condition (i) of Theorem 6.
    if(length(v.cat.nm) > 0) {
        v.strata <- interaction(d.aug[v.cat.nm], drop = TRUE)
    } else {
        v.strata <- knn.strata(d[v.cnt.nm])
    }
    l.strata.rows <- split(seq_len(s.n), v.strata)

    # Bootstrap the stratum labels: each synthetic row is assigned
    # to a stratum drawn with replacement from the empirical
    # stratum distribution.  This matches the stratum prevalences
    # pi_s = n_s / n used in Theorem 5.
    #
    # strata.exact = TRUE (design-Gram projection) draws the labels
    # Without replacement, a permutation of the original labels, so every
    # cross-classification cell keeps its exact original count and the
    # dummy-block Gram D'D is pinned (condition (i) of the design-Gram
    # proposition).  Only defined at synthetic size n.
    if(strata.exact) {
        if(s.n != nrow(d))
            stop("strata.exact = TRUE requires synthetic size == nrow(d)")
        v.strata.sampled <- sample(v.strata, s.n, replace = FALSE)
    } else {
        v.strata.sampled <- sample(v.strata, s.n, replace = TRUE)
    }

    # Reconstruct categorical columns from sampled strata.
    # Each synthetic observation inherits the categorical values of
    # its assigned stratum.
    if(length(v.cat.nm) > 0) {
        l.strata.cat <- Map(function(rows)
                            d.aug[rows[1], v.cat.nm, drop = FALSE],
                            l.strata.rows)
        d.syn.cat <- `rownames<-`(Reduce(rbind,
                                          Map(function(k)
                                              l.strata.cat[[as.character(k)]],
                                              v.strata.sampled)), NULL)
    } else {
        d.syn.cat <- NULL
    }

    # +--------------------------------------------------------------------+
    # | Stage 3: per-stratum CF parameters on h-scale (Theorem 2)          |
    # +--------------------------------------------------------------------+
    # For each stratum s and each continuous variable j, compute:
    #   mu_s    = E[h(Y_j) | S = s]           (stratum mean)
    #   sigma_s = sd(h(Y_j) | S = s)          (stratum sd)
    #   gamma_s = standardised cumulants 3..k  (shape parameters)
    #
    # These parameters define the per-stratum CF quantile function
    # Q_s(p) = mu_s + sigma_s * w(Phi^{-1}(p); gamma_s).
    #
    # By Theorem 2, matching these k cumulants yields score bias
    # bounded by O(kappa_{k+1} / (k+1)!).  cf.safe.order may
    # reduce the effective order if the full expansion is non-monotone.
    l.strata.cf <- Map(function(s.key, v.rows) {
        Map(function(nm) {
            v.raw <- d[[nm]][v.rows]
            l.tr <- l.transforms[[nm]]

            # For zero-inflated variables, condition on Y > 0.
            # The zero mass is handled by the binary indicator in
            # the categorical stratification.
            if(supports[[nm]] == "zeroinfl") {
                v.raw <- v.raw[v.raw > 0]
                if(length(v.raw) == 0)
                    return(list(mu = NA, sigma = 0,
                                gamma = rep(0, max(order - 2, 0)),
                                order = 2, skip = TRUE))
            }

            # Transform to h-scale and compute cumulants.
            v.h <- l.tr[["h"]](v.raw)
            v.h <- v.h[is.finite(v.h)]
            if(length(v.h) < 2)
                return(list(mu = `if`(length(v.h) == 1, v.h, 0),
                            sigma = 0,
                            gamma = rep(0, max(order - 2, 0)),
                            order = 2, skip = FALSE))

            s.mu <- mean(v.h)
            s.sigma <- sd(v.h)
            v.gamma <- std.cumulants(v.h, order)
            s.order <- cf.safe.order(s.mu, s.sigma, v.gamma, order)
            list(mu = s.mu, sigma = s.sigma,
                 gamma = v.gamma, order = s.order, skip = FALSE)
        }, v.cnt.nm)
    }, names(l.strata.rows), l.strata.rows)

    # +--------------------------------------------------------------------+
    # | Stage 4: empirical Bayes shrinkage (Theorem 5)                    |
    # +--------------------------------------------------------------------+
    # Shrink per-stratum means toward the grand mean to reduce the
    # MSE of the synthetic estimator in sparse strata.  See the
    # shrink.moments function above for the full derivation.
    if(shrink)
        l.strata.cf <- shrink.moments(l.strata.cf, l.strata.rows, v.cnt.nm)

    # +--------------------------------------------------------------------+
    # | Stage 5: pooled within-stratum correlation (Theorem 6)            |
    # +--------------------------------------------------------------------+
    # Theorem 6 (direction-agnostic fidelity) requires matching
    # Cov(W | S = s) where W = (h_1(Y_1), ..., h_p(Y_p)) is the
    # vector of support-transformed variables.  The Gaussian copula
    # achieves this by modelling the within-stratum dependence via a
    # single correlation matrix.
    #
    # By Sklar's theorem [sklar-fonction-1959], the joint distribution
    # of W | S = s is determined by:
    #   (a) the marginals F_j(w_j | S = s), and
    #   (b) the copula C(u_1, ..., u_p | S = s).
    # The CF quantile handles (a).  The Gaussian copula handles (b).
    #
    # We pool the within-stratum correlations across strata by
    # computing residuals after removing stratum means, then taking
    # the pairwise correlation of the residuals.  This is the
    # semiparametric approach of Genest et al. (1995)
    # [genest-semipara-1995] adapted to the stratified setting.

    # Transform all continuous variables to h-scale.
    m.h <- mapply(function(nm) {
        l.tr <- l.transforms[[nm]]
        v.raw <- d[[nm]]
        if(supports[[nm]] == "zeroinfl")
            v.raw <- ifelse(v.raw > 0, v.raw, NA)
        v.h <- l.tr[["h"]](v.raw)
        v.h[!is.finite(v.h)] <- NA
        v.h
    }, v.cnt.nm)

    # Compute residuals by subtracting stratum means.
    # The residual correlation is the within-stratum correlation,
    # which enters the Gaussian copula.
    m.resid <- mapply(function(nm, j) {
        v <- m.h[, j]
        v.fitted <- mapply(function(rows) mean(v[rows], na.rm = TRUE),
                           l.strata.rows)
        v - v.fitted[as.integer(v.strata)]
    }, v.cnt.nm, seq_along(v.cnt.nm))
    s.p.cnt <- length(v.cnt.nm)
    m.cor.within <- `if`(s.p.cnt > 1,
                         cor(m.resid, use = "pairwise.complete.obs"),
                         matrix(1, 1, 1))
    m.cor.within[is.na(m.cor.within)] <- 0
    diag(m.cor.within) <- 1

    # Ensure positive definiteness.  The pairwise correlation matrix
    # from incomplete data can have negative eigenvalues.  NearPD
    # projects to the nearest PD correlation matrix.
    v.eigen <- eigen(m.cor.within, only.values = TRUE)[["values"]]
    if(any(v.eigen < 0))
        m.cor.within <- as.matrix(Matrix::nearPD(m.cor.within,
                                                  corr = TRUE)[["mat"]])

    # +--------------------------------------------------------------------+
    # | Stage 6: generate correlated uniforms via Gaussian copula          |
    # +--------------------------------------------------------------------+
    # The Gaussian copula maps N(0, Sigma) draws to Uniform(0,1)
    # via the probability integral transform: U_j = Phi(Z_j).
    # The correlation structure of U mirrors the within-stratum
    # correlation of the h-transformed variables.
    #
    # Two methods are available:
    #   "pseudo" -- standard Gaussian draws (Mersenne Twister)
    #   "sobol"  -- Sobol quasi-random sequence with digital shift
    #              randomisation, giving better space-filling for
    #              small n (lower discrepancy).
    #
    # In both cases, the raw draws are "whitened" (decorrelated by
    # the inverse Cholesky of their sample covariance) then
    # re-correlated using the Cholesky of the target within-stratum
    # correlation matrix.  This ensures the copula correlation
    # matches m.cor.within exactly, not just approximately.
    #
    # Copula joint calibration: the paper notes that
    # the Gaussian copula is naturally suited to direction-agnostic
    # fidelity because it models the joint symmetrically.  At
    # lambda_U = lambda_L = 0 the tail-driven part of the
    # non-elliptical error bound vanishes [nelsen-introduc-2007,
    # Section 5.4], while the body departure that bound does not
    # control remains the operative condition (supplement S9).
    if(method == "sobol") {
        m.u <- qrng::sobol(s.n, s.p.cnt, randomize = "digital.shift",
                           seed = seed)
        if(s.p.cnt > 1) {
            m.z <- qnorm(m.u)
            m.z <- m.z %*% solve(chol(cov(m.z)))
            m.z <- m.z %*% chol(m.cor.within)
            m.u <- pnorm(m.z)
        }
    } else {
        m.z <- matrix(rnorm(s.n * s.p.cnt), nrow = s.n, ncol = s.p.cnt)
        if(s.p.cnt > 1) {
            m.z <- m.z %*% solve(chol(cov(m.z)))
            m.z <- m.z %*% chol(m.cor.within)
        }
        m.u <- pnorm(m.z)
    }

    # +--------------------------------------------------------------------+
    # | Stage 7: CF quantile -> h^{-1} (synthesis)                         |
    # +--------------------------------------------------------------------+
    # For each synthetic observation i assigned to stratum s:
    #   1. Look up the per-stratum CF parameters (mu_s, sigma_s, gamma_s).
    #   2. Apply cf.quantile(u_ij, mu_s, sigma_s, gamma_s, order_s)
    #      to get a value on the h-scale.
    #   3. Apply h^{-1} to map back to the original support.
    #
    # This is where the theoretical pipeline completes:
    #   Correlated Uniform (copula)
    #     -> Phi^{-1} (quantile of standard normal)
    #     -> CF polynomial correction (match k cumulants)
    #     -> mu + sigma * w (location-scale on h-scale)
    #     -> h^{-1} (back to original support)
    #
    # The result satisfies the conditional moment conditions of
    # Theorem 2 within each stratum, and the joint correlation
    # structure of Theorem 6 across variables.
    l.syn.idx <- split(seq_len(s.n), v.strata.sampled)
    m.syn.cnt <- matrix(NA_real_, nrow = s.n, ncol = s.p.cnt)
    m.syn.cnt <- `colnames<-`(m.syn.cnt, v.cnt.nm)

    Map(function(s.key, v.syn.idx) {
        l.cf <- l.strata.cf[[s.key]]
        Map(function(nm, j) {
            l.par <- l.cf[[nm]]
            l.tr  <- l.transforms[[nm]]

            # For zero-inflated variables, set the synthetic value
            # to 0 for observations whose binary indicator is 0.
            # Only the positive-part observations get CF-synthesised.
            if(supports[[nm]] == "zeroinfl") {
                s.ind.nm <- paste0(nm, ".pos")
                b.pos <- d.syn.cat[[s.ind.nm]][v.syn.idx] == "1"
                m.syn.cnt[v.syn.idx[!b.pos], j] <<- 0
                if(!any(b.pos)) return(NULL)
                v.syn.pos <- v.syn.idx[b.pos]
            } else {
                v.syn.pos <- v.syn.idx
            }

            # Degenerate stratum: constant value (sigma = 0) or
            # all-zero positive part (skip = TRUE).
            if(isTRUE(l.par[["skip"]]) || l.par[["sigma"]] == 0) {
                m.syn.cnt[v.syn.pos, j] <<-
                    l.tr[["hinv"]](`if`(is.na(l.par[["mu"]]), 0, l.par[["mu"]]))
                return(NULL)
            }

            # Apply CF quantile using the correlated uniform draws,
            # then invert the support transform.
            v.h <- cf.quantile(m.u[v.syn.pos, j],
                               l.par[["mu"]], l.par[["sigma"]],
                               l.par[["gamma"]], l.par[["order"]])
            m.syn.cnt[v.syn.pos, j] <<- l.tr[["hinv"]](v.h)
        }, v.cnt.nm, seq_along(v.cnt.nm))
    }, names(l.syn.idx), l.syn.idx)

    # +--------------------------------------------------------------------+
    # | Stage 8: assemble output data.frame                                |
    # +--------------------------------------------------------------------+
    # Combine the bootstrapped categorical columns with the
    # CF-synthesised continuous columns, preserving the original
    # column order.  Zero-inflated indicator columns (*.pos) are
    # dropped -- they were internal to the stratification.
    v.orig.cat <- names(v.types[v.types != "cnt"])
    d.syn <- `if`(!is.null(d.syn.cat),
                  cbind(d.syn.cat[v.orig.cat],
                        as.data.frame(m.syn.cnt)),
                  as.data.frame(m.syn.cnt))
    `rownames<-`(d.syn[colnames(d)], NULL)
}
