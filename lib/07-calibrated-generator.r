#!/usr/bin/env r
# +------------------------------------------------------------------------+
# | 07-calibrated-generator.r -- Shrinkage-aware calibrated generator      |
# +------------------------------------------------------------------------+
#
# Implements the calibrated synthetic data generator, the middle ground
# between the raw copula generator (06-gen-syn-family, distributional
# fidelity only) and the projection generator (11-projection, exact
# coefficient reproduction). Where the copula generator sacrifices
# inferential fidelity and the projection generator sacrifices privacy,
# the calibrated generator trades a controlled amount of coefficient
# accuracy for genuine novelty in individual observations.
#
# Architecture -- 5-step shrinkage-calibration pipeline
# (the shrinkage-calibration pipeline):
#   (i)   Compute per-stratum conditional means E[h(Y)|X_cat=s] and
#         variances Var[h(Y)|X_cat=s] on the h-scale for each
#         continuous variable.
#   (ii)  Estimate the between-stratum variance tau^2 by method of
#         moments: tau^2_hat = Var(theta_hat_s) - sigma^2 * mean(1/n_s),
#         floored at zero (non-negative estimator).
#   (iii) Compute per-stratum shrinkage weights:
#             w_s = tau^2 / (tau^2 + sigma^2 / n_s)
#         This is the key formula from Theorem 5 (sparse strata MSE).
#         When n_s is large, w_s -> 1 (trust the stratum).
#         When n_s is small, w_s -> 0 (shrink to grand mean).
#   (iv)  Form shrunk targets:
#             theta^EB_s = w_s * theta_hat_s + (1 - w_s) * theta_grand
#         This is the empirical Bayes posterior mean under the
#         normal-normal hierarchical model.
#   (v)   Generate copula output via gen.syn.family, then calibrate:
#         shift and rescale h(Y*)|s to match the shrunk targets
#         theta^EB_s (mean) and sigma^EB_s (sd), then invert via h^{-1}.
#
# Theoretical connections:
#
# James-Stein [james-estimati-1961]: the positive-part JS estimator
#   dominates the MLE for S >= 3 strata simultaneously. The shrinkage
#   weights w_s implement this: when sum(b.ok) >= 3, we shrink.  Below
#   that threshold, we fall back to the raw per-stratum estimates.
#
# Efron & Morris (1973) [efron-steinse-1973]: the method-of-moments
#   estimator for tau^2 that we use here (floored at zero) comes from
#   their parametric EB framework. Their relative savings loss (RSL)
#   result says EB shrinkage captures a fraction 1 - 2/S of the
#   maximum possible improvement over the MLE when the shrinkage
#   centre is fixed, and 1 - 3/S when that centre is itself
#   estimated, their Section 7, "the effective dimension dropping by
#   one". This file shrinks toward theta_grand, an estimated centre,
#   so 1 - 3/S is the constant that applies here.
#
# Efron & Morris (1977) [efron-steinsp-1977]: the toxoplasmosis
#   example warns that JS shrinkage across truly unrelated problems
#   can mislead. Here, the strata are subgroups of the same DGP,
#   so the hierarchical (exchangeability) assumption is natural.
#
# Cornish-Fisher [cornish-moments-1938]: the CF quantile expansion
#   (in gen-syn-family.r) generates novel values from matched moments
#   up to order 4-6. This file's calibration step is a post-hoc
#   correction on top of the CF output, so the two are complementary:
#   CF provides the shape, calibration corrects the location/scale.
#
# The optimal shrinkage rate (supplement S8): the crossover where raw
#   per-stratum estimation becomes preferable to shrinkage occurs at
#   n ~ S * sigma^2 / tau^2. Below this threshold, shrinkage reduces
#   MSE.  Above it, the bias term dominates. The calibrated generator
#   adapts automatically via the w_s weights.
#
# Theorem 2 and the mean squared error decomposition: matching the
#   first k conditional moments on the h-scale bounds the score bias
#   and so guarantees GLM coefficient consistency. The copula
#   provides the joint dependence structure.  The calibration ensures
#   the per-stratum moments are accurate (shrunk, not noisy).
#
# Role in the paper:
# This is the "calibrated" method in the six-method comparison
# (Section 6). It represents the best a distributional method can do
# when it also borrows strength across strata: better than raw copula
# (which ignores sparse-stratum noise), but still not structurally
# exact like projection (Theorem 7).
#
# +------------------------------------------------------------------------+
# | Setup                                                                  |
# +------------------------------------------------------------------------+
source("lib/00-shared.r")
source("lib/gen-syn-family.r")
# +------------------------------------------------------------------------+
# | Post-hoc calibration on the h-scale                                   |
# +------------------------------------------------------------------------+
# calibrate() adjusts the continuous columns of a synthetic dataset so that
# per-stratum conditional moments of h(Y) match target moments. It operates
# entirely on the h-scale (the link-variance function from the GLM family),
# then inverts via h^{-1}, so the support constraint (positive, count,
# bounded, zero-inflated) is respected by construction.
#
# When shrink=TRUE (the default and the paper's recommended setting):
#   - Targets are the EB-shrunk moments (step iv of the pipeline), not
#     the raw per-stratum sample moments.
#   - This means the calibration and the shrinkage cooperate: the
#     calibration forces the synthetic moments to equal the *shrunk*
#     targets, rather than undoing the shrinkage that gen.syn.family
#     already applied.
#   - The cooperation is essential: without it, calibration would push
#     sparse-stratum moments back toward noisy raw estimates, defeating
#     the purpose of shrinkage.
#
# When shrink=FALSE:
#   - Targets are the raw per-stratum sample moments.
#   - This is the "calibrated but unshrunk" variant used as a
#     comparison in the simulation study.
#
# order controls how many moments to match:
#   order=1: match conditional mean only (location shift).
#   order=2: match conditional mean and variance (location + scale).
calibrate <- function(d.syn, d.orig, order = 2, supports = NULL,
                      shrink = TRUE) {
    # Classify each column as categorical or continuous.
    # id.cat (from gen-syn-family.r) returns "cat", "bin", or "cnt".
    v.types <- mapply(id.cat, d.orig)
    v.cat.nm <- names(v.types[v.types != "cnt"])
    v.cnt.nm <- names(v.types[v.types == "cnt"])

    # If there are no categorical or no continuous variables, calibration
    # is vacuous (nothing to stratify by, or nothing to adjust).
    if(length(v.cat.nm) == 0 || length(v.cnt.nm) == 0) return(d.syn)

    # Detect supports for continuous variables if not supplied.
    # detect.support (from gen-syn-family.r) infers "real", "positive",
    # "zeroinfl", etc. from the data range.
    if(is.null(supports))
        supports <- `names<-`(mapply(detect.support, d.orig[v.cnt.nm]), v.cnt.nm)

    # Build the h / h^{-1} transform pair for each continuous variable.
    # These are the same transforms used by gen.syn.family:
    #   real     -> identity
    #   positive -> log / exp
    #   zeroinfl -> log / exp (applied to positive part only)
    #   count    -> log(x+0.5) / round(exp(x)-0.5)
    #   bounded  -> qlogis / plogis (scaled to [lo, hi])
    l.tr <- Map(function(nm) make.transform(supports[[nm]], d.orig[[nm]]),
                v.cnt.nm)

    # Form strata by crossing all categorical variables.
    # v.strata.orig and v.strata.syn are factor vectors of stratum labels.
    # The strata define the conditioning sets X_cat = s in the paper.
    v.strata.orig <- interaction(d.orig[v.cat.nm], drop = TRUE)
    v.strata.syn  <- interaction(d.syn[v.cat.nm], drop = TRUE)
    l.orig.rows <- split(seq_len(nrow(d.orig)), v.strata.orig)
    l.syn.rows  <- split(seq_len(nrow(d.syn)), v.strata.syn)

    # +--------------------------------------------------------------------+
    # | Pipeline steps (i)-(iv): compute EB-shrunk calibration targets   |
    # +--------------------------------------------------------------------+
    # For each continuous variable, compute per-stratum raw moments on
    # the h-scale, then shrink them toward the grand mean/pooled variance.
    #
    # Step (i): per-stratum means and SDs on the h-scale.
    # Step (ii)-(iii): estimate tau^2 and compute w_s (inside the
    #   shrinkage blocks below).
    # Step (iv): form shrunk targets theta^EB_s.
    #
    # Two shrinkage operations are performed independently:
    #   A) Mean shrinkage: normal-normal EB, same formula as
    #      shrink.moments() in gen-syn-family.r.
    #   B) Variance shrinkage: shrink per-stratum variance toward the
    #      pooled variance, with weight (n_s-1)/(n_s+1). This is a
    #      simpler df-based weight (not full EB for variances), but
    #      sufficient because variance estimation is less noisy than
    #      mean estimation for moderate n_s.
    l.targets <- Map(function(nm) {
        h    <- l.tr[[nm]][["h"]]
        b.zi <- supports[[nm]] == "zeroinfl"

        # Initialise per-stratum storage vectors.
        v.mu.raw <- v.sd.raw <- v.ns <- `names<-`(
            rep(NA_real_, length(l.orig.rows)), names(l.orig.rows))

        # Step (i): compute raw per-stratum means and SDs on the h-scale.
        # For zero-inflated variables, condition on the positive part
        # (the zero indicator is handled separately as a categorical).
        Map(function(s.key) {
            v.raw <- d.orig[[nm]][l.orig.rows[[s.key]]]
            if(b.zi) v.raw <- v.raw[v.raw > 0]
            v.h <- h(v.raw)
            v.h <- v.h[is.finite(v.h)]
            v.ns[[s.key]]     <<- length(v.h)
            if(length(v.h) >= 1) v.mu.raw[[s.key]] <<- mean(v.h)
            if(length(v.h) >= 2) v.sd.raw[[s.key]] <<- sd(v.h)
        }, names(l.orig.rows))

        # If shrinkage is disabled, return raw moments as targets.
        if(!shrink) return(list(mu = v.mu.raw, sd = v.sd.raw))
        # Steps (ii)-(iv) for means: EB shrinkage under the normal-normal
        # hierarchical model.
        #   Level 1 (within-stratum):  Y_is | theta_s ~ N(theta_s, sigma^2)
        #   Level 2 (between-stratum): theta_s ~ N(theta_grand, tau^2)
        #
        # The James-Stein condition S >= 3 is checked here: we only
        # shrink when sum(b.ok) >= 3. Below that threshold, the JS
        # estimator does not dominate the MLE [james-estimati-1961],
        # so we keep raw estimates.
        b.ok <- !is.na(v.mu.raw) & !is.na(v.sd.raw) & v.ns >= 2
        v.mu.sh <- v.mu.raw
        if(sum(b.ok) >= 3) {
            # Step (ii): estimate hyperparameters.
            # s.grand = theta_grand (grand mean across strata).
            # s.sigma2 = pooled within-stratum variance (average of s^2_s).
            # s.tau2 = between-stratum variance by method of moments,
            #   following Efron & Morris (1973) [efron-steinse-1973]:
            #     tau^2_hat = Var(theta_hat_s) - sigma^2 * mean(1/n_s)
            #   floored at zero (non-negative estimator). The subtracted
            #   term corrects for the sampling variability of theta_hat_s.
            s.grand  <- mean(v.mu.raw[b.ok])
            s.sigma2 <- mean(v.sd.raw[b.ok]^2)
            s.tau2   <- max(0, var(v.mu.raw[b.ok]) -
                         s.sigma2 * mean(1 / v.ns[b.ok]))

            # Step (iii): per-stratum shrinkage weights.
            # w_s = tau^2 / (tau^2 + sigma^2 / n_s)  -- Theorem 5.
            # Large n_s -> sigma^2/n_s small -> w_s near 1 (trust stratum).
            # Small n_s -> sigma^2/n_s large -> w_s near 0 (shrink to grand).
            # The crossover (optimal shrinkage rate, supplement S8) is at
            # n ~ S * sigma^2 / tau^2.
            v.w <- mapply(function(ns, ok) {
                if(!ok) return(1)
                s.d <- s.tau2 + s.sigma2 / ns
                `if`(s.d > 0, s.tau2 / s.d, 1)
            }, v.ns, b.ok)

            # Step (iv): form EB posterior means.
            # theta^EB_s = w_s * theta_hat_s + (1 - w_s) * theta_grand
            # For strata with insufficient data (!b.ok), keep the raw
            # estimate (effectively w_s = 1).
            v.mu.sh <- ifelse(b.ok,
                              v.w * v.mu.raw + (1 - v.w) * s.grand,
                              v.mu.raw)
        }
        # Variance shrinkage: shrink per-stratum variance toward the
        # pooled variance. The weight (n_s - 1)/(n_s + 1) is a
        # degrees-of-freedom-based shrinkage: with n_s - 1 df from the
        # data and 1 effective df from the prior (pooled estimate), the
        # weight reflects the relative information. This is simpler than
        # full EB for variances but adequate for the calibration step,
        # which primarily depends on the mean shrinkage for coefficient
        # accuracy (the variance shrinkage stabilises the scale).
        #
        # Same S >= 3 threshold as for means.
        b.ok.sd <- !is.na(v.sd.raw) & v.ns >= 2
        v.sd.sh <- v.sd.raw
        if(sum(b.ok.sd) >= 3) {
            s.var.pool <- mean(v.sd.raw[b.ok.sd]^2)
            v.w.var <- ifelse(b.ok.sd, (v.ns - 1) / (v.ns + 1), 0)
            # Shrunk variance = weighted average of per-stratum and pooled.
            # pmax(0, ...) ensures non-negative variance.
            v.sd.sh <- ifelse(b.ok.sd,
                              sqrt(pmax(0, v.w.var * v.sd.raw^2 +
                                           (1 - v.w.var) * s.var.pool)),
                              sqrt(s.var.pool))
            v.sd.sh[is.na(v.sd.raw)] <- NA_real_
        }
        list(mu = v.mu.sh, sd = v.sd.sh)
    }, v.cnt.nm)
    l.targets <- `names<-`(l.targets, v.cnt.nm)

    # +--------------------------------------------------------------------+
    # | Pipeline step (v): apply calibration to copula output             |
    # +--------------------------------------------------------------------+
    # For each continuous variable and each stratum, shift and rescale
    # the synthetic h(Y*) values to match the shrunk targets from above.
    #
    # The calibration is a two-moment affine transform on the h-scale:
    #   h(Y*)_calibrated = mu^EB_s + (h(Y*) - mean(h(Y*))) * (sd^EB_s / sd(h(Y*)))
    # This preserves the rank order and correlation structure produced by
    # the copula, while correcting the first two conditional moments.
    # After calibration, h^{-1} maps back to the original scale,
    # automatically respecting the support.
    Map(function(nm) {
        h    <- l.tr[[nm]][["h"]]
        hinv <- l.tr[[nm]][["hinv"]]
        b.zi <- supports[[nm]] == "zeroinfl"

        Map(function(s.key) {
            # Skip strata that do not appear in both original and
            # synthetic data (can happen with rare strata).
            if(!s.key %in% names(l.orig.rows)) return(NULL)
            if(!s.key %in% names(l.syn.rows))  return(NULL)
            v.oi <- l.orig.rows[[s.key]]
            v.si <- l.syn.rows[[s.key]]
            # Need at least 2 observations for mean+variance.
            if(length(v.si) < 2 || length(v.oi) < 2) return(NULL)

            v.orig.raw <- d.orig[[nm]][v.oi]
            v.syn.raw  <- d.syn[[nm]][v.si]

            # For zero-inflated variables, only calibrate the positive
            # part. The zero indicator was already handled as a
            # categorical variable (binary bootstrap in gen.syn.family).
            if(b.zi) {
                b.pos.o <- v.orig.raw > 0
                b.pos.s <- v.syn.raw > 0
                if(sum(b.pos.o) < 2 || sum(b.pos.s) < 2) return(NULL)
                v.orig.h <- h(v.orig.raw[b.pos.o])
                v.syn.h  <- h(v.syn.raw[b.pos.s])
            } else {
                v.orig.h <- h(v.orig.raw)
                v.syn.h  <- h(v.syn.raw)
            }

            # Remove non-finite values (e.g. log(0) = -Inf).
            v.orig.h <- v.orig.h[is.finite(v.orig.h)]
            v.syn.h  <- v.syn.h[is.finite(v.syn.h)]
            if(length(v.orig.h) < 2 || length(v.syn.h) < 2) return(NULL)

            # Order-1 calibration: match conditional mean on the h-scale.
            # Uses the shrunk target mu^EB_s from step (iv).
            # Additive shift: h(Y*) <- h(Y*) + (mu^EB_s - mean(h(Y*))).
            # This is the location correction that directly affects the
            # GLM intercept and slope consistency (Theorem 2).
            s.mu.o <- l.targets[[nm]][["mu"]][[s.key]]
            if(is.na(s.mu.o)) s.mu.o <- mean(v.orig.h)
            s.mu.s <- mean(v.syn.h)
            v.syn.h <- v.syn.h + (s.mu.o - s.mu.s)

            if(order >= 2) {
                # Order-2 calibration: match conditional variance.
                # Uses the shrunk target sd^EB_s from the variance
                # shrinkage above. Multiplicative rescaling centred on
                # the (already corrected) mean:
                #   h(Y*) <- mu^EB_s + (h(Y*) - mu^EB_s) * (sd^EB_s / sd_syn)
                # This preserves the mean while correcting the spread.
                # The guard s.sd.s > eps avoids division by zero for
                # degenerate (constant) synthetic strata.
                s.sd.o <- l.targets[[nm]][["sd"]][[s.key]]
                if(is.na(s.sd.o)) s.sd.o <- sd(v.orig.h)
                s.sd.s <- sd(v.syn.h)
                if(s.sd.s > .Machine[["double.eps"]])
                    v.syn.h <- s.mu.o + (v.syn.h - mean(v.syn.h)) *
                               (s.sd.o / s.sd.s)
            }

            # Invert: h-scale -> original scale via h^{-1}.
            # Because the affine transform on the h-scale preserves
            # the range of h, and h^{-1} maps that range back to the
            # support, the calibrated values automatically satisfy the
            # support constraint (positive, count, bounded).
            if(b.zi) {
                v.syn.raw[b.pos.s] <- hinv(v.syn.h)
                d.syn[[nm]][v.si] <<- v.syn.raw
            } else {
                d.syn[[nm]][v.si] <<- hinv(v.syn.h)
            }
        }, names(l.orig.rows))
    }, v.cnt.nm)

    d.syn
}
# +------------------------------------------------------------------------+
# | Combined generator: copula + calibration                               |
# +------------------------------------------------------------------------+
# gen.syn.calibrated is the user-facing function. It runs the full
# 5-step pipeline:
#   Steps (i)-(iv) are inside calibrate() (shrunk target computation).
#   Step (v) is the actual calibration adjustment.
#   The copula + CF quantile generation is step 0 (gen.syn.family).
#
# Parameters:
#   order.cf  -- Cornish-Fisher expansion order (default 4).
#                Controls how many cumulants the CF quantile matches
#                when generating novel values [cornish-moments-1938].
#                Higher orders capture skewness (3) and kurtosis (4).
#   order.cal -- Calibration order (default 2 = mean + variance).
#                order=1 matches only the conditional mean.  Order=2
#                also matches the conditional variance.
#   shrink    -- Whether to use EB-shrunk targets (default TRUE).
#                TRUE = the recommended setting from the paper.
#   seed      -- Master seed for reproducibility (6489 per Section 6).
gen.syn.calibrated <- function(d, order.cf = 4, order.cal = 2,
                               supports = NULL, shrink = TRUE,
                               seed = 6489) {
    # Step 0: generate genuinely new data via copula + CF quantile.
    # gen.syn.family (from gen-syn-family.r) produces novel observations
    # (not copies of originals) with the correct marginal shape and
    # within-stratum correlation structure. However, the per-stratum
    # conditional moments may be noisy due to the CF parametric
    # approximation and finite-sample estimation.
    d.syn <- gen.syn.family(d, order = order.cf, supports = supports,
                            shrink = shrink, seed = seed)

    # Steps (i)-(v): calibrate on the h-scale to match EB-shrunk
    # conditional moments. This corrects the parametric approximation
    # error from step 0, ensuring that E*[h(Y)^k | X_cat=s] matches
    # the target (shrunk) moments. By Theorem 2 and the mean squared
    # error decomposition, this guarantees GLM coefficient consistency.
    calibrate(d.syn, d, order = order.cal, supports = supports,
              shrink = shrink)
}
# +------------------------------------------------------------------------+
# | Test harness                                                           |
# +------------------------------------------------------------------------+
# Compares three methods on a subset of the 10 DGP scenarios:
#   1. Original data (ground truth coefficients).
#   2. Uncalibrated gen.syn.family (copula + CF, no post-hoc correction).
#   3. Calibrated gen.syn.calibrated (copula + CF + EB calibration).
#
# The test prints coefficient values, significance stars, and percentage
# deviation from the original. The calibrated method should show smaller
# deviations, especially for sparse strata where shrinkage helps most.
#
# Guarded by if(!exists(".sobol.main")) so this block does not execute
# when the file is source()'d by the Sobol sensitivity analysis pipeline.
if(!exists(".sobol.main")) {
load("results/d.sim.scenarios.RData")

# Per-scenario support specifications. These must match the DGP:
# Scenarios 1-10 use different response distributions (Gaussian,
# Gamma, Poisson, zero-inflated Gamma, etc.), so the support for
# "response" varies. "age" is always real-valued.  "income" is
# positive in most scenarios.
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

# fit.model: fit the GLM specified by the scenario.
# Uses the same predictor set as the DGP (age, female, race, income).
# For non-Gaussian, non-Poisson families (e.g. Gamma with log link),
# uses quasi(link="log", variance="mu^2") which is the Gamma variance
# function, with explicit starting values to aid convergence.
# For zero-inflated responses with log link, excludes zeros before
# fitting (the zeros are modelled separately via the binary indicator).
fit.model <- function(d, family, link) {
    v.pred <- intersect(c("age", "female", "race", "income"), colnames(d))
    fml <- reformulate(v.pred, response = "response")
    tryCatch(
        `if`(family == "gaussian",
             lm(fml, d),
             `if`(family == "poisson",
                  glm(fml, d, family = poisson),
                  {
                      d.fit <- `if`(link == "log" && any(d[["response"]] == 0),
                                    d[d[["response"]] > 0, ], d)
                      s.np <- ncol(model.matrix(fml, d.fit))
                      glm(fml, d.fit,
                          family = quasi(link = "log", variance = "mu^2"),
                          start = c(1, rep(0, s.np - 1)),
                          control = glm.control(maxit = 100))
                  })),
        error = function(e) NULL)
}

# sig: extract significance stars from a fitted model.
# Used for the side-by-side comparison table.
sig <- function(mo) {
    v.p <- summary(mo)[["coefficients"]][, 4]
    ifelse(v.p < 0.001, "***",
    ifelse(v.p < 0.01,  "** ",
    ifelse(v.p < 0.05,  "*  ", "   ")))
}

# Test scenarios: a subset of the 10 DGPs.
# Scenarios 1,3,4,6,7,9 cover Gaussian, Gamma (positive), zero-inflated,
# Poisson, zero-inflated Gamma, and bounded responses -- a representative
# cross-section of the families in Section 6.
v.test <- c(1, 3, 4, 6, 7, 9)

# For each test scenario: fit the GLM on original, uncalibrated synthetic,
# and calibrated synthetic data, then compare coefficients.
# The percentage deviation (%+.1f%%) quantifies the inferential fidelity
# loss. The calibrated method should show smaller deviations because:
#   - The EB shrinkage (Theorem 5) reduces noise in sparse-stratum moments.
#   - The calibration (step v) forces the synthetic moments to match the
#     shrunk targets exactly (up to the specified order).
#   - Together, these produce coefficient estimates closer to the original.
Map(function(i) {
    s <- l.scenarios[[i]]
    d.orig <- s[["data"]]
    sup <- l.supports[[i]]
    mo.orig <- fit.model(d.orig, s[["family"]], s[["link"]])
    if(is.null(mo.orig)) return(NULL)

    # Generate uncalibrated synthetic data (copula + CF only, no post-hoc).
    d.uncal <- gen.syn.family(d.orig, order = 4, supports = sup)
    mo.uncal <- fit.model(d.uncal, s[["family"]], s[["link"]])

    # Generate calibrated synthetic data (copula + CF + EB calibration).
    d.cal <- gen.syn.calibrated(d.orig, order.cf = 4, order.cal = 2,
                                supports = sup)
    mo.cal <- fit.model(d.cal, s[["family"]], s[["link"]])

    # Extract and compare coefficients.
    v.orig <- coef(mo.orig)
    v.uncal <- `if`(!is.null(mo.uncal), coef(mo.uncal), rep(NA, length(v.orig)))
    v.cal   <- `if`(!is.null(mo.cal), coef(mo.cal), rep(NA, length(v.orig)))

    # Percentage deviation from original coefficients.
    v.pct.u <- round(100 * (v.uncal - v.orig) / abs(v.orig), 1)
    v.pct.c <- round(100 * (v.cal - v.orig) / abs(v.orig), 1)

    v.sig.o <- sig(mo.orig)
    v.sig.u <- `if`(!is.null(mo.uncal), sig(mo.uncal), rep("   ", length(v.orig)))
    v.sig.c <- `if`(!is.null(mo.cal), sig(mo.cal), rep("   ", length(v.orig)))

    cat(sprintf("=== Scenario %d: %s(%s) ===\n", i, s[["family"]], s[["link"]]))
    d.out <- data.frame(
        orig = sprintf("%7.4f%s", v.orig, v.sig.o),
        uncal = sprintf("%7.4f%s (%+.1f%%)", v.uncal, v.sig.u, v.pct.u),
        calib = sprintf("%7.4f%s (%+.1f%%)", v.cal, v.sig.c, v.pct.c))
    print(d.out)
    cat("\n")
}, v.test)
} # end if(!exists(".sobol.main"))
# +------------------------------------------------------------------------+
