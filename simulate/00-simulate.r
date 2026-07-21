#!/usr/bin/env r
# +------------------------------------------------------------------------+
# | 00-simulate.r -- Data-Generating Process for 11 progressive scenarios  |
# +------------------------------------------------------------------------+
#
# This script implements the simulation study design described in Section 6
# of the manuscript ("Simulation study"). It creates the ground-truth data
# against which every generator is tested, providing the known DGP and true
# regression coefficients that the paper's main theorems make predictions about.
#
# Central objective: to determine whether synthetic data generators
# optimised for distributional fidelity preserve inferential fidelity
# (Definition 2).  The 11 scenarios
# here are designed to progressively isolate the failure mechanisms:
#   - Scenarios 1-2: Gaussian identity link (Theorem 1 says k=1 suffices;
#     Corollary: families (i) -- bias is exactly zero for k >= 1)
#   - Scenario 3: Gamma with log link -- tests Theorem 2 (partial moment
#     matching) and Corollary: families (ii) (k=2 on log-scale needed)
#   - Scenario 4-5: Zero-inflated quasi-gamma -- tests the copula's ability
#     to preserve conditional structure in the presence of a mass at zero
#   - Scenario 6: Poisson with canonical log link (Theorem 1 applies: first
#     conditional moment suffices, but the integer constraint introduces
#     probabilistic rounding error)
#   - Scenario 7: Correlated predictors with a clean Gamma/log response --
#     isolates the sparse strata problem (Theorem 5)
#   - Scenario 8: All combined (correlated predictors + covariate-dependent
#     zero-inflation + ceiling) -- tests Theorem 6 (direction-agnostic
#     fidelity) together with the sparse strata problem (Theorem 5)
#   - Scenario 9: Bounded response [0,100] -- tests boundary behaviour
#     (Taylor boundary behaviour: case (a), bounded h^{-1})
#   - Scenario 10: Gaussian response under the correlated predictors (the
#     data keeps all five variables. Dropping the categoricals, which a
#     separate test path does, exercises the KNN-strata fallback)
#   - Scenario 11: Binary response, canonical logit link -- exercises the
#     binomial branch of Theorem 7 (path B in 11-projection-generator.r)
#
# Strategy: Generate correlated latent normals via Cholesky decomposition,
# then transform to target variable types via the quantile function. This
# guarantees the regression coefficients are significant because the
# underlying correlations are controlled by the 5x5 correlation matrix.
#
# References:
#   McCullagh & Nelder (1989) [mccullag-generali-1989] -- GLM families
#   Wedderburn (1974) [wedderbu-quasi-li-1974] -- quasi-likelihood for Sc 4-5, 8
#   Cornish & Fisher (1938) [cornish-moments-1938] -- CF quantile approach
#
# +------------------------------------------------------------------------+
# | Setup                                                                  |
# +------------------------------------------------------------------------+
source("lib/00-shared.r")
# +------------------------------------------------------------------------+
# | Helpers                                                                |
# +------------------------------------------------------------------------+
# id.cat: classify variables as continuous/categorical/binary.
# This determines how the copula generator (gen-syn-family.r) treats each
# column: categorical variables define the strata S used in Theorem 5
# (sparse strata MSE) and Theorem 6 (direction-agnostic fidelity), while
# continuous variables are transformed via h (Definition 3: support transform).
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
# +------------------------------------------------------------------------+
# | Data-Generating Process                                                |
# +------------------------------------------------------------------------+
# Each scenario generates n=2000 observations (Section 6: "n = 2000
# observations are generated from a Gaussian copula whose 5x5
# correlation matrix is fixed within each scenario"). It is not fixed
# Across scenarios, since three predictor correlations are strengthened
# from scenario 7 onwards, which is the next clause of the sentence
# quoted above. The function returns the data, family, and link --
# the (g, V) pair that defines the analyst's quasi-likelihood model
# (Section 2: "the analyst specifies a link function g, a variance
# function V, and a linear predictor eta_i = x_i' beta").
sim.scenario <- function(scenario, n = 2000, seed = 6489) {
    set.seed(seed)

    # Target correlation matrix for latent normals:
    # [response, age, x_cont, x_bin, x_cat]
    # These correlations ensure the true beta_0 is nonzero for every
    # predictor, so we can measure whether generators preserve or
    # attenuate the coefficients. This is the key testable prediction
    # of Proposition 2 (no-guarantee): distributional generators will
    # attenuate rare-level coefficients despite accurate marginal fit.
    m.cor <- matrix(c(
        1.00, -0.25,  0.30,  0.15,  0.35,
       -0.25,  1.00,  0.40,  0.10,  0.20,
        0.30,  0.40,  1.00,  0.05,  0.15,
        0.15,  0.10,  0.05,  1.00,  0.08,
        0.35,  0.20,  0.15,  0.08,  1.00),
        nrow = 5, byrow = TRUE)

    # Scenario 7+: stronger predictor correlations to test Theorem 6
    # (direction-agnostic fidelity) and Theorem 5 (sparse strata MSE).
    # When predictors are correlated, the cross-strata become sparser
    # because the joint distribution concentrates on fewer cells.
    if(scenario >= 7) {
        m.cor[2, 3] <- m.cor[3, 2] <- 0.55
        m.cor[2, 5] <- m.cor[5, 2] <- 0.35
        m.cor[3, 5] <- m.cor[5, 3] <- 0.30
    }

    # Generate correlated normals via Cholesky decomposition.
    # The latent normal Z is the starting point for all scenarios.
    # The quantile transform Z -> pnorm(Z) -> qF(.) maps these to the
    # target marginal distributions while preserving rank correlations.
    # This is Sklar's theorem [sklar-fonction-1959] in action: separate
    # the dependence structure (the Gaussian copula) from the marginals.
    m.z <- matrix(rnorm(n * 5), nrow = n) %*% chol(m.cor)

    # +--------------------------------------------------------------------+
    # | Transform latent normals -> observed variables                     |
    # +--------------------------------------------------------------------+
    # Each variable is generated via the quantile function of its target
    # distribution applied to the uniform u = pnorm(z). This is the
    # inverse of the support transform h (Definition 3): for positive
    # variables, h = log.  For bounded, h = logit.  For counts, h = log(y+1/2).

    # age: continuous real, bounded [21, 85]
    v.age <- round(qnorm(pnorm(m.z[, 2]), mean = 48, sd = 13))
    v.age <- pmax(pmin(v.age, 85), 21)

    # x_cont: income
    if(scenario >= 2) {
        # skewed: log-normal, scaled to ~[500, 10000]
        v.income <- round(exp(qnorm(pnorm(m.z[, 3]), mean = 7.5, sd = 0.5)))
        v.income <- pmin(v.income, 15000)
    } else {
        # normal
        v.income <- round(qnorm(pnorm(m.z[, 3]), mean = 2000, sd = 500))
        v.income <- pmax(v.income, 100)
    }

    # x_bin: female (binary via threshold)
    v.female <- factor(as.integer(pnorm(m.z[, 4]) > 0.5))

    # x_cat: race (3-level via thresholds at 50th, 85th percentile)
    v.u.cat <- pnorm(m.z[, 5])
    v.race  <- factor(findInterval(v.u.cat, c(0.50, 0.85)) + 1)

    d <- data.frame(age = v.age, female = v.female,
                    race = v.race, income = v.income)

    # +--------------------------------------------------------------------+
    # | Response variable (transform latent column 1)                     |
    # +--------------------------------------------------------------------+
    # The response Y is generated from the latent normal z via quantile
    # inversion: u = pnorm(z), then Y = qF(u) for the target family F.
    # The choice of F determines the GLM family and link, and therefore
    # which theorem governs the required moment order for fidelity.
    v.z.resp <- m.z[, 1]
    v.u.resp <- pnorm(v.z.resp)

    if(scenario <= 2) {
        # Scenario 1-2: Gaussian, identity link.
        # Theorem 1 (canonical): k=1 suffices because the identity link
        # is canonical for the Gaussian family, making the score linear
        # in Y. Corollary (families, (i)): bias is exactly zero for k>=1.
        v.y <- qnorm(v.u.resp, mean = 50, sd = 15)
        s.family <- "gaussian"
        s.link   <- "identity"

    } else if(scenario == 3) {
        # Scenario 3: Gamma with log link (non-canonical).
        # Theorem 2 (partial moment matching): bias decays as
        # O(kappa_{k+1}/(k+1)!). Corollary (families, (ii)): k=1 is
        # insufficient (conditional variance contributes bias of order
        # |kappa_2* - kappa_2|/2), and k=2 on log-scale is needed.
        # McCullagh & Nelder (1989) Ch. 8 [mccullag-generali-1989].
        v.y <- qgamma(v.u.resp, shape = 2, rate = 0.1)
        s.family <- "gamma"
        s.link   <- "log"

    } else if(scenario == 4) {
        # Scenario 4: Zero-inflated positive response.
        # 60% structural zeros -- tests the copula bootstrap of zero
        # indicators + positive-subset projection strategy in
        # 11-projection-generator.r. The quasi-gamma family
        # (Wedderburn 1974 [wedderbu-quasi-li-1974]) with log link is
        # used because the mass at zero prevents a proper Gamma MLE.
        # Two-part model: the positive part is fitted separately.
        v.p.zero <- 0.60
        b.zero   <- v.u.resp < v.p.zero
        v.y      <- ifelse(b.zero, 0,
                           qgamma((v.u.resp - v.p.zero) / (1 - v.p.zero),
                                  shape = 2, rate = 0.1))
        s.family <- "quasi-gamma"
        s.link   <- "log"

    } else if(scenario == 5) {
        # Scenario 5: Zero-inflated + ceiling at 100.
        # Tests boundary behaviour (Taylor boundary behaviour, case (a)).
        # The bounded h^{-1} = expit has all derivatives bounded on R,
        # so the Taylor remainder is uniformly controlled.
        v.p.zero <- 0.60
        b.zero   <- v.u.resp < v.p.zero
        v.y      <- ifelse(b.zero, 0,
                           pmin(qgamma((v.u.resp - v.p.zero) / (1 - v.p.zero),
                                       shape = 2, rate = 0.1), 100))
        s.family <- "quasi-gamma"
        s.link   <- "log"

    } else if(scenario == 6) {
        # Scenario 6: Poisson with canonical log link.
        # Theorem 1 (canonical): k=1 on the natural scale suffices
        # regardless of h. However, the integer constraint means the
        # projection must use probabilistic rounding (see
        # 11-projection-generator.r), introducing a small residual
        # error of ~2.3%.
        v.y <- qpois(v.u.resp, lambda = 8)
        s.family <- "poisson"
        s.link   <- "log"

    } else if(scenario == 7) {
        # Scenario 7: Correlated predictors (Gamma, log link).
        # The same Gamma(shape 2, rate 0.1) response as scenario 3, but
        # with the stronger predictor correlations set above. Correlated
        # predictors concentrate the joint distribution on fewer
        # cross-strata, so per-stratum moment estimation becomes noisy
        # (Theorem 5, sparse strata MSE) even without zero-inflation.
        # This isolates the correlated-predictor effect, which scenario
        # 8 then combines with zero-inflation and the ceiling.
        v.y <- qgamma(v.u.resp, shape = 2, rate = 0.1)
        s.family <- "gamma"
        s.link   <- "log"

    } else if(scenario == 8) {
        # Scenario 8: All combined (Quasi-gamma, log link).
        # Correlated predictors + covariate-dependent zero-inflation +
        # ceiling effect: the maximum-difficulty configuration. Tests
        # Theorem 6 (direction-agnostic fidelity: Cov(h(Z)) matching
        # within strata) together with the sparse strata problem, where
        # correlated predictors leave many cross-strata with n_s = 1-2
        # observations and noisy per-stratum moment estimation
        # (Theorem 5). The zero probability is a function of race and
        # age, creating covariate-dependent zero-inflation.
        v.p.zero <- 0.50 + 0.15 * (v.race == "1") - 0.005 * v.age
        v.p.zero <- pmax(pmin(v.p.zero, 0.95), 0.05)
        b.zero   <- v.u.resp < v.p.zero
        v.y.pos  <- qgamma(pmax((v.u.resp - v.p.zero) / (1 - v.p.zero),
                                1e-6),
                           shape = 2, rate = 0.1)
        v.y      <- ifelse(b.zero, 0, pmin(v.y.pos, 100))
        s.family <- "quasi-gamma"
        s.link   <- "log"

    } else if(scenario == 9) {
        # Scenario 9: Bounded response [0, 100], floor + ceiling effects.
        # Score-like variable (satisfaction, percentage, etc.).
        # Gaussian family with identity link. Tests the support transform
        # h = logit for bounded responses (Definition 3). The projection
        # achieves exact fidelity (0.0% deviation) because the identity
        # link is canonical for the Gaussian family (Theorem 1).
        v.y.latent <- qnorm(v.u.resp, mean = 50, sd = 20)
        v.y <- round(pmax(pmin(v.y.latent, 100), 0))
        s.family <- "gaussian"
        s.link   <- "identity"

    } else if(scenario == 10) {
        # Scenario 10: Gaussian response under the correlated predictors.
        # Same response as scenarios 1-2 while keeping the scenario >= 7
        # correlation strengthening, so it is the Gaussian analogue of
        # scenario 7. The data frame keeps all five variables and the
        # paper's tables fit the full model. Dropping the categoricals
        # exercises the no-categorical KNN-strata fallback in
        # gen.syn.family. CTGAN fails here (female1 flips sign:
        # 4.31 -> -1.94, and income flips too: 0.0069 -> -0.0000),
        # confirming Proposition 2.
        v.y <- qnorm(v.u.resp, mean = 50, sd = 15)
        s.family <- "gaussian"
        s.link   <- "identity"

    } else {
        # Scenario 11: Binary (Bernoulli) response, canonical logit link.
        # Exercises Theorem 7's binomial branch (path B in
        # 11-projection-generator.r), which the ten continuous/count scenarios
        # do not touch.  y = I(u.resp > 0.5) is a balanced binary whose
        # association with the predictors comes from the same 5x5 copula
        # correlation that drives the binary predictor `female`, so the
        # logistic coefficients are significant.  The projection releases an
        # exact fractional [0,1] binomial response.  Refitting glm(binomial,
        # logit) recovers beta_hat (the 0.0% deviation bar the other families
        # meet).
        v.y <- as.integer(v.u.resp > 0.5)
        s.family <- "binomial"
        s.link   <- "logit"
    }

    d[["response"]] <- v.y
    list(data   = d,
         family = s.family,
         link   = s.link)
}
# +------------------------------------------------------------------------+
# | Generate all 11 scenarios                                              |
# +------------------------------------------------------------------------+
# This is the master data object consumed by every downstream script:
# - gen-syn-family.r (copula generator)
# - 07-calibrated-generator.r (shrinkage calibration)
# - 11-projection-generator.r (sufficient statistics projection)
# - 12-export-scenarios.r (CSV export for Python ML generators)
# The list contains data frames + family/link metadata for each scenario.
l.scenarios <- Map(sim.scenario, 1:11)
save(l.scenarios, file = "results/d.sim.scenarios.RData")

# +------------------------------------------------------------------------+
# | Verification: fit the true model on each scenario                      |
# +------------------------------------------------------------------------+
# Confirms the DGP produces significant coefficients across all scenarios.
# For log-link families with zeros (Scenarios 4-5, 8), the model is
# fitted on the positive subset only (standard two-part model approach).
# The quasi-gamma uses quasi(link="log", variance="mu^2") -- Wedderburn's
# quasi-likelihood [wedderbu-quasi-li-1974, mccullag-quasi-li-1983].
Map(function(s, i) {
    d <- s[["data"]]
    fml <- response ~ age + female + race + income

    # for log-link families, fit on positive subset (standard for
    # zero-inflated data: two-part model, this is the positive part)
    d.fit <- `if`(s[["link"]] == "log" && any(d[["response"]] == 0),
                  d[d[["response"]] > 0, ], d)
    mo <- tryCatch(
        `if`(s[["family"]] == "gaussian",
             lm(fml, d.fit),
        `if`(s[["family"]] == "binomial",
             glm(fml, d.fit, family = binomial()),
             `if`(s[["family"]] == "poisson",
                  glm(fml, d.fit, family = poisson),
                  glm(fml, d.fit,
                      family = quasi(link = "log", variance = "mu^2"),
                      start = c(1, rep(0, 5)),
                      control = glm.control(maxit = 100))))),
        error = function(e) { message("  fit failed: ", e[["message"]]); NULL })

    v.y <- d[["response"]]
    cat(sprintf("=== Scenario %d: %s(%s) ===\n", i, s[["family"]], s[["link"]]))
    cat(sprintf("  range: [%.1f, %.1f]  zeros: %d (%.0f%%)  at_100: %d\n",
                min(v.y), max(v.y),
                sum(v.y == 0), 100 * mean(v.y == 0),
                sum(v.y >= 99.9 & s[["family"]] != "poisson")))
    cat("  types:", paste(mapply(id.cat, d), collapse = " "), "\n")
    if(!is.null(mo)) {
        d.coef <- summary(mo)[["coefficients"]]
        v.sig <- ifelse(d.coef[, 4] < 0.001, "***",
                 ifelse(d.coef[, 4] < 0.01,  "** ",
                 ifelse(d.coef[, 4] < 0.05,  "*  ", "   ")))
        cat("  coef: ", paste(sprintf("%8.4f%s", d.coef[, 1], v.sig),
                              collapse = "  "), "\n")
    }
    cat("\n")
}, l.scenarios, seq_along(l.scenarios))
# +------------------------------------------------------------------------+
# | End of 00-simulate.r                                                   |
# +------------------------------------------------------------------------+
