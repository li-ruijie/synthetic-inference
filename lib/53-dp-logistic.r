#!/usr/bin/env r
# +------------------------------------------------------------+
# |  DP-Projection for binary logistic regression              |
# +------------------------------------------------------------+
#
# The virtual-y trick requires y'_i = y_i + mu(beta_DP) - mu(beta_hat),
# which is continuous.  For binary y this is not acceptable, since the
# analyst's glm(binomial) requires y in {0, 1}.
#
# Solution: Probabilistic bernoulli sampling.
#
#   1. Compute beta_hat from logistic regression.
#   2. Output-perturb: beta_DP = beta_hat + N(0, sigma^2 I).
#   3. Generate synthetic X via copula (or use original X for public-X).
#   4. For each synthetic obs: p_i = expit(x_i' beta_DP),
#      sample y_syn_i ~ Bernoulli(p_i).
#   5. Analyst fits glm(binomial), gets beta_syn ~ beta_DP + O(1/sqrt(n)).
#
# Exact recovery is not possible (binary y constraint), but for large
# n the sampling error is small relative to the DP noise.  Decision
# agreement should match the OLS case.
#
# +------------------------------------------------------------+
.sobol.main <- TRUE
source("lib/00-shared.r")
source("lib/gen-syn-family.r")
source("lib/07-calibrated-generator.r")
source("lib/09-symmetric-generator.r")
source("lib/11-projection-generator.r")
source("lib/28-dp-projection.r")
# Sub-project B (sensitive-X): reuse A's DP copula helpers from 35-dp-full.r.
.dp.full.no.main <- TRUE
source("lib/35-dp-full.r")
rm(.dp.full.no.main)
rm(.sobol.main)

# +------------------------------------------------------------+
# |  DP-Logistic generator                                     |
# +------------------------------------------------------------+
# Arguments:
#   d       original data frame with binary y
#   y.nm    response name (values in {0, 1})
#   x.nms   predictor names
#   eps, delta  DP parameters
#   seed    RNG seed
#   x.public   if TRUE, use original X (no X-noise).  Else DP-privatise X from
#              x.spec (independent marginals by default, use.copula=TRUE adds A's
#              Mixed-margin Gaussian copula for joint X-fidelity.
#              The categoricals join the copula through their K-1
#              indicator scores exactly as in gen.syn.dp.full, assumptions A1-A6
#              there, score.method = "adaptive" PIT)
#   score.method  copula score transform: "adaptive" (default, PIT through each
#              covariate's released DP histogram, handles skewed marginals),
#              "normal" (public normal-score PIT, recovers correlation at eps>=5
#              on location-scale margins), or "box" (legacy, noise-dominated).
#              See the DP guide
#   back.map   categorical back-map for the mixed copula, "probit" (default,
#              recommended) or "logit" (attenuates, see A2 on gen.syn.dp.full)
#
# Returns: list with d.syn, beta.hat, beta.dp, sigma, Delta.
gen.syn.dp.logistic <- function(d, y.nm, x.nms = NULL,
                                  eps = 10, delta = 1e-6,
                                  x.public = TRUE, x.spec = NULL,
                                  sens.method = c("subagg", "local"),
                                  clip.lo = NULL, clip.hi = NULL, m = NULL,
                                  mu.split = c(marg = 0.05, corr = 0.05, coef = 0.90),
                                  eps.split = c(X = 0.2, beta = 0.8),
                                  eps.corr.frac = 0.5, use.copula = FALSE,
                                  score.method = "adaptive",
                                  back.map = c("probit", "logit"), seed = 42) {
    sens.method <- match.arg(sens.method)
    back.map <- match.arg(back.map)
    if(!x.public && is.null(x.spec))
        stop("x.public = FALSE requires x.spec (per-column type/breaks/levels)")
    if(is.null(x.nms))
        x.nms <- if(!is.null(x.spec)) names(x.spec) else setdiff(names(d), y.nm)
    if(!all(d[[y.nm]] %in% c(0, 1, TRUE, FALSE)))
        stop("y must be binary (0/1 or TRUE/FALSE)")
    d[[y.nm]] <- as.integer(d[[y.nm]])
    # Sensitive-X: code categorical predictors as factors before the fit so
    # beta_hat's dummy names match the sampled X_syn model matrix (else the
    # analyst-side refit silently drops that coefficient).
    if(!x.public)
        for(nm in x.nms) {
            spec <- x.spec[[nm]]
            if(spec[["type"]] == "categorical")
                d[[nm]] <- factor(as.character(d[[nm]]), levels = spec[["levels"]])
        }
    fml <- reformulate(x.nms, response = y.nm)
    score.method <- match.arg(score.method, c("box", "normal", "adaptive"))
    n <- nrow(d)

    # copula eligibility, hoisted because the budget split needs it.
    v.cont.nms <- if(x.public) character(0) else
        x.nms[vapply(x.nms, function(nm)
                     x.spec[[nm]][["type"]] == "continuous", logical(1))]
    do.copula <- !x.public && isTRUE(use.copula) && length(v.cont.nms) >= 2L

    # ---- budget allocation --------------------------------------------------
    # subagg: every channel is a Gaussian mechanism, so the release composes in
    # One mu-GDP quadrature.  With public X there is no X channel and the whole
    # mu goes to the coefficient.  The Bernoulli draw from expit(x* beta_DP) is
    # post-processing and costs nothing (supplement S6, third extension).
    if(sens.method == "subagg") {
        if(is.null(clip.lo) || is.null(clip.hi))
            stop("sens.method = \"subagg\" requires the public per-coefficient ",
                 "bounds clip.lo and clip.hi")
        if(abs(sum(mu.split) - 1) > 1e-8) stop("mu.split must sum to 1")
        s.mu.total <- gdp.mu.from.eps.delta(eps, delta)
        s.mu.coef <- if(x.public) s.mu.total
                     else s.mu.total * sqrt(mu.split[["coef"]])
        s.mu.corr <- if(do.copula) s.mu.total * sqrt(mu.split[["corr"]]) else 0
        s.mu.marg <- if(x.public) 0
                     else if(do.copula) s.mu.total * sqrt(mu.split[["marg"]])
                     else gdp.compose(s.mu.total *
                              sqrt(c(mu.split[["marg"]], mu.split[["corr"]])))
        s.mu.each  <- if(x.public) 0 else s.mu.marg / sqrt(length(x.nms))
        s.eps.each <- NULL; s.eps.corr <- 0; s.eps.X <- NA_real_
        s.eps.beta <- NA_real_
    } else {
        warning("sens.method = \"local\" calibrates beta to the DFBETA LOCAL ",
                "sensitivity, which is a function of the realised data. The ",
                "release therefore carries NO formal (eps, delta)-DP guarantee ",
                "Retained only for reproducibility. ",
                "runs. Use sens.method = \"subagg\".")
        if(abs(sum(eps.split) - 1) > 1e-8) stop("eps.split must sum to 1")
        # Public X spends the whole budget on beta.  Sensitive X splits it.
        if(x.public) { s.eps.X <- 0; s.eps.beta <- eps }
        else { s.eps.X <- eps * eps.split[["X"]]
               s.eps.beta <- eps * eps.split[["beta"]] }
        s.eps.corr <- if(do.copula) s.eps.X * eps.corr.frac else 0
        s.eps.each <- (s.eps.X - s.eps.corr) / length(x.nms)
        s.mu.total <- NA_real_; s.mu.coef <- NA_real_
        s.mu.corr  <- NA_real_; s.mu.marg <- NA_real_; s.mu.each <- NULL
    }

    dp.marg.cont <- function(x, breaks)
        if(sens.method == "subagg") dp.hist.continuous(x, breaks, mu = s.mu.each)
        else dp.hist.continuous(x, breaks, eps = s.eps.each)
    dp.marg.cat <- function(x, levels)
        if(sens.method == "subagg") dp.cat(x, levels, mu = s.mu.each)
        else dp.cat(x, levels, eps = s.eps.each)

    # Step 1: fit logistic regression (non-private reference, evaluation only)
    mo <- glm(fml, d, family = binomial(link = "logit"))
    v.beta.hat <- coef(mo)
    m.X <- model.matrix(fml, d)

    # Step 2: DP beta
    if(sens.method == "subagg") {
        # Subsample-and-aggregate, binomial/logit.  Worst-case per-coordinate
        # sensitivity s_j = (hi_j - lo_j)/m over a data-independent partition,
        # so the release is mu_coef-GDP and is sound under full-row replacement,
        # which is the relation that matters once X is sensitive.  The binomial
        # has no dispersion, so the aggregate stays p-dimensional.
        l.sa <- subsample.aggregate.beta.loglink(
            d, y.nm, x.nms, "binomial", "logit",
            eps = NA_real_, delta = delta, m = m,
            clip.lo = clip.lo, clip.hi = clip.hi,
            seed = seed + 10L, mu = s.mu.coef)
        v.beta.dp <- l.sa[["beta.dp"]]
        s.sigma <- NA_real_; s.LS.std <- NA_real_; v.s <- NULL
    } else {
        l.sa <- NULL
        # Deprecated DFBETA local sensitivity.  Data-dependent, so no formal
        # (eps, delta)-DP.  Kept only for reproducibility.
        v.sd <- apply(m.X, 2, stats::sd)
        v.s <- ifelse(v.sd > 1e-8, v.sd, 1)
        m.dfb <- dfbeta(mo)
        s.LS.std <- 2 * max(sqrt(rowSums(sweep(m.dfb, 2, v.s, "*")^2)))
        s.sigma <- s.LS.std * dp.gm.multiplier(s.eps.beta, delta)
        set.seed(seed + 1L)
        v.eta.std <- rnorm(length(v.beta.hat), 0, s.sigma)
        s.norm <- sqrt(sum(v.eta.std^2))
        s.max <- 3 * sqrt(length(v.beta.hat)) * s.sigma
        if(s.norm > s.max) v.eta.std <- v.eta.std * (s.max / s.norm)
        v.beta.dp <- `names<-`((v.beta.hat * v.s + v.eta.std) / v.s,
                               names(v.beta.hat))
    }

    # Step 3: construct X_syn
    if(x.public) {
        d.syn   <- d
        m.X.syn <- m.X
        m.R.dp  <- NULL
    } else {
        set.seed(seed)
        l.dp.X <- list()
        if(do.copula) {
            # Mixed-margin copula.  The categoricals join
            # the copula through their K-1 indicator scores, exactly as V4's
            # Stage 1 in 35-dp-full.r (whose helpers this file sources).  The
            # released vech(S) covers p_total columns at the same M^2 p / n
            # sensitivity form, and sampling runs through the calibrated
            # back-map with the T25b de-attenuation.  Assumptions A1-A6 on
            # gen.syn.dp.full apply verbatim.
            v.cat.nms <- setdiff(x.nms, v.cont.nms)
            for(nm in v.cont.nms)
                l.dp.X[[nm]] <- dp.marg.cont(d[[nm]], x.spec[[nm]][["breaks"]])
            for(nm in v.cat.nms)
                l.dp.X[[nm]] <- dp.marg.cat(d[[nm]], x.spec[[nm]][["levels"]])
            s.tau <- if(score.method == "box") 1 / (2 * n) else 0.005
            s.M   <- qnorm(1 - s.tau)
            m.W.cont <- vapply(v.cont.nms, function(nm) {
                v.bx <- range(x.spec[[nm]][["breaks"]])
                dp.copula.scores(d[[nm]], v.bx[1], v.bx[2], s.tau,
                                 method = score.method, dp.hist = l.dp.X[[nm]])
            }, numeric(n))
            m.W <- cbind(m.W.cont,
                         do.call(cbind, lapply(v.cat.nms, function(nm)
                             cat.indicator.scores(d[[nm]], l.dp.X[[nm]], s.M))))
            m.R.dp <- if(sens.method == "subagg")
                          dp.copula.corr(m.W, NULL, delta, s.M, mu = s.mu.corr)
                      else dp.copula.corr(m.W, s.eps.corr, delta, s.M)
            d.syn <- sample.from.dp.mixed.copula(m.R.dp, l.dp.X, v.cont.nms,
                                                 v.cat.nms, back.map, n)
        } else {
            m.R.dp <- NULL
            d.syn  <- data.frame(.row = seq_len(n))
            for(nm in x.nms) {
                spec <- x.spec[[nm]]
                if(spec[["type"]] == "continuous") {
                    l.dp.X[[nm]] <- dp.marg.cont(d[[nm]], spec[["breaks"]])
                    d.syn[[nm]]  <- sample.from.dp.hist(l.dp.X[[nm]], n)
                } else {
                    l.dp.X[[nm]] <- dp.marg.cat(d[[nm]], spec[["levels"]])
                    d.syn[[nm]]  <- sample.from.dp.cat(l.dp.X[[nm]], n)
                }
            }
            d.syn[[".row"]] <- NULL
        }
        d.syn   <- d.syn[x.nms]
        m.X.syn <- model.matrix(reformulate(x.nms), d.syn)
    }

    # Step 4: sample y_syn ~ Bernoulli(p_i)
    v.beta.aligned <- rep(0, ncol(m.X.syn))
    names(v.beta.aligned) <- colnames(m.X.syn)
    v.common <- intersect(names(v.beta.dp), colnames(m.X.syn))
    v.beta.aligned[v.common] <- v.beta.dp[v.common]
    v.eta <- as.numeric(m.X.syn %*% v.beta.aligned)
    v.p <- plogis(v.eta)
    set.seed(seed + 2L)
    v.y.syn <- rbinom(nrow(m.X.syn), 1, v.p)
    d.syn[[y.nm]] <- v.y.syn

    # Not private, evaluation only: beta.hat is the full-data fit, and subagg
    # carries its own non-private intermediates.  The DP release is d.syn.
    list(d.syn       = d.syn,
         beta.hat    = v.beta.hat,
         beta.dp     = v.beta.dp,
         sens.method = sens.method,
         eps         = eps,
         delta       = delta,
         # mu-GDP accounting (subagg).  mu.check re-composes in quadrature.
         mu.total    = s.mu.total,
         mu.channels = c(marg = s.mu.marg, corr = s.mu.corr, coef = s.mu.coef),
         mu.check    = if(sens.method == "subagg")
                           gdp.compose(c(s.mu.marg, s.mu.corr, s.mu.coef))
                       else NA_real_,
         subagg      = l.sa,
         # legacy (sens.method = "local") diagnostics
         sigma.beta  = s.sigma,
         Delta       = s.LS.std,
         col.scale   = v.s,
         x.public    = x.public,
         use.copula  = do.copula,
         dp.corr     = m.R.dp,
         eps.allocation = if(sens.method == "local")
                              c(X = s.eps.X, beta = s.eps.beta, corr = s.eps.corr)
                          else NULL,
         # Under subagg the whole release is calibrated to one (eps, delta) via
         # mu, so delta no longer accumulates across channels.
         delta.total = if(sens.method == "subagg") delta
                       else if(do.copula) 2 * delta else delta)
}

# +------------------------------------------------------------+
# |  Monte Carlo driver                                         |
# +------------------------------------------------------------+
# This file defines a mechanism.  The Monte Carlo demonstration lives in
# 36-mc-all-variants.r, the single DP-variant driver, which takes every bound
# from dp-public-bounds.r and refuses to run on an invented one.
