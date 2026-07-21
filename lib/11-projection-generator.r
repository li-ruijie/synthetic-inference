#!/usr/bin/env r
# +------------------------------------------------------------------------+
# | 11-projection-generator.r -- Sufficient statistics projection          |
# +------------------------------------------------------------------------+
#
# Main result.  This file implements the constructive proof of Theorem 7
# (sufficient statistics projection).  The idea: instead of training a
# generative model that approximately matches moments of Y, we *directly
# project* copula-generated Y_cop onto the (p+1)-dimensional affine
# constraint surface where every GLM inference quantity is exact:
#
#   (C1)  X' w (y* - mu_hat)  =  0          [score equation, p constraints]
#   (C2)  l_sat(y*)           =  l_sat(y)    [saturated log-likelihood, 1]
#
# Theorem 7 shows that, for X of full column rank, any y* satisfying
# (C1)-(C2) gives:
#   - exact beta_hat*  = beta_hat
#   - exact dispersion phi* = phi
#   - exact standard errors, t-statistics, and p-values
# for every submodel S of the saturated model.  Only p+1 constraints
# are imposed on n unknowns, so the constraint surface has dimension
# n - p - 1 -- the same as the model's residual degrees of freedom.
#
# The paper's Table 3 confirms: projection achieves 0.0% deviation
# (sd = 0.0) across 100 Monte Carlo replications for all continuous GLMs.
# This is the numerical verification of the constructive proof.
#
# Strategy by family:
#
#   Gaussian / identity (OLS):
#     Corollary (closed-form projection).
#     y* = X beta_hat + c (I - H) y_cop,  c = ||( I-H)y|| / ||(I-H)y_cop||
#     The fitted component comes from the original data (encodes X'y),
#     the residual component carries n-p-1 df of genuine novelty from the
#     copula, and the scalar c enforces RSS* = RSS (constraint C2).
#
#   Positive-support GLM with log link (Gamma, quasi-gamma, Poisson):
#     The log-scale Newton step.  Work on z = log(y), solve the score
#     equation X' w (exp(z*) - mu_hat) = 0 by Newton iteration.
#     Because y = exp(z) > 0 by construction, the support constraint
#     is automatically satisfied.  Minimum-norm updates lie in col(X),
#     so the copula residual component is untouched.
#
#   Quasi-gamma (non-canonical, positive):
#     Alternating y-scale projection: project score to zero, then clip
#     negative values, repeat.  Newton on log-scale would collapse
#     Pearson dispersion for quasi families.
#
#   Poisson / negative binomial (integer-valued):
#     Apply log-scale Newton for continuous projection, then
#     probabilistic rounding: floor(y) with probability 1-frac(y),
#     ceil(y) with probability frac(y).  E[round(y)] = y, so the
#     score constraint holds in expectation.  Residual error is ~2.3%.
#
#   Binomial / logit (binary):
#     Score projection with a [0, 1] box.  Matching X' (y* - mu_hat) = 0
#     recovers beta_hat, its SEs and p-values exactly (dispersion fixed at 1),
#     so the copula residual is shrunk by the largest scalar that keeps y* a
#     valid fractional response in [0, 1].  A Bernoulli-rounded {0, 1} release
#     is the approximate integer-support alternative (see V5 DP-logistic).
#
#   Zero-inflated (ZI / two-part):
#     The copula bootstrap already handles zero indicators via
#     categorical strata (gen-syn-family.r).  project.glm() subsets
#     to positive observations and projects with Newton on the
#     conditional model.  A standalone project.twopart() with hurdle
#     logistic resampling is retained for reference but unused.
#
# Privacy (Remark, projection privacy):
#   The fitted component X beta_hat embeds X'y (p numbers -- same risk
#   as publishing regression output).  The residual component carries
#   n - p - 1 df of genuine novelty from the copula bootstrap.
#
# High p/n regime (Remark):
#   As p/n -> 1, H -> I, residuals -> 0, and the synthetic data
#   collapses to the fitted values.  The generator degrades gracefully
#   but offers no privacy in this limit.
#
# Numerical conditioning:
#   The inverse link can produce negative mu from extrapolated eta
#   (e.g. Gamma with inverse link when eta crosses zero).  Observations
#   with invalid mu are detected and replaced with median(mu_valid).

# +------------------------------------------------------------------------+
# | Setup                                                                  |
# +------------------------------------------------------------------------+
source("lib/00-shared.r")
source("lib/gen-syn-family.r")
# +------------------------------------------------------------------------+
# | OLS projection (Corollary to Theorem 7)                              |
# +------------------------------------------------------------------------+
# Closed-form Gaussian projection.  Decomposes y* into:
#   fitted:    X* beta_hat  (from original data -- encodes X'y)
#   residual:  c (I - H*) y_cop  (from copula -- n-p-1 df of novelty)
# where c = sqrt(phi_orig * df / RSS_cop) rescales so that RSS* = RSS.
# This satisfies both score (C1) and saturated LL (C2) exactly.
#
# QR decomposition is used for the hat matrix because the synthetic
# design X* may be rank-deficient (sparse strata can drop factor levels).
# Theorem 7 requires the (weighted) design to have
# full column rank.  A rank-deficient design is outside its hypotheses, so the
# response is a named refusal and not a ridge: regularising it would emit a
# release the theorem does not cover, quietly, and which branch was taken would
# depend on the units of the covariates.  The 1e-10 added below the call sites
# stays absolute, a numerical nicety for the admissible case.  Scaling it to
# max(abs(.)) was measured to cost the projection its machine-zero exactness on
# the binomial branch and on the zero-inflated quasi-gamma designs.
.assert.full.rank <- function(m.wX, s.where) {
    s.rk <- qr(m.wX)[["rank"]]
    if(s.rk < ncol(m.wX))
        stop(sprintf(paste0("projection (%s): the weighted design is rank ",
                            "deficient, rank %d of %d columns.  Theorem ",
                            "'sufficient statistics projection' requires full ",
                            "column rank."), s.where, s.rk, ncol(m.wX)))
    invisible(TRUE)
}

project.ols <- function(d.syn, d.orig, y.nm, x.nms) {
    fml <- reformulate(x.nms, response = y.nm)

    # fit the original model to obtain beta_hat and dispersion phi
    mo.orig <- lm(fml, d.orig)
    v.beta  <- coef(mo.orig)
    s.phi   <- sum(residuals(mo.orig)^2) /
               (nrow(d.orig) - length(v.beta))

    # synthetic design matrix and copula response
    m.Xs    <- model.matrix(fml, d.syn)
    v.y.cop <- d.syn[[y.nm]]

    # align beta_hat to synthetic X* columns.  The copula bootstrap
    # may generate strata with fewer factor levels than the original,
    # so some columns in X* may be absent from beta_hat.  Coefficients
    # for missing levels are set to zero (conservative).
    v.nms.Xs <- colnames(m.Xs)
    v.nms.b  <- names(v.beta)
    v.common <- intersect(v.nms.Xs, v.nms.b)
    v.beta.aligned <- rep(0, ncol(m.Xs))
    v.beta.aligned <- `names<-`(v.beta.aligned, v.nms.Xs)
    v.beta.aligned[v.common] <- v.beta[v.common]

    # fitted component: X* beta_hat (carries the X'y information)
    v.fit.s <- as.numeric(m.Xs %*% v.beta.aligned)

    # residual component: (I - H*) y_cop via QR factorisation.
    # H* = X* (X*'X*)^{-1} X*' is the hat matrix of the synthetic design.
    # QR handles rank deficiency (rank < ncol) gracefully.
    l.qr   <- qr(m.Xs)
    s.rank <- l.qr[["rank"]]
    v.Hy   <- qr.fitted(l.qr, v.y.cop)
    v.e.cop <- v.y.cop - v.Hy              # (I - H*) y_cop

    # scale factor c so that RSS* = phi * df  (constraint C2).
    # c = sqrt(phi_orig * df / RSS_cop), giving:
    #   y* = X* beta_hat + c * (I - H*) y_cop
    # Since (I-H*) is idempotent and orthogonal to col(X*), we have
    # X*' y* = X*' X* beta_hat => beta_hat* = beta_hat  (constraint C1).
    s.rss.cop <- sum(v.e.cop^2)
    s.df      <- nrow(m.Xs) - s.rank
    s.scale   <- `if`(s.rss.cop > 0 && s.df > 0,
                       sqrt(s.phi * s.df / s.rss.cop), 1)

    v.fit.s + v.e.cop * s.scale
}

# +------------------------------------------------------------------------+
# | Canonical link detection                                              |
# +------------------------------------------------------------------------+
# A link g is canonical for exponential family f iff V(mu) g'(mu) = 1,
# which makes the score equation X' (y - mu) = 0 linear in y.  When the
# link is canonical, the score constraint (C1) is a linear system in y*,
# so projection is a single linear solve.  For non-canonical links the
# score is nonlinear in y*, requiring iterative (Newton) projection.
is.canonical <- function(family, link) {
    switch(family,
        gaussian       = link == "identity",
        poisson        = link == "log",
        binomial       = link == "logit",
        gamma          = link == "inverse",
        `neg.binomial` = link == "log",
        `inv.gaussian` = link == "1/mu^2",
        FALSE)
}
# +------------------------------------------------------------------------+
# | GLM family object constructor                                         |
# +------------------------------------------------------------------------+
# Maps the paper's family/link string pairs to R family objects.
# The fallback for unrecognised families is quasi(log, mu^2), which
# gives quasi-likelihood with V(mu) = mu^2 (gamma-like variance).
make.glm.family <- function(family, link) {
    switch(family,
        gaussian       = gaussian(link = link),
        poisson        = poisson(link = link),
        binomial       = binomial(link = link),
        gamma          = Gamma(link = link),
        `neg.binomial` = MASS::negative.binomial(theta = 1, link = link),
        `inv.gaussian` = inverse.gaussian(link = link),
        quasi(link = link, variance = "mu^2"))
}
# +------------------------------------------------------------------------+
# | GLM projection (Theorem 7, general case)                             |
# +------------------------------------------------------------------------+
# Core projection engine.  Given copula output y_cop and the original
# data, find y* satisfying:
#   (C1) X*' w (y* - mu_hat) = 0     (score equation at beta_hat)
#   (C2) dispersion phi* = phi        (saturated LL / Pearson scaling)
#
# Three code paths, selected by family and link:
#   1. Real support (Gaussian non-identity link): weighted linear projection
#      on y-scale + RSS rescaling.
#   2. Positive support, quasi-gamma: alternating y-scale projection with
#      positivity clipping (Newton on log-scale would collapse Pearson phi).
#   3. Positive support, standard families: Newton on z = log(y) scale
#      (the log-scale Newton step).  y = exp(z) > 0 by construction.
#      Poisson/NB get probabilistic rounding afterwards.
project.glm <- function(d.syn, d.orig, y.nm, x.nms, family, link) {
    fml <- reformulate(x.nms, response = y.nm)

    # dispatch Gaussian/identity to the closed-form OLS path
    if(family == "gaussian" && link == "identity")
        return(project.ols(d.syn, d.orig, y.nm, x.nms))

    # for links that require y > 0 (log, inverse, 1/mu^2), exclude
    # zero observations from the original fit (ZI zeros are handled
    # separately by the copula's categorical bootstrap)
    b.has.zeros <- any(d.orig[[y.nm]] == 0)
    d.orig.fit <- `if`(link %in% c("log", "inverse", "1/mu^2") && b.has.zeros,
                       d.orig[d.orig[[y.nm]] > 0, ], d.orig)
    # fit the original model to extract beta_hat and phi.
    # quasi families need an explicit start (intercept = 1, slopes = 0)
    # because glm() cannot auto-initialise quasi variance functions.
    # standard families try default initialisation first, then fall back
    # to an eta-based start (intercept = g(y_bar), slopes = 0).
    mo.orig <- tryCatch({
        fam.obj <- make.glm.family(family, link)
        b.quasi <- inherits(fam.obj, "family") &&
                   fam.obj[["family"]] == "quasi"
        if(b.quasi) {
            s.np <- ncol(model.matrix(fml, d.orig.fit))
            glm(fml, d.orig.fit, family = fam.obj,
                start = c(1, rep(0, s.np - 1)),
                control = glm.control(maxit = 200))
        } else {
            tryCatch(
                glm(fml, d.orig.fit, family = fam.obj,
                    control = glm.control(maxit = 200)),
                error = function(e) {
                    s.np <- ncol(model.matrix(fml, d.orig.fit))
                    s.eta0 <- fam.obj[["linkfun"]](
                        mean(d.orig.fit[[y.nm]]))
                    glm(fml, d.orig.fit, family = fam.obj,
                        start = c(s.eta0, rep(0, s.np - 1)),
                        control = glm.control(maxit = 200))
                })
        }
    }, error = function(e) NULL)
    # if the original fit fails, return copula values unmodified
    if(is.null(mo.orig)) return(d.syn[[y.nm]])

    # extract beta_hat (p coefficients) and dispersion phi from
    # the original fit.  These are the target quantities: projection
    # constructs y* so that refitting gives these back exactly.
    v.beta     <- coef(mo.orig)
    s.phi.orig <- summary(mo.orig)[["dispersion"]]

    # build synthetic design matrix.  For ZI data with log/inverse link,
    # subset to positive observations (zeros handled by copula bootstrap).
    # v.idx.fit is reused by the final write-back, which must target
    # exactly the rows the projection ran on.
    v.idx.fit <- `if`(link %in% c("log", "inverse", "1/mu^2") && b.has.zeros,
                      which(d.syn[[y.nm]] > 0), seq_len(nrow(d.syn)))
    d.syn.fit <- d.syn[v.idx.fit, , drop = FALSE]
    m.Xs    <- model.matrix(fml, d.syn.fit)
    v.y.cop <- d.syn.fit[[y.nm]]
    s.n     <- nrow(m.Xs)
    s.p     <- ncol(m.Xs)

    # support flags: b.positive means y > 0 required (Gamma, Poisson, etc.);
    # b.integer means y must be a non-negative integer (Poisson, NB).
    b.positive <- family %in% c("gamma", "inv.gaussian", "poisson",
                                "neg.binomial") ||
                  (grepl("quasi", family) && link %in% c("log", "inverse"))
    b.integer  <- family %in% c("poisson", "neg.binomial")
    # b.binary means y in [0, 1] required (binomial).  Projected onto its own
    # [0, 1] box path below, not the real-support or positive-support paths.
    b.binary   <- family == "binomial"

    # compute mu_hat(X*) = g^{-1}(X* beta_hat) -- the fitted values
    # from the original model evaluated at synthetic covariates.
    # Fallback: if predict() fails (e.g. missing factor levels), compute
    # mu manually via aligned beta and the inverse link.
    v.mu.s <- tryCatch(
        predict(mo.orig, newdata = d.syn.fit, type = "response"),
        error = function(e) {
            v.nms.Xs <- colnames(m.Xs)
            v.nms.b  <- names(v.beta)
            v.common <- intersect(v.nms.Xs, v.nms.b)
            v.beta.a <- `names<-`(rep(0, ncol(m.Xs)), v.nms.Xs)
            v.beta.a[v.common] <- v.beta[v.common]
            family(mo.orig)[["linkinv"]](as.numeric(m.Xs %*% v.beta.a))
        })

    # Numerical conditioning: for positive-support families,
    # the inverse link can produce non-positive mu when extrapolated
    # covariates push eta outside the training range (e.g. Gamma with
    # inverse link when eta crosses zero).  Replace invalid mu with
    # the median of the valid fitted values.
    b.valid <- rep(TRUE, s.n)
    if(b.positive) {
        b.valid <- is.finite(v.mu.s) & v.mu.s > 1e-6
        if(!all(b.valid)) {
            v.mu.s[!b.valid] <- median(v.mu.s[b.valid])
        }
    }
    # linear predictor eta = g(mu) for score weight computation
    v.eta.s <- family(mo.orig)[["linkfun"]](pmax(v.mu.s, 1e-10))

    # score weights w_i = (dmu/deta)_i / V(mu_i), evaluated at beta_hat.
    # The GLM score equation is X' diag(w) (y - mu) = 0.
    # For canonical links w_i simplifies: e.g. Poisson(log) gives w_i = 1,
    # Gamma(inverse) gives w_i = 1/mu^2.
    v.dmu.deta <- family(mo.orig)[["mu.eta"]](v.eta.s)
    v.V        <- family(mo.orig)[["variance"]](pmax(v.mu.s, 1e-10))
    v.w        <- v.dmu.deta / v.V

    # quasi families (V(mu) = mu^2 with log link) give w_i = 1/mu_i and
    # score = X'(y-mu)/mu.  For these, Newton on log-scale collapses
    # Pearson dispersion, so we use y-scale alternating projection instead.
    # Standard families (Gamma, Poisson, NB) use log-scale Newton, which
    # respects the support constraint y > 0 automatically.
    b.quasi.gamma <- grepl("quasi", family)

    # ---- path B: Binomial (project onto [0, 1] with exact score match) ----
    # For binomial, matching the score X*' W (y* - mu_hat) = 0 recovers
    # beta_hat exactly, and the standard errors follow from the Fisher
    # information X*' W X* alone, so no dispersion rescaling is needed
    # (phi = 1 is fixed).  The response must stay in [0, 1] (a valid binomial
    # proportion).  One weighted score correction zeroes the score, leaving a
    # score-orthogonal residual r.  Scaling r by any scalar keeps the score at
    # zero, so shrink it by the largest s in [0, 1] that keeps mu + s r inside
    # [0, 1] element-wise.  The result is an exact, valid fractional response;
    # a Bernoulli-rounded {0, 1} release would trade this exactness for integer
    # support (the irreducible sampling variance of the V5 DP-logistic mode).
    if(b.binary) {
        m.wX    <- sweep(m.Xs, 1, v.w, `*`)
        .assert.full.rank(m.wX, "binomial")
        m.XtwX  <- crossprod(m.wX)
        diag(m.XtwX) <- diag(m.XtwX) + 1e-10
        v.score <- as.numeric(crossprod(m.wX, v.y.cop - v.mu.s))
        v.corr  <- as.numeric(m.wX %*% solve(m.XtwX, v.score))
        v.r     <- (v.y.cop - v.corr) - v.mu.s        # score-orthogonal residual
        # largest s in [0, 1] keeping mu + s r within [0, 1] element-wise
        v.s.up  <- ifelse(v.r >  1e-12, (1 - v.mu.s) / v.r, Inf)
        v.s.lo  <- ifelse(v.r < -1e-12, (0 - v.mu.s) / v.r, Inf)
        s.shr   <- max(0, min(1, min(v.s.up), min(v.s.lo)))
        v.y.proj <- pmin(pmax(v.mu.s + s.shr * v.r, 0), 1)
    # ---- path 1: Real support (Gaussian with non-identity link) ----------
    # Score constraint X*' W (y* - mu) = 0 is linear in y*, so one
    # weighted least-squares correction suffices.  Then rescale
    # residuals to match RSS (constraint C2: saturated LL).
    } else if(!b.positive) {
        # weighted design matrix and score correction
        m.wX    <- sweep(m.Xs, 1, v.w, `*`)
        .assert.full.rank(m.wX, "non-canonical")
        m.XtwX  <- crossprod(m.wX)
        diag(m.XtwX) <- diag(m.XtwX) + 1e-10  # Tikhonov regularisation
        v.score <- as.numeric(crossprod(m.wX, v.y.cop - v.mu.s))
        # minimum-norm correction: project score residual back onto col(wX)
        v.corr  <- as.numeric(m.wX %*% solve(m.XtwX, v.score))
        v.y.proj <- v.y.cop - v.corr

        # constraint C2: rescale residuals so RSS* = RSS_orig
        v.e.proj   <- v.y.proj - v.mu.s
        s.rss.orig <- sum(residuals(mo.orig)^2)
        s.rss.proj <- sum(v.e.proj^2)
        if(s.rss.proj > 0) {
            s.scale  <- sqrt(s.rss.orig / s.rss.proj)
            v.y.proj <- v.mu.s + v.e.proj * s.scale
        }
    # ---- path 2: Quasi-gamma (y-scale alternating projection) -----------
    # For quasi(log, mu^2), the score X' diag(1/mu) (y - mu) = 0 is
    # linear in y but positivity y > 0 is not guaranteed after a single
    # projection step.  Alternate: (a) project score to zero, (b) clip
    # negatives to 1e-6, until all y > 0.  This preserves Pearson
    # dispersion better than Newton on log-scale for quasi families.
    } else if(b.quasi.gamma) {
        m.wX    <- sweep(m.Xs, 1, v.w, `*`)
        .assert.full.rank(m.wX, "quasi-gamma")
        m.XtwX  <- crossprod(m.wX)
        diag(m.XtwX) <- diag(m.XtwX) + 1e-10
        m.XtwX.inv <- solve(m.XtwX)
        v.y.cur <- v.y.cop

        s.maxiter <- 50
        s.iter    <- 0
        repeat {
            s.iter <- s.iter + 1
            # score correction step: project onto C1
            v.score <- as.numeric(crossprod(m.wX, v.y.cur - v.mu.s))
            v.corr  <- as.numeric(m.wX %*% (m.XtwX.inv %*% v.score))
            v.y.cur <- v.y.cur - v.corr

            # positivity enforcement: clip negatives, then re-check score
            s.n.neg <- sum(v.y.cur <= 0)
            if(s.n.neg == 0 || s.iter >= s.maxiter) break
            v.y.cur <- pmax(v.y.cur, 1e-6)
        }

        # constraint C2: rescale Pearson residuals to match original phi.
        # Pearson residual r_i = (y_i - mu_i) / sqrt(V(mu_i)).
        # phi_P = sum(r^2) / (n-p), so scaling r by sqrt(phi_orig/phi_cur)
        # enforces phi* = phi exactly.
        v.pearson  <- (v.y.cur - v.mu.s) /
                      sqrt(family(mo.orig)[["variance"]](pmax(v.mu.s, 1e-10)))
        s.phi.cur  <- sum(v.pearson^2) / max(s.n - s.p, 1)
        if(s.phi.cur > 1e-10) {
            s.c <- sqrt(s.phi.orig / s.phi.cur)
            v.y.cur <- v.mu.s + (v.y.cur - v.mu.s) * s.c
            v.y.cur <- pmax(v.y.cur, 1e-6)
        }
        v.y.proj <- v.y.cur
    # ---- path 3: Positive support (log-scale Newton iteration) ----------
    # The log-scale Newton step.  Work on z = log(y), solve
    #   F(z) := X*' w (exp(z) - mu_hat) = 0
    # by Newton's method.  Since y = exp(z) > 0 for all real z, the
    # support constraint is satisfied automatically at every iterate.
    #
    # Jacobian: J = X*' diag(w * exp(z)) = X*' diag(w * y).
    # Minimum-norm update: dz = diag(wy) X (X' diag(wy)^2 X)^{-1} score.
    # Because dz lies in col(diag(wy) X), the component of z orthogonal
    # to this column space (the copula residuals on log-scale) is
    # untouched.  This preserves n-p-1 df of copula novelty.
    #
    # Convergence is fast (typically 3-5 iterations) because the copula
    # output is already close to the constraint surface.
    } else {
        v.z <- log(pmax(v.y.cop, 1e-10))

        s.maxiter <- 30
        s.tol     <- 1e-8
        s.iter    <- 0
        repeat {
            s.iter <- s.iter + 1
            v.y <- exp(v.z)

            # evaluate the score: X*' w (y - mu_hat)
            v.score <- as.numeric(crossprod(m.Xs, v.w * (v.y - v.mu.s)))
            if(max(abs(v.score)) < s.tol || s.iter > s.maxiter) break

            # Gauss-Newton step: solve J J' lambda = score for the
            # minimum-norm update dz = J' lambda
            v.wy  <- v.w * v.y
            m.WYX <- sweep(m.Xs, 1, v.wy, `*`)
            m.JJt <- crossprod(m.WYX)
            diag(m.JJt) <- diag(m.JJt) * (1 + 1e-10)  # regularise

            v.lambda <- tryCatch(solve(m.JJt, v.score),
                                 error = function(e) rep(0, s.p))
            v.dz <- as.numeric(m.WYX %*% v.lambda)

            # damped step: cap max |dz| at 3 to prevent overshooting
            s.step <- min(1, 3 / (max(abs(v.dz)) + 1e-16))
            v.z <- v.z - s.step * v.dz
        }

        # constraint C2: rescale Pearson residuals to match phi_orig.
        # Pearson residual: r_i = (y_i - mu_i) / sqrt(V(mu_i)).
        # After rescaling: y*_i = mu_i + sqrt(V(mu_i)) * r_i * c,
        # where c = sqrt(phi_orig / phi_newton).
        v.y.newton <- exp(v.z)
        v.sd.s     <- sqrt(family(mo.orig)[["variance"]](pmax(v.mu.s, 1e-10)))
        v.pearson  <- (v.y.newton - v.mu.s) / v.sd.s
        s.phi.newton <- sum(v.pearson^2) / max(s.n - s.p, 1)
        if(s.phi.newton > 1e-10) {
            s.c <- sqrt(s.phi.orig / s.phi.newton)
            v.y.proj <- v.mu.s + v.sd.s * v.pearson * s.c
            # enforce y > 0 after Pearson rescaling (rescaling can push
            # observations with small mu below zero)
            if(b.positive) v.y.proj <- pmax(v.y.proj, 1e-10)
        } else {
            v.y.proj <- v.y.newton
        }

        # Poisson / NB: probabilistic rounding to integers.
        # P(ceil) = frac(y), P(floor) = 1-frac(y), so E[round(y)] = y.
        # This introduces ~2.3% residual error in the score constraint
        # (unavoidable for integer-valued families).
        if(b.integer) {
            v.y.cont <- v.y.proj
            v.frac <- v.y.cont - floor(v.y.cont)
            set.seed(42)
            v.y.proj <- pmax(ifelse(runif(s.n) < v.frac,
                                    ceiling(v.y.cont),
                                    floor(v.y.cont)), 0)
        }

        # for observations with invalid mu (numerical conditioning),
        # keep the original copula values unchanged
        if(!all(b.valid))
            v.y.proj[!b.valid] <- v.y.cop[!b.valid]
    }

    # write projected values back into the full synthetic vector, onto
    # the rows the projection ran on (v.idx.fit).  For ZI data the rows
    # outside v.idx.fit keep their zeros from the copula's categorical
    # bootstrap.
    v.result <- d.syn[[y.nm]]
    v.result[v.idx.fit] <- v.y.proj
    v.result
}
# +------------------------------------------------------------------------+
# | ZI projection (conditional on positive)                               |
# +------------------------------------------------------------------------+
# For zero-inflated data, the copula bootstrap preserves the zero
# structure via categorical strata (gen-syn-family.r).  This function
# projects only the positive-subset values onto the conditional model's
# score constraint.  It duplicates logic now handled inside project.glm()
# and is retained for reference but currently unused by the main pipeline.
project.zi <- function(d.syn, d.orig, y.nm, x.nms, family, link) {
    v.y.orig <- d.orig[[y.nm]]
    v.y.syn  <- d.syn[[y.nm]]

    b.zero.syn  <- v.y.syn == 0
    b.zero.orig <- v.y.orig == 0

    # degenerate: nothing to project if all zeros on either side
    if(all(b.zero.syn) || all(b.zero.orig))
        return(v.y.syn)

    # subset to positive observations for conditional projection
    v.idx.pos.syn  <- which(!b.zero.syn)
    v.idx.pos.orig <- which(!b.zero.orig)

    d.orig.pos <- d.orig[v.idx.pos.orig, , drop = FALSE]
    d.syn.pos  <- d.syn[v.idx.pos.syn, , drop = FALSE]

    # fit the conditional (positive-only) model on original data
    fml <- reformulate(x.nms, response = y.nm)
    mo.pos <- tryCatch({
        s.np <- length(x.nms) + 1
        `if`(family == "poisson",
             glm(fml, d.orig.pos, family = poisson()),
             glm(fml, d.orig.pos,
                 family = quasi(link = link, variance = "mu^2"),
                 start = c(1, rep(0, s.np - 1)),
                 control = glm.control(maxit = 100)))
    }, error = function(e) NULL)
    if(is.null(mo.pos)) return(v.y.syn)

    v.beta.pos <- coef(mo.pos)
    s.phi.pos  <- summary(mo.pos)[["dispersion"]]

    # design matrices for synthetic and original positive subsets
    m.X.syn <- model.matrix(fml, d.syn.pos)
    m.X.orig <- model.matrix(fml, d.orig.pos)

    # fitted values from the conditional model at synthetic covariates
    v.eta.syn  <- as.numeric(m.X.syn %*% v.beta.pos)
    v.mu.syn   <- `if`(link == "log", exp(v.eta.syn), v.eta.syn)

    # decompose into fitted + residual (the projection structure)
    v.mu.orig.pos <- fitted(mo.pos)
    v.e.orig <- d.orig.pos[[y.nm]] - v.mu.orig.pos
    v.e.syn  <- d.syn.pos[[y.nm]] - v.mu.syn

    # variance function V(mu) for Pearson residuals:
    # log link => V(mu)=mu^2 for quasi-gamma, V(mu)=mu for Poisson
    if(link == "log") {
        v.V.orig <- v.mu.orig.pos^2
        v.V.syn  <- v.mu.syn^2
        if(family == "poisson") {
            v.V.orig <- v.mu.orig.pos
            v.V.syn  <- v.mu.syn
        }
    } else {
        v.V.orig <- rep(1, length(v.mu.orig.pos))
        v.V.syn  <- rep(1, length(v.mu.syn))
    }

    # rescale copula residuals to match original dispersion (C2)
    s.phi.syn <- `if`(length(v.e.syn) > length(v.beta.pos),
                       sum(v.e.syn^2 / v.V.syn) /
                           (length(v.e.syn) - length(v.beta.pos)),
                       s.phi.pos)
    s.scale <- `if`(s.phi.syn > 0, sqrt(s.phi.pos / s.phi.syn), 1)

    # assemble: original fitted + scaled novel residuals
    v.y.proj.pos <- v.mu.syn + v.e.syn * s.scale

    # enforce y > 0 for positive-support families
    v.y.proj.pos <- pmax(v.y.proj.pos, .Machine[["double.eps"]])

    v.y.syn[v.idx.pos.syn] <- v.y.proj.pos
    v.y.syn
}
# +------------------------------------------------------------------------+
# | Two-part (hurdle) projection                                          |
# +------------------------------------------------------------------------+
# Alternative ZI handler using explicit hurdle logistic resampling.
# Part 1: fit logistic regression I(Y>0) ~ X on original data, then
#   assign synthetic indicators by ranking fitted probabilities pi_hat(X*)
#   and selecting the top n_pos (deterministic to preserve zero proportion).
# Part 2: project positive values via project.glm() on the conditional model.
#
# This function is retained for reference but unused -- the main pipeline
# uses project.glm() directly, which subsets to positives and projects
# with Newton.  The copula's categorical bootstrap already handles zero
# indicators via strata.
project.twopart <- function(d.syn, d.orig, y.nm, x.nms,
                            family, link) {
    # Part 1: Hurdle -- project zero/non-zero indicators
    v.y.orig <- d.orig[[y.nm]]
    v.y.syn  <- d.syn[[y.nm]]
    v.ind.orig <- as.integer(v.y.orig > 0)

    # fit logistic hurdle: I(Y>0) ~ X on original
    fml <- reformulate(x.nms, response = "ind")
    d.hurdle.orig <- cbind(d.orig, ind = v.ind.orig)
    mo.hurdle <- tryCatch(
        glm(fml, d.hurdle.orig, family = binomial(link = "logit")),
        error = function(e) NULL)

    if(!is.null(mo.hurdle)) {
        v.alpha <- coef(mo.hurdle)
        m.Xs <- model.matrix(fml, cbind(d.syn, ind = 0))
        v.pi.hat <- binomial()[["linkinv"]](as.numeric(m.Xs %*% v.alpha))

        # deterministic assignment: rank by pi_hat, top-n_pos become
        # positive.  This preserves zero proportion exactly.
        s.n.pos.orig <- sum(v.ind.orig)
        v.ind.syn <- rep(0L, nrow(d.syn))
        v.ind.syn[order(v.pi.hat, decreasing = TRUE)[
            seq_len(s.n.pos.orig)]] <- 1L
    } else {
        # fallback: keep copula's zero/non-zero pattern
        v.ind.syn <- as.integer(v.y.syn > 0)
    }

    # Part 2: project positive values via conditional GLM
    v.idx.pos.orig <- which(v.y.orig > 0)
    v.idx.pos.syn  <- which(v.ind.syn == 1)

    # degenerate: fewer than 2 positives on either side
    if(length(v.idx.pos.syn) < 2 || length(v.idx.pos.orig) < 2) {
        v.y.syn[v.ind.syn == 0] <- 0
        return(v.y.syn)
    }

    d.orig.pos <- d.orig[v.idx.pos.orig, , drop = FALSE]
    d.syn.pos  <- d.syn[v.idx.pos.syn, , drop = FALSE]

    # ensure copula values designated as positive are actually > 0
    # (the bootstrap may assign zero to a "positive" slot, so replace
    # with the 1st percentile of the original positive distribution)
    if(any(d.syn.pos[[y.nm]] <= 0)) {
        v.pos.orig <- d.orig.pos[[y.nm]]
        s.small    <- quantile(v.pos.orig, 0.01)
        d.syn.pos[[y.nm]] <- pmax(d.syn.pos[[y.nm]], s.small)
    }

    # project positive values using the conditional GLM
    v.y.pos <- project.glm(d.syn.pos, d.orig.pos, y.nm, x.nms,
                           family, link)

    # assemble: zeros from hurdle + projected positives
    v.result <- rep(0, nrow(d.syn))
    v.result[v.idx.pos.syn] <- v.y.pos
    v.result
}
# +------------------------------------------------------------------------+
# | Combined generator: copula + projection (main entry point)            |
# +------------------------------------------------------------------------+
# Two-step procedure:
#   Step 1: generate a copula-based synthetic dataset (gen-syn-family.r).
#           This produces y_cop that approximately preserves the marginal
#           and joint distributions.  Categorical variables (including ZI
#           zero indicators) are reproduced via stratified bootstrap.
#   Step 2: project y_cop onto the sufficient statistics constraint
#           surface (C1)+(C2).  After projection, y* satisfies Theorem 7
#           and all GLM inference quantities are exact.
#
# The response variable y.nm is projected.  All predictors remain as
# generated by the copula (their joint distribution with y is preserved
# by the copula, and the projection only modifies y).
gen.syn.projected <- function(d, y.nm = "response", family = "gaussian",
                              link = "identity", supports = NULL,
                              order = 4, seed = 6489,
                              design.exact = FALSE) {
    # Step 1: copula-based starting point.  design.exact draws the stratum
    # labels as a permutation (exact cell counts, condition (i) of the
    # design-Gram proposition).
    d.syn <- gen.syn.family(d, order = order, supports = supports,
                            shrink = TRUE, seed = seed,
                            strata.exact = design.exact)

    # identify predictor columns (everything except the response)
    v.types <- mapply(id.cat, d)
    v.x.nms <- setdiff(names(d), y.nm)

    # Step 1b: pin the design Gram before the y-projection, so the
    # refit information matches and the full triple is exact
    # (the design-Gram proposition).
    if(design.exact) {
        d.syn <- project.designgram(d.syn, d, y.nm, v.x.nms,
                                    family = family, link = link)
    }

    # Step 2: project onto the sufficient statistics constraint.
    # For ZI data, project.glm() automatically subsets to positive
    # observations and projects with Newton on the conditional model.
    # The copula's categorical bootstrap already handles zero indicators.
    v.y.proj <- project.glm(d.syn, d, y.nm, v.x.nms, family, link)

    d.syn[[y.nm]] <- v.y.proj
    d.syn
}
# +------------------------------------------------------------------------+
# | Two-part (hurdle) projection generator (zero-pattern preserving)       |
# +------------------------------------------------------------------------+
# gen.syn.projected above draws a fresh zero pattern from the copula, so for a
# zero-inflated (hurdle) response project.glm() refits the positive block on
# the synthetic positive design.  The coefficients still recover, but the
# positive-block information X_syn+' W X_syn+ differs from the original, so the
# standard errors do not (the 3.7% and 4.3% of the fixed-design standard-error
# table on Sc.4 and Sc.8).
#
# The hurdle likelihood separates into a zero block (a GLM for 1{y>0} on the
# full design X) and a positive block (a GLM on the positive rows, design X_+),
# and the two are variation-independent, so the information is block-diagonal.
# This generator preserves the original zero pattern and projects the positive
# component onto the original positive design.  The zero block is then exact by
# identity (same indicator, same X), and the positive block is an ordinary
# sufficient-statistics projection on the shared positive design, exact in
# coefficients and standard errors.  This is the "Proposition (two-part
# projection)" of the supplement, demonstrated by 18c-twopart-se.r.
#
# The design is released unchanged (covariates treated as public, as in the
# fixed-design demonstration 18b), so this is an opt-in mode for hurdle data,
# and gen.syn.projected remains the default.  family and link describe the
# positive block and must use a positive-support link (log, inverse, or
# 1/mu^2) so that project.glm() subsets to the positive rows.  The postcondition
# below fails loudly otherwise.
gen.syn.projected.twopart <- function(d, y.nm = "response", family = "gaussian",
                                      link = "identity", supports = NULL,
                                      order = 4, seed = 6489,
                                      design.exact = FALSE) {
    # The two-part full-synthesis corollary.  Full synthesis with the
    # design pinned per block.  The copula's joined stratification (the .pos
    # indicator) under strata.exact gives exact per-cell positive counts, so the
    # synthetic zero pattern is a permutation of the original one, and the
    # disjoint block recolouring makes the positive-block and all-rows Grams
    # exact.
    if(design.exact) {
        d.syn <- gen.syn.family(d, order = order, supports = supports,
                                shrink = TRUE, seed = seed,
                                strata.exact = TRUE)
        v.x.nms <- setdiff(names(d), y.nm)
        b.pos.syn <- d.syn[[y.nm]] > 0
        d.syn <- .designgram.twopart(d.syn, d, y.nm, v.x.nms)

        # positive block: ordinary y-projection on the pinned positive
        # design (project.glm subsets to the synthetic positive rows)
        v.y.proj <- project.glm(d.syn, d, y.nm, v.x.nms, family, link)
        d.syn[[y.nm]] <- v.y.proj
        stopifnot(all((d.syn[[y.nm]] > 0) == b.pos.syn))
        return(d.syn)
    }

    # original zero/positive partition, released unchanged
    b.pos <- d[[y.nm]] > 0

    # Step 1. Copula starting point.  Only the positive rows carry its novelty.
    d.syn <- gen.syn.family(d, order = order, supports = supports,
                            shrink = TRUE, seed = seed)

    # Step 2. Two-part mixture on the original design.  Inject copula novelty on
    # the positive rows, falling back to the original value wherever the copula
    # drew a non-positive (which would break the hurdle support), and force the
    # zero rows to exactly zero.
    d.mix <- d
    v.pos.cop <- d.syn[[y.nm]][b.pos]
    v.pos.cop[v.pos.cop <= 0] <- d[[y.nm]][b.pos][v.pos.cop <= 0]
    d.mix[[y.nm]][!b.pos] <- 0
    d.mix[[y.nm]][b.pos]  <- v.pos.cop

    # Step 3. Project the positive component onto the original positive design.
    # project.glm() subsets to the positive rows internally, and d.mix carries the
    # original zero pattern, so those are exactly the original positive rows.
    v.x.nms  <- setdiff(names(d), y.nm)
    v.y.proj <- project.glm(d.mix, d, y.nm, v.x.nms, family, link)

    d.out <- d
    d.out[[y.nm]] <- v.y.proj
    stopifnot(all((d.out[[y.nm]] > 0) == b.pos))   # zero pattern preserved
    d.out
}
# +------------------------------------------------------------------------+
# | Direction-agnostic projection (OLS cross-product matching)            |
# +------------------------------------------------------------------------+
# When Y is unknown (Scenario 1: any variable could be the response),
# we need y* such that every OLS regression among the variables gives
# exact beta_hat.  For OLS, beta_hat depends only on the cross-product
# matrix Z'Z where Z = [1, D, C] (intercept, dummies, continuous).
#
# The projection leaves D (categorical columns) untouched and moves only
# the continuous columns C so that:
#   (C1) D' C* = D_o' C_o     (dummy-continuous cross-products)
#   (C2) C*' C* = C_o' C_o    (continuous-continuous covariance)
#
# Precondition: the full Gram also needs D'D = D_o'D_o, which this
# function cannot enforce (it does not move D).  The caller must supply a
# base draw whose categorical cell counts equal the original's, i.e. the
# strata.exact permutation draw of gen.syn.family (condition (i) of the
# design-Gram projection proposition).  gen.syn.projected.symmetric does
# so.  On a free bootstrap draw, (C1) forces original column sums onto
# perturbed cell counts, distorting per-stratum means and making
# dummy-involved directions worse than the unprojected copula.
#
# Solution: two-step projection.
#   Step 1: project C* so D' C* = target, using H_D = D(D'D)^{-1}D'.
#           Fitted part = D (D'D)^{-1} D_o' C_o.  Residual = (I-H_D) C_cop.
#   Step 2: Cholesky-rotate the residual so C*'C* = target.
#           Since (I-H_D) is orthogonal to col(D), the rotation does not
#           affect D'C*, so both constraints are satisfied simultaneously.
project.crossproduct <- function(d.syn, d.orig) {
    v.types <- mapply(id.cat, d.orig)
    v.cnt.nm <- names(v.types[v.types == "cnt"])
    if(length(v.cnt.nm) == 0) return(d.syn)

    # partition the full model matrix Z = [1, D, C] into:
    #   D: intercept + dummy columns (fixed by bootstrap)
    #   C: continuous columns (free to project)
    fml <- reformulate(names(d.orig))
    m.Z.orig <- model.matrix(fml, d.orig)
    m.Z.syn  <- model.matrix(fml, d.syn)
    v.cnt.idx <- which(colnames(m.Z.orig) %in% v.cnt.nm)
    v.cat.idx <- setdiff(seq_len(ncol(m.Z.orig)), v.cnt.idx)
    s.p.cnt <- length(v.cnt.idx)

    m.D.syn  <- m.Z.syn[, v.cat.idx, drop = FALSE]
    m.D.orig <- m.Z.orig[, v.cat.idx, drop = FALSE]
    m.C.syn  <- m.Z.syn[, v.cnt.idx, drop = FALSE]
    m.C.orig <- m.Z.orig[, v.cnt.idx, drop = FALSE]

    # Step 1: project C* so D' C* = D_o' C_o (constraint C1).
    # H_D = D_syn (D_syn'D_syn)^{-1} D_syn' is the hat matrix for D_syn.
    # Fitted part of C* comes from the original cross-product target;
    # residual part comes from the copula (untouched by projection).
    m.DtD.syn <- crossprod(m.D.syn)
    m.DtD.inv <- tryCatch(solve(m.DtD.syn), error = function(e) NULL)
    if(is.null(m.DtD.inv)) return(d.syn)

    m.H.D <- m.D.syn %*% m.DtD.inv %*% t(m.D.syn)
    m.target.DtC <- crossprod(m.D.orig, m.C.orig)

    m.fitted <- m.D.syn %*% m.DtD.inv %*% m.target.DtC
    m.resid  <- m.C.syn - m.H.D %*% m.C.syn
    m.C.step1 <- m.fitted + m.resid

    # Step 2: Cholesky rotation of residuals so C*'C* = C_o'C_o (C2).
    # The residuals lie in null(D'), so rotating them preserves D'C*.
    # Target residual cross-product = C_o'C_o - fitted'fitted.
    # Rotation: R_new = R_old L_old^{-1} L_target, where L = chol(R'R).
    m.resid.orig <- m.C.orig - m.H.D %*% m.C.orig
    m.CtC.target <- crossprod(m.C.orig) - crossprod(m.fitted)
    m.CtC.resid  <- crossprod(m.resid)

    m.L.target <- tryCatch(chol(m.CtC.target), error = function(e) NULL)
    m.L.resid  <- tryCatch(chol(m.CtC.resid), error = function(e) NULL)

    if(!is.null(m.L.target) && !is.null(m.L.resid)) {
        m.resid.new <- m.resid %*% solve(m.L.resid) %*% m.L.target
        m.C.final <- m.fitted + m.resid.new
    } else {
        # fallback: skip Cholesky rotation if either matrix is singular
        m.C.final <- m.C.step1
    }

    # write projected continuous columns back into the data frame
    Reduce(function(d.acc, j) {
        d.acc[[v.cnt.nm[j]]] <- m.C.final[, j]
        d.acc
    }, seq_len(s.p.cnt), init = d.syn)
}

# +------------------------------------------------------------------------+
# | Design-Gram projection: pin X*'X* = X'X for a designated Y            |
# +------------------------------------------------------------------------+
# Generalises project.crossproduct to the designated-response setting: the
# response column is left to the y-projection, and the predictor block
# X = [D C] (intercept + dummies, continuous) is pinned so X*'X* = X'X
# (the design-Gram proposition).  Step 1 pins D'C (per-stratum sums), step 2
# Cholesky-rotates the D-orthogonal residual so C'C matches.  The stratum
# permutation of gen.syn.family(strata.exact = TRUE) has already pinned D'D.
#
# family/link != gaussian/identity pins the Fisher-weighted Gram
# X*' W(X* beta) X* instead, by fixed point.  Clause (b) of the design-Gram
# proposition governs that pin for a canonical link, the Poisson channel among
# the designs here.  Its response condition is the canonical score identity
# X*'y* = X*'mu, which zeroes the score only then.  The gamma and quasi-gamma
# log channels sit outside the clause, and project.glm reaches them by a
# Stronger condition than it states, zeroing the general weighted score
# X*' diag((dmu/deta)/V(mu)) (y* - mu*).  What costs exactness there is the
# alternation with the positivity clip, and on Poisson the response rounding.
# The information weights (dmu/deta)^2 / V(mu) are not project.glm()'s score
# weights (dmu/deta) / V(mu): the former enter summary.glm's covariance.
#
# Every failure stops loudly.  A silent fallback here is the project's
# error: the caller would carry an identity nothing delivers.

# One weighted two-step recolouring pass: pin (sqrt(w)D)'(sqrt(w)C) and
# (sqrt(w)C)'(sqrt(w)C) to the original weighted targets.  w = 1 gives the
# unweighted OLS case.  chol() returns upper-triangular L with L'L = A, so
# the rotation M = L.res^{-1} L.tgt right-multiplies the residual block.
.designgram.recolour <- function(m.D.syn, m.D.orig, m.C.orig, m.C.cur, v.w,
                                 v.w.orig = v.w) {
    v.rw.s <- sqrt(v.w)
    v.rw.o <- sqrt(v.w.orig)
    m.Dw.s <- m.D.syn  * v.rw.s
    m.Dw.o <- m.D.orig * v.rw.o
    m.Cw.o <- m.C.orig * v.rw.o
    m.Cw.c <- m.C.cur  * v.rw.s

    m.DtD <- crossprod(m.Dw.s)
    m.DtD.inv <- tryCatch(solve(m.DtD),
                          error = function(e)
                              stop("design-gram: weighted D'D singular"))
    m.tgt.DtC <- crossprod(m.Dw.o, m.Cw.o)
    m.fit  <- m.Dw.s %*% (m.DtD.inv %*% m.tgt.DtC)
    m.res  <- m.Cw.c - m.Dw.s %*% (m.DtD.inv %*% crossprod(m.Dw.s, m.Cw.c))

    m.G.tgt <- crossprod(m.Cw.o) - crossprod(m.fit)
    m.G.res <- crossprod(m.res)
    m.L.tgt <- tryCatch(chol(m.G.tgt),
                        error = function(e)
                            stop("design-gram: target within Gram not PD"))
    m.L.res <- tryCatch(chol(m.G.res),
                        error = function(e)
                            stop("design-gram: residual Gram not PD"))
    m.Cw.new <- m.fit + m.res %*% solve(m.L.res, m.L.tgt)
    m.Cw.new / v.rw.s
}

project.designgram <- function(d.syn, d.orig, y.nm, x.nms,
                               family = "gaussian", link = "identity",
                               tol = 1e-8, maxiter = 100L) {
    fml <- reformulate(x.nms)
    m.X.orig <- model.matrix(fml, d.orig)
    m.X.syn  <- model.matrix(fml, d.syn)
    if(!identical(colnames(m.X.orig), colnames(m.X.syn)))
        stop("design-gram: original and synthetic model matrices disagree")

    v.types   <- mapply(id.cat, d.orig[x.nms])
    v.cnt.nm  <- names(v.types[v.types == "cnt"])
    if(length(v.cnt.nm) == 0) return(d.syn)
    v.cnt.idx <- which(colnames(m.X.orig) %in% v.cnt.nm)
    v.cat.idx <- setdiff(seq_len(ncol(m.X.orig)), v.cnt.idx)

    m.D.syn  <- m.X.syn[, v.cat.idx, drop = FALSE]
    m.D.orig <- m.X.orig[, v.cat.idx, drop = FALSE]
    m.C.orig <- m.X.orig[, v.cnt.idx, drop = FALSE]
    m.C.cur  <- m.X.syn[, v.cnt.idx, drop = FALSE]

    # the permutation draw must have pinned the counts already
    if(max(abs(crossprod(m.D.syn) - crossprod(m.D.orig))) > 1e-9)
        stop("design-gram: D'D not pinned; generate with strata.exact = TRUE")

    b.ols <- family == "gaussian" && link == "identity"
    if(b.ols) {
        m.C.new <- .designgram.recolour(m.D.syn, m.D.orig, m.C.orig, m.C.cur,
                                        v.w = rep(1, nrow(m.D.syn)))
    } else {
        # Fisher-weighted fixed point (clause (b) of the design-Gram
        # proposition, which is stated for a canonical link).  The exact
        # solution exists (a row permutation of the original attains it).  The
        # loop starts at the copula draw and stops when the full weighted Gram,
        # including the D'WD block that only the weights can move, is within tol
        # relative deviation.  The original fit reuses the project.glm() idiom
        # (quasi families need an explicit start).
        fml.y <- reformulate(x.nms, response = y.nm)
        fam.obj <- make.glm.family(family, link)
        b.quasi <- inherits(fam.obj, "family") &&
                   fam.obj[["family"]] == "quasi"
        d.orig.fit <- `if`(link %in% c("log", "inverse", "1/mu^2") &&
                           any(d.orig[[y.nm]] == 0),
                           d.orig[d.orig[[y.nm]] > 0, ], d.orig)
        mo.orig <- `if`(b.quasi,
            glm(fml.y, d.orig.fit, family = fam.obj,
                start = c(1, rep(0, ncol(m.X.orig) - 1)),
                control = glm.control(maxit = 200)),
            glm(fml.y, d.orig.fit, family = fam.obj,
                control = glm.control(maxit = 200)))
        v.beta <- coef(mo.orig)[colnames(m.X.orig)]
        if(any(is.na(v.beta)))
            stop("design-gram: original fit rank-deficient on x.nms")
        w.info <- function(m.X) {
            v.eta <- as.numeric(m.X %*% v.beta)
            fam.obj[["mu.eta"]](v.eta)^2 /
                fam.obj[["variance"]](fam.obj[["linkinv"]](v.eta))
        }
        v.w.orig <- w.info(m.X.orig)
        m.I.tgt  <- crossprod(m.X.orig * sqrt(v.w.orig))
        s.scl.w  <- max(abs(m.I.tgt))

        s.iter <- 0L
        repeat {
            s.iter <- s.iter + 1L
            m.X.cur <- m.X.syn
            m.X.cur[, v.cnt.idx] <- m.C.cur
            v.w.cur <- w.info(m.X.cur)
            m.I.cur <- crossprod(m.X.cur * sqrt(v.w.cur))
            if(max(abs(m.I.cur - m.I.tgt)) < tol * s.scl.w) break
            if(s.iter > maxiter)
                stop("design-gram: weighted fixed point failed to converge")
            m.C.cur <- .designgram.recolour(m.D.syn, m.D.orig, m.C.orig,
                                            m.C.cur, v.w = v.w.cur,
                                            v.w.orig = v.w.orig)
        }
        m.C.new <- m.C.cur
    }

    d.out <- Reduce(function(d.acc, j) {
        d.acc[[v.cnt.nm[j]]] <- m.C.new[, j]
        d.acc
    }, seq_along(v.cnt.nm), init = d.syn)

    # postcondition: the identity the caller relies on, at matrix scale.
    # OLS asserts the raw Gram.  The weighted branch asserts the Fisher
    # Gram it pinned, at its own stopping tolerance.
    m.X.chk <- model.matrix(fml, d.out)
    if(b.ols) {
        s.scl <- max(abs(crossprod(m.X.orig)))
        if(max(abs(crossprod(m.X.chk) - crossprod(m.X.orig))) > 1e-8 * s.scl)
            stop("design-gram: postcondition X*'X* = X'X violated")
    } else {
        m.I.chk <- crossprod(m.X.chk * sqrt(w.info(m.X.chk)))
        if(max(abs(m.I.chk - m.I.tgt)) > tol * s.scl.w)
            stop("design-gram: postcondition weighted Gram violated")
    }
    d.out
}

# Two-part design pinning (the two-part full-synthesis corollary): disjoint row
# blocks.  The positive rows are recoloured to the original positive-block
# targets and frozen, so the positive-model Gram X_+'X_+ is exact and the
# subsequent y-projection delivers the exact positive triple (log links with V =
# mu^2 have Fisher weight 1, so the information is the unweighted Gram).  The
# zero rows are then recoloured to the original zero-rows targets, and the
# all-rows Gram is exact as the sum of the two blocks.  The zero-model score at
# the original zero fit is not pinned by these second moments.  It carries the
# within-cell higher-order mismatch, measured by the driver and bounded through
# clause (b) of the design-Gram proposition at the achieved deviation.  A dummy
# level absent from a block is absent from both trees by the count pinning.  Its
# columns are dropped from that block's pin.
.designgram.block <- function(m.X.orig, m.X.syn, v.cnt.idx, v.rows.o,
                              v.rows.s) {
    m.Xo <- m.X.orig[v.rows.o, , drop = FALSE]
    m.Xs <- m.X.syn[v.rows.s, , drop = FALSE]
    v.cat.idx <- setdiff(seq_len(ncol(m.X.orig)), v.cnt.idx)
    m.D.o <- m.Xo[, v.cat.idx, drop = FALSE]
    m.D.s <- m.Xs[, v.cat.idx, drop = FALSE]

    # drop dummy columns empty in this block (same set on both sides,
    # asserted, since the cell counts are pinned)
    v.keep <- colSums(abs(m.D.o)) > 0
    if(!identical(v.keep, colSums(abs(m.D.s)) > 0))
        stop("design-gram twopart: block dummy support differs")
    m.D.o <- m.D.o[, v.keep, drop = FALSE]
    m.D.s <- m.D.s[, v.keep, drop = FALSE]

    .designgram.recolour(m.D.s, m.D.o,
                         m.Xo[, v.cnt.idx, drop = FALSE],
                         m.Xs[, v.cnt.idx, drop = FALSE],
                         v.w = rep(1, nrow(m.D.s)))
}

.designgram.twopart <- function(d.syn, d.orig, y.nm, x.nms) {
    b.pos.o <- d.orig[[y.nm]] > 0
    b.pos.s <- d.syn[[y.nm]] > 0
    if(sum(b.pos.s) != sum(b.pos.o))
        stop("design-gram twopart: positive count not preserved")

    fml <- reformulate(x.nms)
    m.X.orig <- model.matrix(fml, d.orig)
    m.X.syn  <- model.matrix(fml, d.syn)
    if(!identical(colnames(m.X.orig), colnames(m.X.syn)))
        stop("design-gram twopart: model matrices disagree")
    v.types   <- mapply(id.cat, d.orig[x.nms])
    v.cnt.nm  <- names(v.types[v.types == "cnt"])
    if(length(v.cnt.nm) == 0) return(d.syn)
    v.cnt.idx <- which(colnames(m.X.orig) %in% v.cnt.nm)

    m.C.pos <- .designgram.block(m.X.orig, m.X.syn, v.cnt.idx,
                                 which(b.pos.o), which(b.pos.s))
    m.C.zer <- .designgram.block(m.X.orig, m.X.syn, v.cnt.idx,
                                 which(!b.pos.o), which(!b.pos.s))

    # Zero-model score refinement (the two-part full-synthesis corollary's score
    # condition).  The block pins fix the counts, the per-margin sums, and both
    # block Grams, which pins the linear part of the zero-model score at the
    # original fit b0 exactly (the z-weighted sums are the positive block's
    # pinned totals).  The remaining misfit is the nonlinear functional sum_i
    # x*_i p(x*_i'b0), large when the zero inflation is covariate-dependent
    # (Sc.4 and Sc.8 read 5.8% and 7.7% zero-block coefficient deviation without
    # this step).  Damped minimum-norm Gauss-Newton steps on the continuous
    # block, restricted per block to the span orthogonal to that block's dummy
    # columns (preserving every sum pin), drive the full score to tolerance, and
    # the block Grams are re-pinned exactly after each step.  The exact solution
    # exists: the row permutation of the original attains score zero and both
    # Grams simultaneously.  Loud stop on non-convergence.
    v.z <- as.numeric(b.pos.s)
    fml.z <- reformulate(x.nms, response = "z.designgram.ind")
    d.z.fit <- d.orig
    d.z.fit[["z.designgram.ind"]] <- as.numeric(b.pos.o)
    mo.z0 <- glm(fml.z, d.z.fit, family = binomial(),
                 control = glm.control(maxit = 200))
    v.b0 <- coef(mo.z0)[colnames(m.X.orig)]
    if(any(is.na(v.b0)))
        stop("design-gram twopart: zero model rank-deficient")
    v.b0.C <- v.b0[v.cnt.idx]
    v.cat.idx <- setdiff(seq_len(ncol(m.X.orig)), v.cnt.idx)

    # per-block projection onto the orthogonal complement of the block's
    # (non-empty) dummy columns, preserving the D'C sum pins
    .dg.orth <- function(m.E, v.rows, m.D.all) {
        m.D <- m.D.all[v.rows, , drop = FALSE]
        m.D <- m.D[, colSums(abs(m.D)) > 0, drop = FALSE]
        m.E[v.rows, ] <- m.E[v.rows, , drop = FALSE] -
            m.D %*% solve(crossprod(m.D),
                          crossprod(m.D, m.E[v.rows, , drop = FALSE]))
        m.E
    }
    m.D.all <- m.X.syn[, v.cat.idx, drop = FALSE]
    v.rows.p <- which(b.pos.s)
    v.rows.z <- which(!b.pos.s)

    m.C.cur <- matrix(NA_real_, nrow(m.X.syn), length(v.cnt.idx))
    m.C.cur[v.rows.p, ] <- m.C.pos
    m.C.cur[v.rows.z, ] <- m.C.zer

    # The step is constrained: it reduces the score while holding both
    # block Grams to first order (a per-iteration exact re-pin trades
    # first-order moves with the score step and oscillates).  The dual
    # space stacks the p score equations with the r(r+1)/2 symmetric
    # Gram-change equations per block.  The minimum-norm step solves the
    # small normal system K lambda = -[s; 0; 0] through the forward map
    # and its adjoint, and the exact re-pin runs only at convergence
    # checkpoints, where its perturbation is second order.
    s.r <- length(v.cnt.idx)
    m.vech <- which(lower.tri(diag(s.r), diag = TRUE), arr.ind = TRUE)
    s.nG <- nrow(m.vech)
    .sym.basis <- function(k) {
        m.B <- matrix(0, s.r, s.r)
        m.B[m.vech[k, 1], m.vech[k, 2]] <- 1
        m.B[m.vech[k, 2], m.vech[k, 1]] <- 1
        m.B
    }
    .gram.rows <- function(m.S) {
        mapply(function(k) `if`(m.vech[k, 1] == m.vech[k, 2],
                                m.S[m.vech[k, 1], m.vech[k, 2]],
                                2 * m.S[m.vech[k, 1], m.vech[k, 2]]),
               seq_len(s.nG))
    }
    .repin <- function(m.C) {
        d.tmp <- Reduce(function(d.acc, j) {
            d.acc[[v.cnt.nm[j]]] <- m.C[, j]
            d.acc
        }, seq_along(v.cnt.nm), init = d.syn)
        m.X.tmp <- model.matrix(fml, d.tmp)
        m.C[v.rows.p, ] <- .designgram.block(m.X.orig, m.X.tmp, v.cnt.idx,
                                             which(b.pos.o), v.rows.p)
        m.C[v.rows.z, ] <- .designgram.block(m.X.orig, m.X.tmp, v.cnt.idx,
                                             which(!b.pos.o), v.rows.z)
        m.C
    }

    # Stopping rule: the coefficient-relevant scale.  The refit's
    # zero-block coefficient moves by the Newton increment I(b0)^{-1} s
    # to first order, so that increment is what the tolerance governs.
    # The descent is monotone with a fast phase and a slow tail (the
    # residual rotates into the near-null space of the first-order
    # constraint-preserving map, where progress is second order), so
    # the loop takes the fast phase and accepts at the iteration cap
    # with the achieved increment, which the driver measures and the
    # companion table reports.  Convergence to 1e-6 or a stall below
    # 1e-4 exits early, and an achieved increment above the order-one
    # sanity bound stops loudly (a broken loop, not a slow one: over
    # the section-5 two-part seeds the achieved increments range to
    # 0.20 with no divergence).
    s.damp   <- 0.7
    s.p.dim  <- ncol(m.X.syn)
    s.q.dim  <- s.p.dim + 2 * s.nG
    s.maxit  <- 40L
    s.outer  <- 0L
    s.prev   <- Inf
    repeat {
        s.outer <- s.outer + 1L
        m.X.cur <- m.X.syn
        m.X.cur[, v.cnt.idx] <- m.C.cur
        v.eta <- as.numeric(m.X.cur %*% v.b0)
        v.p   <- plogis(v.eta)
        v.s   <- as.numeric(crossprod(m.X.cur, v.z - v.p))
        v.w0.chk <- v.p * (1 - v.p)
        m.I0  <- crossprod(m.X.cur * sqrt(v.w0.chk))
        s.db  <- max(abs(solve(m.I0, v.s)))
        b.stall <- abs(s.prev - max(abs(v.s))) < 1e-3 * max(abs(v.s))
        s.prev  <- max(abs(v.s))
        if(s.db < 1e-6 || (b.stall && s.db < 1e-4)) {
            m.C.cur <- .repin(m.C.cur)
            m.X.chk2 <- m.X.syn
            m.X.chk2[, v.cnt.idx] <- m.C.cur
            v.p2 <- plogis(as.numeric(m.X.chk2 %*% v.b0))
            v.s2 <- as.numeric(crossprod(m.X.chk2, v.z - v.p2))
            m.I2 <- crossprod(m.X.chk2 * sqrt(v.p2 * (1 - v.p2)))
            if(max(abs(solve(m.I2, v.s2))) < 1e-4) break
        }
        if(s.outer >= s.maxit) {
            if(s.db > 1)
                stop("design-gram twopart: score refinement failed to converge")
            m.C.cur <- .repin(m.C.cur)
            break
        }

        v.w0 <- v.p * (1 - v.p)
        m.C.p <- m.C.cur[v.rows.p, , drop = FALSE]
        m.C.z <- m.C.cur[v.rows.z, , drop = FALSE]
        # forward map: score change and the two blocks' Gram changes
        L.fwd <- function(m.E) {
            v.sc <- -as.numeric(crossprod(m.X.cur, v.w0 *
                                          as.numeric(m.E %*% v.b0.C)))
            v.sc[v.cnt.idx] <- v.sc[v.cnt.idx] +
                as.numeric(crossprod(m.E, v.z - v.p))
            m.Ep <- m.E[v.rows.p, , drop = FALSE]
            m.Ez <- m.E[v.rows.z, , drop = FALSE]
            m.Sp <- crossprod(m.Ep, m.C.p)
            m.Sz <- crossprod(m.Ez, m.C.z)
            c(v.sc, .gram.rows(m.Sp + t(m.Sp)), .gram.rows(m.Sz + t(m.Sz)))
        }
        # adjoint: score dual a plus one symmetric dual per block
        L.adj <- function(v.dual) {
            v.a  <- v.dual[seq_len(s.p.dim)]
            v.gp <- v.dual[s.p.dim + seq_len(s.nG)]
            v.gz <- v.dual[s.p.dim + s.nG + seq_len(s.nG)]
            m.Lp <- Reduce(`+`, Map(function(k) v.gp[k] * .sym.basis(k),
                                    seq_len(s.nG)), matrix(0, s.r, s.r))
            m.Lz <- Reduce(`+`, Map(function(k) v.gz[k] * .sym.basis(k),
                                    seq_len(s.nG)), matrix(0, s.r, s.r))
            m.E <- -(v.w0 * as.numeric(m.X.cur %*% v.a)) %o% v.b0.C
            m.E <- m.E + (v.z - v.p) %o% v.a[v.cnt.idx]
            m.E[v.rows.p, ] <- m.E[v.rows.p, , drop = FALSE] +
                2 * m.C.p %*% m.Lp
            m.E[v.rows.z, ] <- m.E[v.rows.z, , drop = FALSE] +
                2 * m.C.z %*% m.Lz
            m.E <- .dg.orth(m.E, v.rows.p, m.D.all)
            .dg.orth(m.E, v.rows.z, m.D.all)
        }
        m.K <- mapply(function(k) {
            L.fwd(L.adj(`[<-`(numeric(s.q.dim), k, 1)))
        }, seq_len(s.q.dim))
        v.rhs <- c(-v.s, numeric(2 * s.nG))
        v.lambda <- tryCatch(
            solve(m.K + diag(1e-10 * max(abs(m.K)), s.q.dim), v.rhs),
            error = function(e)
                stop("design-gram twopart: score system singular"))
        m.C.cur <- m.C.cur + s.damp * L.adj(v.lambda)
    }
    m.C.pos <- m.C.cur[v.rows.p, , drop = FALSE]
    m.C.zer <- m.C.cur[v.rows.z, , drop = FALSE]

    d.out <- Reduce(function(d.acc, j) {
        v <- d.acc[[v.cnt.nm[j]]]
        v[b.pos.s]  <- m.C.pos[, j]
        v[!b.pos.s] <- m.C.zer[, j]
        d.acc[[v.cnt.nm[j]]] <- v
        d.acc
    }, seq_along(v.cnt.nm), init = d.syn)

    # postconditions: both block Grams, hence the all-rows Gram, exact
    m.X.chk <- model.matrix(fml, d.out)
    s.scl <- max(abs(crossprod(m.X.orig)))
    if(max(abs(crossprod(m.X.chk[b.pos.s, , drop = FALSE]) -
               crossprod(m.X.orig[b.pos.o, , drop = FALSE]))) > 1e-8 * s.scl)
        stop("design-gram twopart: positive-block Gram violated")
    if(max(abs(crossprod(m.X.chk) - crossprod(m.X.orig))) > 1e-8 * s.scl)
        stop("design-gram twopart: all-rows Gram violated")
    d.out
}

# +------------------------------------------------------------------------+
# | Direction-agnostic generator (Scenario 1: unknown Y)                  |
# +------------------------------------------------------------------------+
# When any variable could be the response, we need the synthetic data to
# give exact inference regardless of which variable is designated Y.
#
# For OLS (Gaussian/identity): project.crossproduct pins D'C* and C*'C*,
# and the base draw must pin the third block D'D as well, so the label
# draw is the strata.exact permutation (condition (i) of the design-Gram
# projection proposition).  With all three blocks pinned, Z*'Z* = Z'Z
# exactly and every OLS regression among the variables is exact
# simultaneously, coefficients and standard errors.  On a free bootstrap
# draw the cell counts move (multinomial noise), the enforced original
# column sums then distort the per-stratum means by the count ratio, and
# dummy-involved coefficients absorb level-scale errors, so exactness
# Requires the permutation draw, and this function always uses it.
#
# For non-OLS GLMs: sequential per-variable projection.  Each continuous
# variable is projected in turn, treating it as Y with the rest as X.
# This is approximate, since projecting variable j re-perturbs the score
# for variable k, and no exactness claim is made for this path.
gen.syn.projected.symmetric <- function(d, family = "gaussian",
                                        link = "identity",
                                        supports = NULL,
                                        order = 4, seed = 6489) {
    d.syn <- gen.syn.family(d, order = order, supports = supports,
                            shrink = TRUE, seed = seed,
                            strata.exact = TRUE)

    if(family == "gaussian" && link == "identity") {
        # exact: cross-product matching gives all OLS regressions at once
        return(project.crossproduct(d.syn, d))
    }

    # non-OLS: project each continuous variable sequentially
    v.types <- mapply(id.cat, d)
    v.cnt.nm <- names(v.types[v.types == "cnt"])

    Reduce(function(d.acc, y.nm) {
        v.x.nms <- setdiff(names(d), y.nm)
        v.y.proj <- project.glm(d.acc, d, y.nm, v.x.nms, family, link)
        d.acc[[y.nm]] <- v.y.proj
        d.acc
    }, v.cnt.nm, init = d.syn)
}
# +------------------------------------------------------------------------+
# | Test harness                                                          |
# +------------------------------------------------------------------------+
# Runs only in standalone mode (not when sourced by the Sobol sensitivity
# analysis in the main simulation pipeline).  Verifies that projection
# achieves 0.0% deviation for beta_hat and SE across all GLM families.
if(!exists(".sobol.main")) {
load("results/d.sim.scenarios.RData")

# support declarations per scenario (controls copula transforms):
#   real     = (-inf, +inf), no transform
#   positive = (0, +inf),    log transform before copula
#   zeroinfl = {0} U (0,+inf), zero indicator + log(positive part)
#   count    = {0,1,2,...},  treated as zeroinfl for copula
#   bounded  = (0,1),       logit transform
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

# helper: fit the analysis model on a dataset (original or synthetic).
# dispatches to lm / poisson glm / quasi-gamma glm as appropriate.
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

# significance stars for display
sig <- function(mo) {
    v.p <- summary(mo)[["coefficients"]][, 4]
    ifelse(v.p < 0.001, "***",
    ifelse(v.p < 0.01,  "** ",
    ifelse(v.p < 0.05,  "*  ", "   ")))
}

# percentage deviation of synthetic from original (per coefficient)
pct.dev <- function(v.syn, v.orig)
    round(100 * (v.syn - v.orig) / abs(v.orig), 1)

# maximum absolute percentage deviation across all coefficients
max.abs.pct <- function(v.syn, v.orig)
    max(abs(pct.dev(v.syn, v.orig)), na.rm = TRUE)

# +------------------------------------------------------------------------+
# | Scenario 2: known model + known Y (projection)                       |
# +------------------------------------------------------------------------+
# The analyst specifies both the GLM family/link and which variable is Y.
# This is the main use case.  For each scenario, compare:
#   copula  = uncalibrated copula baseline (distributional fidelity only)
#   proj    = copula + projection (should give 0.0% deviation)
cat("=== SCENARIO 2: known model + known Y (projection) ===\n\n")

# All ten scenarios, so every per-coefficient block the paper reports has this
# block as its emitting script.
v.test <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
Map(function(i) {
    s <- l.scenarios[[i]]
    d.orig <- s[["data"]]
    sup <- l.supports[[i]]
    mo.orig <- fit.model(d.orig, s[["family"]], s[["link"]])
    if(is.null(mo.orig)) return(NULL)
    v.orig <- coef(mo.orig)

    # projected generator
    d.proj <- gen.syn.projected(d.orig, y.nm = "response",
                                family = s[["family"]], link = s[["link"]],
                                supports = sup)
    mo.proj <- fit.model(d.proj, s[["family"]], s[["link"]])

    # copula baseline (uncalibrated)
    d.cop <- gen.syn.family(d.orig, order = 4, supports = sup)
    mo.cop <- fit.model(d.cop, s[["family"]], s[["link"]])

    v.proj <- `if`(!is.null(mo.proj), coef(mo.proj), rep(NA, length(v.orig)))
    v.cop  <- `if`(!is.null(mo.cop), coef(mo.cop), rep(NA, length(v.orig)))

    cat(sprintf("--- Scenario %d: %s(%s)  max|%%|: copula=%.1f%%  proj=%.1f%% ---\n",
                i, s[["family"]], s[["link"]],
                max.abs.pct(v.cop, v.orig),
                max.abs.pct(v.proj, v.orig)))
    d.out <- data.frame(
        orig  = sprintf("%8.4f%s", v.orig, sig(mo.orig)),
        copul = sprintf("%8.4f (%+.1f%%)", v.cop, pct.dev(v.cop, v.orig)),
        proj  = sprintf("%8.4f (%+.1f%%)", v.proj, pct.dev(v.proj, v.orig)))
    print(d.out)

    # compare SEs
    v.se.orig <- summary(mo.orig)[["coefficients"]][, 2]
    v.se.proj <- `if`(!is.null(mo.proj),
                      summary(mo.proj)[["coefficients"]][, 2],
                      rep(NA, length(v.se.orig)))
    v.se.cop  <- `if`(!is.null(mo.cop),
                      summary(mo.cop)[["coefficients"]][, 2],
                      rep(NA, length(v.se.orig)))
    cat(sprintf("  SE max|%%|: copula=%.1f%%  proj=%.1f%%\n\n",
                max.abs.pct(v.se.cop, v.se.orig),
                max.abs.pct(v.se.proj, v.se.orig)))
}, v.test)

# +------------------------------------------------------------------------+
# | Canonical vs non-canonical link tests                                 |
# +------------------------------------------------------------------------+
# Verify that projection works for both canonical links (where the score
# is linear in y, so one-step projection suffices) and non-canonical links
# (where Newton iteration is needed).
cat("\n=== CANONICAL LINK TESTS ===\n\n")

# Gamma(inverse): canonical link for Gamma family.  V(mu) g'(mu) = 1
# with V(mu) = mu^2, g(mu) = 1/mu, g'(mu) = -1/mu^2.
cat("--- Gamma(inverse) on Scenario 3 data ---\n")
d.sc3 <- l.scenarios[[3]][["data"]]
sup.sc3 <- l.supports[[3]]
d.proj.gi <- gen.syn.projected(d.sc3, y.nm = "response",
                               family = "gamma", link = "inverse",
                               supports = sup.sc3)
d.cop.sc3 <- gen.syn.family(d.sc3, order = 4, supports = sup.sc3)

fit.gi <- function(d) tryCatch(
    glm(response ~ age + female + race + income, d,
        family = Gamma(link = "inverse"),
        start = c(0.05, rep(0, 5)),
        control = glm.control(maxit = 100)),
    error = function(e) NULL)
mo.gi.o <- fit.gi(d.sc3)
mo.gi.p <- fit.gi(d.proj.gi)
mo.gi.c <- fit.gi(d.cop.sc3)
if(!is.null(mo.gi.o) && !is.null(mo.gi.p) && !is.null(mo.gi.c)) {
    v.o <- coef(mo.gi.o)
    v.p <- coef(mo.gi.p)
    v.c <- coef(mo.gi.c)
    cat(sprintf("  coef max|%%|: copula=%.1f%%  proj=%.1f%%\n",
                max.abs.pct(v.c, v.o), max.abs.pct(v.p, v.o)))
    v.se.o <- summary(mo.gi.o)[["coefficients"]][, 2]
    v.se.p <- summary(mo.gi.p)[["coefficients"]][, 2]
    cat(sprintf("  SE   max|%%|: proj=%.1f%%\n", max.abs.pct(v.se.p, v.se.o)))
    d.out <- data.frame(
        orig  = sprintf("%10.6f%s", v.o, sig(mo.gi.o)),
        copul = sprintf("%10.6f (%+.1f%%)", v.c, pct.dev(v.c, v.o)),
        proj  = sprintf("%10.6f (%+.1f%%)", v.p, pct.dev(v.p, v.o)))
    print(d.out)
} else cat("  fit failed\n")

# Poisson(log): canonical link.  Integer-valued, so projection uses
# log-scale Newton + probabilistic rounding (~2.3% residual error).
cat("\n--- Poisson(log) [canonical] on Scenario 6 ---\n")
d.sc6 <- l.scenarios[[6]][["data"]]
sup.sc6 <- l.supports[[6]]
d.proj.p <- gen.syn.projected(d.sc6, y.nm = "response",
                              family = "poisson", link = "log",
                              supports = sup.sc6)
d.cop.sc6 <- gen.syn.family(d.sc6, order = 4, supports = sup.sc6)

fit.pois <- function(d) tryCatch(
    glm(response ~ age + female + race + income, d, family = poisson()),
    error = function(e) NULL)
mo.p.o <- fit.pois(d.sc6)
mo.p.p <- fit.pois(d.proj.p)
mo.p.c <- fit.pois(d.cop.sc6)
if(!is.null(mo.p.o) && !is.null(mo.p.p) && !is.null(mo.p.c)) {
    v.o <- coef(mo.p.o)
    v.p <- coef(mo.p.p)
    v.c <- coef(mo.p.c)
    cat(sprintf("  coef max|%%|: copula=%.1f%%  proj=%.1f%%\n",
                max.abs.pct(v.c, v.o), max.abs.pct(v.p, v.o)))
    v.se.o <- summary(mo.p.o)[["coefficients"]][, 2]
    v.se.p <- summary(mo.p.p)[["coefficients"]][, 2]
    cat(sprintf("  SE   max|%%|: proj=%.1f%%\n", max.abs.pct(v.se.p, v.se.o)))
    d.out <- data.frame(
        orig  = sprintf("%8.4f%s", v.o, sig(mo.p.o)),
        copul = sprintf("%8.4f (%+.1f%%)", v.c, pct.dev(v.c, v.o)),
        proj  = sprintf("%8.4f (%+.1f%%)", v.p, pct.dev(v.p, v.o)))
    print(d.out)
} else cat("  fit failed\n")

# Gamma(log): non-canonical link.  Newton iteration needed because
# the score X' diag(w) (y - mu) = 0 is nonlinear in y when g != g_canonical.
cat("\n--- Gamma(log) [non-canonical] on Scenario 3 ---\n")
d.proj.gl <- gen.syn.projected(d.sc3, y.nm = "response",
                               family = "gamma", link = "log",
                               supports = sup.sc3)
fit.gl <- function(d) tryCatch(
    glm(response ~ age + female + race + income, d,
        family = Gamma(link = "log")),
    error = function(e) NULL)
mo.gl.o <- fit.gl(d.sc3)
mo.gl.p <- fit.gl(d.proj.gl)
mo.gl.c <- fit.gl(d.cop.sc3)
if(!is.null(mo.gl.o) && !is.null(mo.gl.p) && !is.null(mo.gl.c)) {
    v.o <- coef(mo.gl.o)
    v.p <- coef(mo.gl.p)
    v.c <- coef(mo.gl.c)
    cat(sprintf("  coef max|%%|: copula=%.1f%%  proj=%.1f%%\n",
                max.abs.pct(v.c, v.o), max.abs.pct(v.p, v.o)))
    v.se.o <- summary(mo.gl.o)[["coefficients"]][, 2]
    v.se.p <- summary(mo.gl.p)[["coefficients"]][, 2]
    cat(sprintf("  SE   max|%%|: proj=%.1f%%\n\n", max.abs.pct(v.se.p, v.se.o)))
} else cat("  fit failed\n\n")

# +------------------------------------------------------------------------+
# | Scenario 1: known model + unknown Y (direction-agnostic)             |
# +------------------------------------------------------------------------+
# The analyst specifies the GLM family but not which variable is Y.
# Every continuous variable is tried as the response.  For OLS, the
# cross-product projection on the strata.exact permutation draw pins the
# full Gram Z*'Z* = Z'Z, so every direction is exact simultaneously,
# coefficients and standard errors (the copula column keeps its free
# draw and is the comparison).
cat("=== SCENARIO 1: known model + unknown Y (direction-agnostic) ===\n\n")

# helper: fit OLS with y.nm as response, rest as predictors
fit.ols <- function(d, y.nm) {
    v.pred <- setdiff(intersect(
        c("response", "age", "female", "race", "income"), colnames(d)),
        y.nm)
    fml <- reformulate(v.pred, response = y.nm)
    tryCatch(lm(fml, d), error = function(e) NULL)
}

# Gaussian/identity: test all three continuous directions
d.orig <- l.scenarios[[1]][["data"]]
sup <- l.supports[[1]]
d.proj.sym <- gen.syn.projected.symmetric(d.orig, family = "gaussian",
                                          link = "identity", supports = sup)
d.cop <- gen.syn.family(d.orig, order = 4, supports = sup)

Map(function(y.nm) {
    mo.o <- fit.ols(d.orig, y.nm)
    mo.p <- fit.ols(d.proj.sym, y.nm)
    mo.c <- fit.ols(d.cop, y.nm)
    if(is.null(mo.o) || is.null(mo.p) || is.null(mo.c)) return(NULL)
    v.o <- coef(mo.o)
    v.p <- coef(mo.p)
    v.c <- coef(mo.c)
    cat(sprintf("--- Y = %s  max|%%|: copula=%.1f%%  proj=%.1f%% ---\n",
                y.nm, max.abs.pct(v.c, v.o), max.abs.pct(v.p, v.o)))
    d.out <- data.frame(
        orig  = sprintf("%8.4f", v.o),
        copul = sprintf("%8.4f (%+.1f%%)", v.c, pct.dev(v.c, v.o)),
        proj  = sprintf("%8.4f (%+.1f%%)", v.p, pct.dev(v.p, v.o)))
    print(d.out)
    cat("\n")
}, c("response", "age", "income"))

# +------------------------------------------------------------------------+
# | ZI test: Scenario 4 (zero-inflated quasi-gamma)                      |
# +------------------------------------------------------------------------+
# Verifies ZI handling: copula bootstrap preserves zero proportion,
# projection on the positive subset gives exact conditional inference.
cat("=== ZI TEST: Scenario 4 (zero-inflated quasi-gamma) ===\n\n")
d.orig <- l.scenarios[[4]][["data"]]
sup <- l.supports[[4]]
d.proj <- gen.syn.projected(d.orig, y.nm = "response",
                            family = "quasi-gamma", link = "log",
                            supports = sup)
mo.orig <- fit.model(d.orig, "quasi-gamma", "log")
mo.proj <- fit.model(d.proj, "quasi-gamma", "log")
d.cop <- gen.syn.family(d.orig, order = 4, supports = sup)
mo.cop <- fit.model(d.cop, "quasi-gamma", "log")

v.orig <- coef(mo.orig)
v.proj <- coef(mo.proj)
v.cop  <- coef(mo.cop)

cat(sprintf("Zero proportions -- orig: %.1f%%  proj: %.1f%%  copula: %.1f%%\n",
            100 * mean(d.orig[["response"]] == 0),
            100 * mean(d.proj[["response"]] == 0),
            100 * mean(d.cop[["response"]] == 0)))

cat(sprintf("Coef max|%%| -- copula: %.1f%%  proj: %.1f%%\n",
            max.abs.pct(v.cop, v.orig), max.abs.pct(v.proj, v.orig)))

v.se.orig <- summary(mo.orig)[["coefficients"]][, 2]
v.se.proj <- summary(mo.proj)[["coefficients"]][, 2]
v.se.cop  <- summary(mo.cop)[["coefficients"]][, 2]
cat(sprintf("SE max|%%|   -- copula: %.1f%%  proj: %.1f%%\n\n",
            max.abs.pct(v.se.cop, v.se.orig),
            max.abs.pct(v.se.proj, v.se.orig)))

d.out <- data.frame(
    orig = sprintf("%8.4f%s", v.orig, sig(mo.orig)),
    copul = sprintf("%8.4f (%+.1f%%)", v.cop, pct.dev(v.cop, v.orig)),
    proj  = sprintf("%8.4f (%+.1f%%)", v.proj, pct.dev(v.proj, v.orig)))
print(d.out)
cat("\n")
} # end if(!exists(".sobol.main"))
# +------------------------------------------------------------------------+
