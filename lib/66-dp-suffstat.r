#!/usr/bin/env r
# +------------------------------------------------------------------------+
# |  66-dp-suffstat.r -- Route 3: the whitened sufficient-statistic release |
# +------------------------------------------------------------------------+
#
# What this is.  A joint Gaussian release of the OLS
# sufficient statistics on the whitened public scale,
#
#     release   vech(X_w' X_w),   X_w' y_s,   y_s' y_s
#     solve     beta_w = (S_dp + lambda I)^{-1} b_dp        (post-processing)
#
# with w = W x, W = Sigma_pop^{-1/2} the public whitening map, and y_s the
# publicly standardised, publicly clipped response.  X is private here: the
# three sensitivities are functionals of (W, the public box, the public
# clip), not of the realised design, which is what Route 1 cannot offer
# (dp.public.B.x.design's warning).  The survey (65-dp-suffstat-survey.r)
# measured maxpct 29.8% at eps = 10 on Scenario 1 under a more conservative
# vech constant, against subagg's 66.7% at kappa = 3, and this file is the
# mechanism that survey approved.
#
# The guarantee.  Three Gaussian channels in one mu quadrature,
#
#     mu_total^2 = mu_vech^2 + mu_xty^2 + mu_yy^2,
#
# each calibrated at its L2 sensitivity from dp.public.suffstat.sens():
# sqrt(2) B_x^2 for the vech (a proven, strict upper bound on the whitened
# box), 2 B_x B_y for X'y, and B_y^2 for y'y (attained).  Everything after
# the release --
# the lambda, the solve, the PD repair, the back-map, the synthetic draw --
# reads only released and public objects, so it is post-processing and
# spends nothing.
#
# The lambda is not part of the guarantee.  Any lambda that reads only
# public and released objects preserves mu-GDP.  Its value moves utility
# only.  Two rules are provided:
#
#   "public-curve"       lambda* minimises the population error curve, a
#                        deterministic functional of the published DGP
#                        computed at a fixed public seed (the same premise
#                        as every public box here).  The default, and it spends no budget.
#   "released-eigenmin"  the AdaSSP-style rule for deployments without a
#                        published DGP: spend a budget share releasing
#                        lambda_min(S) (Weyl gives replace-one sensitivity
#                        B_x^2) and set lambda from the released value.
#                        Opt-in, costs budget, heuristic in its constant
#                        (the guarantee does not depend on it).
#
# What the analyst refits.  On the raw scale the solve returns beta_DP with
# a known error law around the ridge target: conditional on S_dp,
#
#     Cov(beta_w | S_dp) ~= M [ sigma_xty^2 I
#                               + sigma_vech^2 (|beta|^2 I + beta beta') ] M,
#     M = (S_dp + lambda I)^{-1},
#
# using Var((E beta)_i) = sigma_vech^2 |beta|^2 and Cov((E beta)_i, (E beta)_k)
# = sigma_vech^2 beta_i beta_k for symmetric E with iid upper-triangle noise.
# The plug-in at beta_w_dp is returned as Cov.beta.dp (raw scale), with the
# Exact de-shrink bias.fold = beta.infer - beta.dp beside it.
#
# Calibrated inference recentres.  Without it the driver reads coverage as low
# as 66% at the strong-shrinkage cells
# (eps <= 2, where the public curve picks lambda = n) for the ridge
# interval beta.dp +/- 1.96 sqrt(SE_syn^2 + diag Cov + bias^2), for two
# reasons.  The first-order fold lambda M beta_w_dp reads half the true
# shrinkage at lambda = n, and any fold proportional to the released
# beta_w shrinks exactly when the noise has pushed a coefficient toward
# zero, which is when the error is largest.  The calibrated interval is
# instead centred at the lambda-0 solve beta.infer = S_dp^{-1} b_dp (an
# exact post-processing identity, beta_w0 = beta_w + lambda S_dp^{-1}
# beta_w), with the same covariance plug-in evaluated at the de-shrunk
# coefficient, returned as beta.infer and Cov.beta.infer.  beta.dp stays
# the point estimate (the ridge MSE choice), and the recentred interval
# reports the uncertainty the release actually carries, so at a
# near-singular S_dp it is very large.
#
# Scope.  Gaussian identity link, scenarios with a public B_y (1, 2, 9, 10,
# and the PHI application).  The unbounded-response scenarios have no public
# clip and hence no y channel: subagg remains the only route there, and the
# the re-scoped two-stage work targets them.
#
# Dependencies: 28-dp-projection.r (gdp.mu.from.eps.delta, synth.response),
#               dp-public-bounds.r (the three public-constant constructors).
# This file defines a mechanism and nothing else.  The driver is
# 66b-mc-dp-suffstat.r.
# +------------------------------------------------------------------------+

if(!exists("gdp.mu.from.eps.delta")) {
    .dp.suppress.demo <- TRUE; source("lib/28-dp-projection.r")
}

# +------------------------------------------------------------------------+
# |  The released-object algebra, shared by the mechanism and its verify    |
# |  tasks                                                                  |
# +------------------------------------------------------------------------+
# Statistics of the whitened, publicly clipped data.  Exposed separately so
# verify-dp-soundness.R can drive adversarial neighbouring rows through the
# Exact map the mechanism releases and measure the move of each statistic.
dp.suffstat.stats <- function(d, y.nm, l.white, l.ystd) {
    v.x.nms <- setdiff(names(d), y.nm)
    m.X <- model.matrix(reformulate(v.x.nms), dp.public.clip.x(d))
    m.Xw <- m.X %*% l.white[["W"]]                      # W symmetric
    v.yc <- pmin(pmax(d[[y.nm]], l.ystd[["y.clip"]][[1L]]),
                 l.ystd[["y.clip"]][[2L]])
    v.ys <- (v.yc - l.ystd[["centre"]]) / l.ystd[["scale"]]
    list(S = crossprod(m.Xw), b = as.numeric(crossprod(m.Xw, v.ys)),
         q = sum(v.ys^2), m.X = m.X, n = nrow(m.X), p = ncol(m.X))
}

# vech extraction and symmetric reassembly, the release's own layout.
dp.suffstat.vech <- function(m.S) m.S[upper.tri(m.S, diag = TRUE)]
dp.suffstat.unvech <- function(v, s.p) {
    m <- matrix(0, s.p, s.p)
    m[upper.tri(m, diag = TRUE)] <- v
    m + t(m) - diag(diag(m))
}

# The public-curve lambda: the argmin of the population error curve over a
# public grid, computed by a fixed-public-seed MC from published constants
# only (beta_pop through the same whitening algebra, S_pop = n I exactly by
# the whitening's construction).  Deterministic in its arguments, reads no
# sample statistic, so it is public in the same sense as every box from
# dp-public-bounds.r.
dp.suffstat.lambda.public <- function(sc, s.n, eps, delta, mu.split,
                                      l.white, l.ystd, l.sens,
                                      v.grid = c(0, 1e-4, 3e-4, 1e-3, 3e-3,
                                                 0.01, 0.03, 0.1, 0.3, 1),
                                      B.pop = 400L, pub.seed = 20260713L) {
    l.k <- readRDS("results/dp-design-constants.rds")[[paste0("sc", sc)]]
    v.beta.pop <- l.k[["beta"]]
    # beta_pop on the (whitened x, standardised y) scale: undo the y scale,
    # move the y centre into the intercept, then beta_w = W^{-1} beta = the
    # solve against W (W symmetric).
    v.b.std <- v.beta.pop / l.ystd[["scale"]]
    v.b.std[[1L]] <- (v.beta.pop[[1L]] - l.ystd[["centre"]]) / l.ystd[["scale"]]
    v.beta.w <- solve(l.white[["W"]], v.b.std)
    s.p <- length(v.beta.w)

    v.mu <- sqrt(mu.split) * gdp.mu.from.eps.delta(eps, delta)
    s.sig.S <- l.sens[["Delta.vech"]] / v.mu[["vech"]]
    s.sig.b <- l.sens[["Delta.xty"]]  / v.mu[["xty"]]

    m.S.pop <- diag(s.p) * s.n                          # whitened population design
    v.b.pop <- as.numeric(m.S.pop %*% v.beta.w)
    set.seed(pub.seed)
    v.obj <- numeric(length(v.grid))
    for(i in seq_len(B.pop)) {
        m.E <- dp.suffstat.unvech(rnorm(s.p * (s.p + 1) / 2, sd = s.sig.S), s.p)
        v.bn <- v.b.pop + rnorm(s.p, sd = s.sig.b)
        l.e <- eigen(m.S.pop + m.E, symmetric = TRUE)
        v.proj <- as.numeric(crossprod(l.e[["vectors"]], v.bn))
        for(k in seq_along(v.grid)) {
            v.d <- l.e[["values"]] + v.grid[[k]] * s.n
            if(min(v.d) <= 1e-8) { v.obj[[k]] <- v.obj[[k]] + 1e6; next }
            v.bw <- as.numeric(l.e[["vectors"]] %*% (v.proj / v.d))
            v.obj[[k]] <- v.obj[[k]] +
                sum(((v.bw - v.beta.w) / pmax(abs(v.beta.w), 1e-10))^2)
        }
    }
    v.grid[[which.min(v.obj)]] * s.n
}

# +------------------------------------------------------------------------+
# |  gen.syn.dp.suffstats: the mechanism                                    |
# +------------------------------------------------------------------------+
# Arguments:
#   d             original data frame (private).
#   y.nm          response column.
#   sc            scenario id, for the public lambda curve only.
#   eps, delta    the privacy budget.  Converted once to mu_total.
#   mu.split      named fractions of mu^2 for (vech, xty, yy), summing to 1.
#   whitening     dp.public.whitening(sc).   Required, no default: a generator
#   y.std         dp.public.y.std(sc).       must refuse to invent its own
#   sens          dp.public.suffstat.sens(). bounds, by design.
#   lambda.method "public-curve" (default) or "released-eigenmin" (opt-in,
#                 diverts eig.split of mu^2 from the vech channel to release
#                 lambda_min(S) at the Weyl sensitivity B_x^2).
#   dispersion    NULL (default) privatises the variance through the y'y
#                 channel.  A dp.dispersion.public object instead declares
#                 sigma^2 public (the Route 2 free mode, 28-dp-projection.r):
#                 the y'y channel is dropped, its mu^2 share renormalised
#                 over the two remaining channels, which is the like-for-like
#                 footing against subagg's recorded numbers.
#   x.public      TRUE additionally emits d.syn on the original design via
#                 synth.response (X public, the Route 1 premise).  FALSE (the
#                 default, X private) emits the statistics-only release for
#                 the V4-style pipeline.
#   seed          RNG seed for the noise draws.
gen.syn.dp.suffstats <- function(d, y.nm = "response", sc,
                                 eps = 1, delta = 1e-6,
                                 mu.split = c(vech = 0.45, xty = 0.45, yy = 0.10),
                                 whitening = NULL, y.std = NULL, sens = NULL,
                                 lambda.method = c("public-curve",
                                                   "released-eigenmin"),
                                 eig.split = 0.10,
                                 dispersion = NULL,
                                 x.public = FALSE, seed = 6489) {
    lambda.method <- match.arg(lambda.method)
    if(is.null(whitening) || is.null(y.std) || is.null(sens))
        stop("whitening, y.std, and sens are REQUIRED, from dp.public.whitening, ",
             "dp.public.y.std, and dp.public.suffstat.sens. A generator does ",
             "not invent its own bounds.")
    if(abs(sum(mu.split) - 1) > 1e-12 ||
       !identical(sort(names(mu.split)), c("vech", "xty", "yy")))
        stop("mu.split must be named fractions (vech, xty, yy) summing to 1")
    b.disp.public <- inherits(dispersion, "dp.dispersion")
    if(b.disp.public && dispersion[["mode"]] != "public")
        stop("dispersion must be NULL (private y'y channel) or a ",
             "dp.dispersion.public object; the private-box mode belongs to ",
             "the subagg channel, not here")
    if(b.disp.public) {
        mu.split[["yy"]] <- 0
        mu.split <- mu.split / sum(mu.split)
    }

    l.st <- dp.suffstat.stats(d, y.nm, whitening, y.std)
    s.n <- l.st[["n"]]; s.p <- l.st[["p"]]

    # The mu accounting.  The released-eigenmin rule takes its share out of
    # the vech channel (the matrix side), leaving xty and yy untouched.
    s.mu.total <- gdp.mu.from.eps.delta(eps, delta)
    v.f <- mu.split
    s.f.eig <- 0
    if(lambda.method == "released-eigenmin") {
        s.f.eig <- eig.split * v.f[["vech"]]
        v.f[["vech"]] <- v.f[["vech"]] - s.f.eig
    }
    s.sig.S  <- sens[["Delta.vech"]] / (sqrt(v.f[["vech"]]) * s.mu.total)
    s.sig.b  <- sens[["Delta.xty"]]  / (sqrt(v.f[["xty"]])  * s.mu.total)
    s.sig.q  <- `if`(b.disp.public, NA_real_,
                     sens[["Delta.yy"]] / (sqrt(v.f[["yy"]]) * s.mu.total))

    # The release: the Gaussian channels, one quadrature.
    set.seed(seed + 1L)
    m.S.dp <- l.st[["S"]] +
        dp.suffstat.unvech(rnorm(s.p * (s.p + 1) / 2, sd = s.sig.S), s.p)
    set.seed(seed + 2L)
    v.b.dp <- l.st[["b"]] + rnorm(s.p, sd = s.sig.b)
    set.seed(seed + 3L)
    s.q.dp <- `if`(b.disp.public, NA_real_,
                   l.st[["q"]] + rnorm(1L, sd = s.sig.q))

    # The lambda.  Both rules read only public and released objects.
    s.lambda <- if(lambda.method == "public-curve") {
        dp.suffstat.lambda.public(sc, s.n, eps, delta, v.f,
                                  whitening, y.std, sens)
    } else {
        s.sig.eig <- whitening[["B.x"]]^2 / (sqrt(s.f.eig) * s.mu.total)
        set.seed(seed + 4L)
        s.eigmin.dp <- min(eigen(m.S.dp, symmetric = TRUE,
                                 only.values = TRUE)[["values"]]) +
            rnorm(1L, sd = s.sig.eig)
        # Heuristic constant (utility only, not the guarantee): raise the
        # solve's floor to the scale of the matrix noise if the released
        # eigenvalue sits below it.
        max(0, 2 * s.sig.S * sqrt(s.p) - max(s.eigmin.dp, 0))
    }

    # PD repair, post-processing of the released matrix: the solve needs
    # S_dp + lambda I positive definite, so the effective lambda is raised
    # by any released eigenvalue deficit plus a small public floor.
    s.eig.min <- min(eigen(m.S.dp + s.lambda * diag(s.p), symmetric = TRUE,
                           only.values = TRUE)[["values"]])
    s.lambda.eff <- s.lambda + max(0, 1e-6 * s.n - s.eig.min)

    m.M <- solve(m.S.dp + s.lambda.eff * diag(s.p))
    v.beta.w <- as.numeric(m.M %*% v.b.dp)

    # Back to the analyst's raw scale: beta_xstd = W beta_w (W symmetric),
    # then undo the y standardisation, the centre landing on the intercept.
    v.beta.dp <- as.numeric(whitening[["W"]] %*% v.beta.w) * y.std[["scale"]]
    v.beta.dp[[1L]] <- v.beta.dp[[1L]] + y.std[["centre"]]
    names(v.beta.dp) <- colnames(l.st[["m.X"]])

    # sigma^2: from the released statistics alone (floored like Route 1), or
    # the declared public value, whose provenance the object itself carries.
    s.sigma2.dp <- if(b.disp.public) {
        dispersion[["sigma2"]]
    } else {
        s.rss.w <- s.q.dp - 2 * sum(v.beta.w * v.b.dp) +
            as.numeric(t(v.beta.w) %*% m.S.dp %*% v.beta.w)
        max(s.rss.w * y.std[["scale"]]^2 / (s.n - s.p), 1e-6)
    }

    # The analyst-side error law (plug-in), on the raw scale.  Var(E beta) =
    # sigma_S^2 (|beta|^2 I + beta beta') for the symmetric vech noise, so
    # conditional on S_dp the coefficient noise covariance and everything
    # below are functions of released and public objects only.
    m.cov.w <- m.M %*% (s.sig.b^2 * diag(s.p) +
                        s.sig.S^2 * (sum(v.beta.w^2) * diag(s.p) +
                                     tcrossprod(v.beta.w))) %*% m.M
    m.cov.raw <- (y.std[["scale"]]^2) *
        whitening[["W"]] %*% m.cov.w %*% whitening[["W"]]

    # The inference pair (see the header).  The lambda-0 solve through the
    # released matrix, eigenvalues floored at the same public 1e-6 n as the
    # PD repair, so the pair stays defined when S_dp is not PD, in which
    # case the covariance is very large, which is the true state of affairs.
    # The de-shrink identity beta_w0 = beta_w + lambda.eff S_dp^{-1} beta_w
    # holds exactly wherever the floor does not bind, so bias.fold =
    # beta.infer - beta.dp coordinate-wise and the pair is post-processing.
    l.e0 <- eigen(m.S.dp, symmetric = TRUE)
    m.S0inv <- l.e0[["vectors"]] %*%
        (t(l.e0[["vectors"]]) / pmax(l.e0[["values"]], 1e-6 * s.n))
    v.beta.w0 <- as.numeric(m.S0inv %*% v.b.dp)
    v.beta.infer <- as.numeric(whitening[["W"]] %*% v.beta.w0) * y.std[["scale"]]
    v.beta.infer[[1L]] <- v.beta.infer[[1L]] + y.std[["centre"]]
    names(v.beta.infer) <- colnames(l.st[["m.X"]])
    v.bias.w <- v.beta.w0 - v.beta.w
    v.bias.raw <- as.numeric(whitening[["W"]] %*% v.bias.w) * y.std[["scale"]]
    m.cov.w0 <- m.S0inv %*% (s.sig.b^2 * diag(s.p) +
                             s.sig.S^2 * (sum(v.beta.w0^2) * diag(s.p) +
                                          tcrossprod(v.beta.w0))) %*% m.S0inv
    m.cov.infer <- (y.std[["scale"]]^2) *
        whitening[["W"]] %*% m.cov.w0 %*% whitening[["W"]]

    # Optional synthetic data under the public-X premise: the release is a
    # function of (beta_DP, sigma2_DP, X, seed) only, via the audited
    # synth.response constructor.
    d.syn <- NULL
    if(isTRUE(x.public)) {
        v.x.nms <- setdiff(names(d), y.nm)
        d.syn <- synth.response(dp.public.clip.x(d), y.nm, v.x.nms,
                                "gaussian", "identity",
                                v.beta.dp, s.sigma2.dp, seed)
    }

    list(d.syn       = d.syn,
         beta.dp     = v.beta.dp,
         sigma2.dp   = s.sigma2.dp,
         beta.w      = v.beta.w,
         S.dp        = m.S.dp,
         b.dp        = v.b.dp,
         q.dp        = s.q.dp,
         Cov.beta.dp = m.cov.raw,
         bias.fold   = v.bias.raw,
         beta.infer  = v.beta.infer,
         Cov.beta.infer = m.cov.infer,
         lambda      = s.lambda,
         lambda.eff  = s.lambda.eff,
         mu.total    = s.mu.total,
         mu.channels = sqrt(c(v.f, eig = s.f.eig)) * s.mu.total,
         sigma.channels = c(vech = s.sig.S, xty = s.sig.b, yy = s.sig.q),
         sens        = sens,
         n           = s.n, p = s.p,
         lambda.method = lambda.method)
}
