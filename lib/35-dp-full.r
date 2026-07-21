#!/usr/bin/env r
# +------------------------------------------------------------+
# |  Full (eps, delta)-DP synthetic data for sensitive X       |
# +------------------------------------------------------------+
#
# Extends the public-X DP-Projection variant (33-dp-projection-
# public-x.r) to settings where X itself is sensitive.  Releases
# DP marginals for each column of X (independent sampling), then
# uses the same full-DP pipeline for y:
#
#   1. DP marginal histograms for continuous X columns (B bins each).
#   2. DP proportions for categorical X columns.
#   3. Sample X_syn.  The continuous
#      block is drawn from a DP Gaussian copula (use.copula, default
#      TRUE), which restores its correlation structure.  Categorical
#      columns are still drawn independently.  use.copula = FALSE
#      recovers the legacy independent-marginal sampling.
#   4. DP beta_hat via Gaussian mechanism (as before).
#   5. DP sigma^2 via the residual-norm release: Laplace on
#      g = ||(I - H)y||_2 at sensitivity 2 B_y, debiased square.
#   6. e_syn = projected fresh Gaussian residuals.
#   7. y_syn = X_syn * beta_DP + e_syn.
#   8. Release (X_syn, y_syn).
#
# Privacy budget:
#   eps_total = sum(eps_X_marginals) + eps_beta + eps_sigma
#   via sequential composition.  For simulation scenarios with 4
#   X variables, the default split is 10% on X marginals, 80% on
#   beta, 10% on sigma.  Since beta_syn = beta_DP exactly, X_syn's
#   joint structure does not affect inferential fidelity, so eps_X
#   is starved in favour of eps_beta.  For dual-purpose releases
#   that also need descriptive fidelity of X, override to
#   eps.split = c(X = 0.5, beta = 0.4, sigma = 0.1).
#
# Trade-offs:
#   - The DP copula restores the continuous block's correlations, but
#     the categorical columns are still drawn independently, so the
#     categorical-to-continuous joint structure is not preserved.
#   - Either way inferential fidelity is unaffected: beta_syn = beta_DP
#     exactly, whatever joint structure X_syn carries.
#   - p-dimensional histograms scale: for high-p, the independence
#     assumption prevents curse-of-dimensionality in DP.

# Depends on: base R, plus dp.gm.multiplier() (28-dp-projection.r).
if(!exists("dp.gm.multiplier")) {
    .dp.suppress.demo <- TRUE; source("lib/28-dp-projection.r")
}

rlaplace <- function(n, location = 0, scale = 1) {
    v.u <- runif(n) - 0.5
    location - scale * sign(v.u) * log(1 - 2 * abs(v.u))
}

# +------------------------------------------------------------+
# |  DP marginal for one continuous variable                   |
# +------------------------------------------------------------+
# DP histogram: B bins, noisy counts, clipped to >= 0, renormalised.
# Sampling: uniform within each sampled bin.
#
# Two calibrations, exactly one of eps or mu:
#   mu   Gaussian mechanism, mu-GDP.  A replace-one row swap decreases one
#        count and increases another, so the count vector moves by (-1, +1)
#        and its L2 sensitivity is sqrt(2).  The noise sd is sqrt(2)/mu.  This
#        lets the marginal channel join the mu-GDP quadrature with the
#        correlation and coefficient releases.
#   eps  legacy Laplace(2/eps).  The L1 sensitivity is 2 under replace-one
#        (review C1).  Pure eps-DP, so it does not enter a mu-GDP quadrature.
#        Retained for the legacy sens.method = "local" path only.
#
# Args:
#   x          vector of observations
#   breaks     public bin boundaries (required for DP)
#   eps, mu    privacy budget for this marginal (supply exactly one)
dp.hist.continuous <- function(x, breaks, eps = NULL, mu = NULL) {
    if(is.null(eps) == is.null(mu))
        stop("supply exactly one of eps (Laplace) or mu (Gaussian, mu-GDP)")
    B <- length(breaks) - 1
    # The count vector must be the histogram, indexed by bin.  table() on an
    # integer vector drops the bins that no row happens to occupy, which would
    # left-pack the counts of the occupied bins and pad the zeros onto the end,
    # so an empty interior bin would slide every bin above the gap down one slot.
    #
    # That would be a privacy failure and not only a fidelity one, since the
    # sensitivity would then be computed for a function the code does not
    # release.  The sqrt(2) below is the L2 sensitivity of a
    # Histogram (a replace-one swap decreases one bin's count and increases
    # another's).  The left-packed vector has no such bound: moving one row into
    # an empty interior bin shifts every count above it, so its L2 sensitivity
    # is O(n).  Measured on an adversarial neighbouring pair at n = 2000 it is
    # 1413, against the 1.414 assumed here, a thousandfold under-noising of
    # this channel.
    #
    # factor(levels = seq_len(B)) keeps the empty bins, so the vector is the
    # histogram and the sqrt(2) bound is the true one.  Guarded by
    # verify-dp-soundness.R Task 5.
    v.counts <- as.numeric(table(factor(findInterval(x, breaks[-c(1, B + 1)]) + 1L,
                                        levels = seq_len(B))))
    v.noise <- if(is.null(mu)) rlaplace(B, 0, 2 / eps)
               else rnorm(B, 0, sqrt(2) / mu)
    v.counts.dp <- pmax(0, v.counts + v.noise)          # clip, then renormalise
    v.probs <- if(sum(v.counts.dp) > 0) v.counts.dp / sum(v.counts.dp)
               else rep(1 / B, B)
    list(breaks = breaks, probs = v.probs, B = B)
}

sample.from.dp.hist <- function(dp.hist, n) {
    B <- dp.hist[["B"]]; breaks <- dp.hist[["breaks"]]
    v.bin <- sample.int(B, n, replace = TRUE, prob = dp.hist[["probs"]])
    # uniform within bin
    v.lo <- breaks[v.bin]; v.hi <- breaks[v.bin + 1L]
    runif(n, v.lo, v.hi)
}

# Quantile function (inverse CDF) of a DP histogram: maps uniforms to
# values, so the copula sampler can attach each column's DP marginal.
qdphist <- function(u, dp.hist) {
    v.probs  <- dp.hist[["probs"]]
    v.breaks <- dp.hist[["breaks"]]
    B        <- dp.hist[["B"]]
    v.cum    <- cumsum(v.probs)
    u        <- pmin(pmax(u, 1e-12), 1 - 1e-12)
    v.b      <- pmin(findInterval(u, v.cum) + 1L, B)
    v.lo.cum <- c(0, v.cum)[v.b]
    v.frac   <- ifelse(v.probs[v.b] > 0, (u - v.lo.cum) / v.probs[v.b], 0.5)
    v.breaks[v.b] + v.frac * (v.breaks[v.b + 1L] - v.breaks[v.b])
}

# Forward CDF of a DP histogram: the piecewise-linear inverse of qdphist.  Maps
# a value to F_hat(x) in [0, 1] from the released DP bin probabilities, flat
# outside [breaks[1], breaks[B + 1]].  Backs the data-adaptive copula PIT, which
# transforms each covariate through its own already-released DP marginal.
pdphist <- function(x, dp.hist) {
    v.probs  <- dp.hist[["probs"]]
    v.breaks <- dp.hist[["breaks"]]
    B        <- dp.hist[["B"]]
    v.cum    <- c(0, cumsum(v.probs))                   # left-edge cum prob, length B+1
    v.xc     <- pmin(pmax(x, v.breaks[1]), v.breaks[B + 1L])
    v.b      <- pmax(pmin(findInterval(v.xc, v.breaks, rightmost.closed = TRUE), B), 1L)
    v.frac   <- (v.xc - v.breaks[v.b]) / (v.breaks[v.b + 1L] - v.breaks[v.b])
    v.cum[v.b] + v.probs[v.b] * v.frac
}

# +------------------------------------------------------------+
# |  DP marginal for one categorical variable                  |
# +------------------------------------------------------------+
# Same two calibrations as dp.hist.continuous: a replace-one row swap moves one
# unit between two levels, so the count vector has L1 sensitivity 2 (Laplace)
# and L2 sensitivity sqrt(2) (Gaussian, mu-GDP).
dp.cat <- function(x, levels, eps = NULL, mu = NULL) {
    if(is.null(eps) == is.null(mu))
        stop("supply exactly one of eps (Laplace) or mu (Gaussian, mu-GDP)")
    L <- length(levels)
    v.counts <- as.numeric(table(factor(x, levels = levels)))
    v.noise <- if(is.null(mu)) rlaplace(L, 0, 2 / eps)
               else rnorm(L, 0, sqrt(2) / mu)
    v.counts.dp <- pmax(0, v.counts + v.noise)
    v.probs <- if(sum(v.counts.dp) > 0) v.counts.dp / sum(v.counts.dp)
               else rep(1 / L, L)
    list(levels = levels, probs = v.probs)
}

sample.from.dp.cat <- function(dp.cat, n) {
    factor(sample(dp.cat[["levels"]], n, replace = TRUE,
                  prob = dp.cat[["probs"]]),
           levels = dp.cat[["levels"]])
}

# +------------------------------------------------------------+
# |  Gaussian-copula helpers (continuous covariate dependence) |
# +------------------------------------------------------------+
# Bounded normal scores for the Gaussian copula.  Three public PIT variants, all
# row-local and bounded by +- qnorm(1 - tau) so the correlation sensitivity is
# unchanged (it tracks the score bound, not the score spread):
#   "box"      box-uniform PIT (x - lo)/(hi - lo).  Distribution-agnostic, but
#              compresses concentrated data so the second-moment signal falls far
#              below the DP noise.
#   "normal"   location-scale standardisation with a public scale from the range,
#              sigma_pub = (hi - lo)/8.  Spreads the scores (variance ~1), lifting
#              the signal ~9x at no privacy cost.  Best for near location-scale
#              marginals.
#   "adaptive" normal-score PIT through the covariate's own released DP histogram,
#              score = qnorm(F_hat(x)).  Drops "normal"'s location-scale
#              assumption, so it also handles skewed and bounded marginals.  F_hat
#              is already public (the released marginal), so no extra privacy cost.
#              Requires dp.hist.  The default for the DP generators.
dp.copula.scores <- function(x, lo, hi, tau, method = c("box", "normal", "adaptive"),
                             dp.hist = NULL) {
    method <- match.arg(method)
    if(method == "adaptive") {
        if(is.null(dp.hist)) stop("score.method = \"adaptive\" requires dp.hist")
        return(qnorm(pmin(pmax(pdphist(x, dp.hist), tau), 1 - tau)))
    }
    if(method == "normal") {
        s.M <- qnorm(1 - tau)
        return(pmin(pmax((x - (lo + hi) / 2) / ((hi - lo) / 8), -s.M), s.M))
    }
    v.xc <- pmin(pmax(x, lo), hi)
    v.u  <- pmin(pmax((v.xc - lo) / (hi - lo), tau), 1 - tau)
    qnorm(v.u)
}

# DP release of the continuous covariates' correlation matrix via a Gaussian
# mechanism on the second-moment matrix of bounded scores.  Raw second-moment,
# not Fisher-z (whose sensitivity diverges near +-1).  S is not mean-centred,
# so recovered correlations are biased toward zero.
#
# Sensitivity.  What is released is vech(S), the upper triangle including the
# diagonal (see the rnorm below), not the full p x p matrix.  A Gaussian
# mechanism is mu-GDP at mu = Delta / sigma with Delta the L2 sensitivity of the
# Released vector, so the constant here is the sensitivity of vech(S):
#
#     Delta_2 = s.M^2 * p / n           (replace-one-row)
#
# and it is sharp: at the maximiser both swapped rows sit on the clamp and are
# orthogonal, which zeroes the diagonal of the move and attains the bound.
#
# The triangle-inequality bound on the Frobenius norm of the full matrix S is a
# different object from the one released, and it over-noises the channel by a
# factor of two.  Frobenius counts each off-diagonal entry twice where vech
# counts it once, so it can only ever over-noise.  A larger Delta means a larger
# sigma and hence a delivered mu smaller than the claimed one, which costs
# utility.
#
# Calibration: supply mu for the mu-GDP form (sigma = Delta_2 / mu), which is
# what composes in quadrature with the marginal and coefficient channels.
# Supplying (eps, delta) instead reproduces the legacy classical calibration.

# The single definition of this channel's sensitivity.  dp.copula.corr calibrates
# to it and verify-dp-soundness.R Task 5 measures against it, so the calibration
# and the test that guards it cannot drift apart.
dp.copula.corr.delta <- function(n, p, s.M) s.M^2 * p / n

dp.copula.corr <- function(m.W, eps, delta, s.M, seed = NA, mu = NULL) {
    if(!is.na(seed)) set.seed(seed)
    n <- nrow(m.W); p <- ncol(m.W)
    m.S <- crossprod(m.W) / n
    s.Delta <- dp.copula.corr.delta(n, p, s.M)
    s.sigma <- if(is.null(mu)) s.Delta * dp.gm.multiplier(eps, delta)
               else s.Delta / mu
    m.E <- matrix(0, p, p)
    v.idx <- upper.tri(m.S, diag = TRUE)
    m.E[v.idx] <- rnorm(sum(v.idx), 0, s.sigma)
    m.E[lower.tri(m.E)] <- t(m.E)[lower.tri(m.E)]
    m.S.dp <- m.S + m.E
    l.eig <- eigen((m.S.dp + t(m.S.dp)) / 2, symmetric = TRUE)
    v.lam <- pmax(l.eig[["values"]], 1e-8)
    m.psd <- l.eig[["vectors"]] %*% (v.lam * t(l.eig[["vectors"]]))
    v.d <- 1 / sqrt(diag(m.psd))
    m.R <- v.d * m.psd * rep(v.d, each = p)
    (m.R + t(m.R)) / 2
}

# Sample the continuous block from a DP Gaussian copula: correlated
# latent normals -> uniforms -> per-column DP-histogram quantiles.
sample.from.dp.copula <- function(m.R, l.dp.hist, n, seed = NA) {
    if(!is.na(seed)) set.seed(seed)
    p <- ncol(m.R)
    m.L <- chol(m.R + 1e-10 * diag(p))
    m.Z <- matrix(rnorm(n * p), n, p) %*% m.L
    m.U <- pnorm(m.Z)
    l.col <- Map(function(j) qdphist(m.U[, j], l.dp.hist[[j]]), seq_len(p))
    names(l.col) <- names(l.dp.hist)
    as.data.frame(l.col)
}

# +------------------------------------------------------------+
# |  Mixed-margin copula: categorical covariates                |
# +------------------------------------------------------------+
# Brings categorical covariates into the same Gaussian copula as the
# continuous block, so X_syn preserves their categorical-to-continuous joint
# structure.  This is a descriptive-fidelity feature only: beta_syn = beta_DP
# exactly regardless of X_syn's joint structure, so inferential fidelity is
# unaffected either way.  Assumptions in application are documented on
# gen.syn.dp.full below and mirrored in the manuscript (S6) and the package.
#
# Each categorical is encoded by the K-1 non-reference level indicators (level
# 1 is the reference).  Each indicator gets a bounded van der Waerden
# normal-score computed from the released (public) DP proportions, so it adds
# no privacy cost and lands in [-s.M, s.M].  These columns append to the
# continuous score matrix, so dp.copula.corr covers the whole block at
# sensitivity s.M^2 * p_total / n without change (it reads p = ncol(m.W), and
# the vech bound above holds for any clamped score matrix).  The bound is
# Attained only in the continuous sub-block: a categorical column takes one of
# K fixed winsorised values strictly inside the clamp, so its channel is
# over-noised in the safe direction (larger Delta => delivered mu <= claimed).
cat.indicator.scores <- function(x, dp.cat.obj, s.M) {
    v.levels <- dp.cat.obj[["levels"]]
    v.probs  <- dp.cat.obj[["probs"]]
    K <- length(v.levels)
    if(K < 2L) return(matrix(numeric(0), length(x), 0L))
    v.xc <- as.character(x)
    m.sc <- vapply(seq(2L, K), function(k) {
        p.k  <- v.probs[k]
        s.hi <- min(max(qnorm(1 - p.k / 2),   -s.M), s.M)   # rows AT level k
        s.lo <- min(max(qnorm((1 - p.k) / 2), -s.M), s.M)   # rows NOT at level k
        ifelse(v.xc == as.character(v.levels[k]), s.hi, s.lo)
    }, numeric(length(x)))
    matrix(m.sc, nrow = length(x))
}

# Back-map from a categorical's K-1 latents to a level, preserving the released
# marginal (via calibrated shifts) and inducing the copula association.
#   "probit"  level = argmax(0, Z_k + shift_k).  Deterministic, so it transfers
#             the full released latent correlation.  The natural inverse of a
#             Gaussian copula and consistent with the normal-scores used to fit
#             it.  Recommended.
#   "logit"   level ~ Categorical(softmax(0, Z_k + shift_k)).  Adds independent
#             choice noise, so the recovered association is attenuated (the
#             softmax step is orthogonal to the continuous latents).  Exposed as
#             an option.  See the assumptions note.
# Reference level = utility 0 (the leading column of zeros).  Neither back-map
# touches the DP guarantee: both are post-processing of the released noisy
# correlation.
.mixed.cat.props <- function(m.Z, v.shift, back.map) {
    m.U <- cbind(0, sweep(m.Z, 2L, v.shift, "+"))
    if(back.map == "probit") {
        tabulate(max.col(m.U, ties.method = "first"), ncol(m.U)) / nrow(m.U)
    } else {
        m.P <- exp(m.U - apply(m.U, 1L, max)); m.P <- m.P / rowSums(m.P)
        colMeans(m.P)
    }
}

.mixed.cat.calibrate <- function(m.Z.cal, v.probs, back.map) {
    v.shift <- log(pmax(v.probs[-1L], 1e-6) / max(v.probs[1L], 1e-6))
    for(it in seq_len(100L)) {
        v.p <- .mixed.cat.props(m.Z.cal, v.shift, back.map)
        if(max(abs(v.p - v.probs)) < 1e-4) break
        v.shift <- v.shift + 0.7 * log(pmax(v.probs[-1L], 1e-6) / pmax(v.p[-1L], 1e-6))
    }
    v.shift
}

.mixed.cat.assign <- function(m.Z, v.shift, v.levels, back.map) {
    m.U <- cbind(0, sweep(m.Z, 2L, v.shift, "+"))
    v.k <- if(back.map == "probit") {
        max.col(m.U, ties.method = "first")
    } else {
        m.P   <- exp(m.U - apply(m.U, 1L, max)); m.P <- m.P / rowSums(m.P)
        m.cum <- t(apply(m.P, 1L, cumsum))
        1L + rowSums(m.cum < runif(nrow(m.U)))
    }
    factor(v.levels[v.k], levels = v.levels)
}

# De-attenuation of the categorical-continuous cross entries.
# The fit measures the point-biserial correlation of a two-point indicator score
# against the continuous score, while the argmax back-map needs the latent
# (biserial) correlation.  For a normal threshold at released probability p_k the
# two differ by exactly
#
#     c_k = dnorm(qnorm(p_k)) / sqrt(p_k * (1 - p_k))
#
# (the truncated-first-moment identity E[1{U > tau} Y] = rho * phi(tau)), so the cross
# entries are divided by c_k before sampling.  c_k is a function of the released
# proportions only, so this is public post-processing: the released vech(S), its
# sensitivity and Task 7 are untouched (assumption A6 on gen.syn.dp.full).
# Applied to the categorical-continuous cross entries only.  Within- and
# cross-categorical blocks are left as fitted, since their exact inversion is
# the multinomial-probit MLE (follow-on work).  Validated on Scenario 1 at
# B = 100: the eps = 10 recovery of eta(race, income) rises from 63% to 90% of
# the real association, with the group-mean spread at 93%.
.mixed.deattenuate <- function(m.R, l.dp.X, v.cont.nms, v.cat.nms) {
    n.cont <- length(v.cont.nms)
    if(n.cont == 0L || length(v.cat.nms) == 0L) return(m.R)
    v.c <- unlist(Map(function(nm) {
        v.p <- l.dp.X[[nm]][["probs"]]
        if(length(v.p) < 2L) return(numeric(0))
        vapply(seq_len(length(v.p) - 1L) + 1L, function(k) {
            p.k <- min(max(v.p[k], 0.01), 0.99)
            dnorm(qnorm(p.k)) / sqrt(p.k * (1 - p.k))
        }, numeric(1))
    }, v.cat.nms))
    for(j in seq_along(v.c)) {
        col <- n.cont + j
        m.R[seq_len(n.cont), col] <- pmin(pmax(
            m.R[seq_len(n.cont), col] / v.c[[j]], -0.99), 0.99)
        m.R[col, seq_len(n.cont)] <- m.R[seq_len(n.cont), col]
    }
    # nearest-PSD + unit diagonal, the same post-processing dp.copula.corr ends on
    p <- ncol(m.R)
    l.eig <- eigen((m.R + t(m.R)) / 2, symmetric = TRUE)
    v.lam <- pmax(l.eig[["values"]], 1e-8)
    m.psd <- l.eig[["vectors"]] %*% (v.lam * t(l.eig[["vectors"]]))
    v.d <- 1 / sqrt(diag(m.psd))
    m.R2 <- v.d * m.psd * rep(v.d, each = p)
    (m.R2 + t(m.R2)) / 2
}

# Sample the whole X block from the mixed DP copula: correlated latents ->
# continuous through their DP-histogram quantiles, categorical through the
# calibrated probit/logit back-map.  Column order matches the extended m.W:
# continuous (in v.cont.nms order), then each categorical's K-1 columns.
# The latents are drawn from the de-attenuated correlation (T25b above).
sample.from.dp.mixed.copula <- function(m.R, l.dp.X, v.cont.nms, v.cat.nms,
                                        back.map, n, seed = NA) {
    if(!is.na(seed)) set.seed(seed)
    m.R <- .mixed.deattenuate(m.R, l.dp.X, v.cont.nms, v.cat.nms)
    p.total <- ncol(m.R)
    m.L <- chol(m.R + 1e-10 * diag(p.total))
    m.Z <- matrix(rnorm(n * p.total), n, p.total) %*% m.L
    d.out <- data.frame(.row = seq_len(n))
    j <- 0L
    for(nm in v.cont.nms) {
        j <- j + 1L
        d.out[[nm]] <- qdphist(pnorm(m.Z[, j]), l.dp.X[[nm]])
    }
    for(nm in v.cat.nms) {
        v.levels <- l.dp.X[[nm]][["levels"]]
        v.probs  <- l.dp.X[[nm]][["probs"]]
        K <- length(v.levels)
        v.cols <- seq.int(j + 1L, j + K - 1L)
        m.Z.cat <- m.Z[, v.cols, drop = FALSE]
        j <- j + (K - 1L)
        # calibrate shifts on a pool drawn from the released sub-correlation
        m.R.cc  <- m.R[v.cols, v.cols, drop = FALSE]
        m.Zc    <- matrix(rnorm(20000L * (K - 1L)), 20000L, K - 1L) %*%
                       chol(m.R.cc + 1e-10 * diag(K - 1L))
        v.shift <- .mixed.cat.calibrate(m.Zc, v.probs, back.map)
        d.out[[nm]] <- .mixed.cat.assign(m.Z.cat, v.shift, v.levels, back.map)
    }
    d.out[[".row"]] <- NULL
    d.out
}

# +------------------------------------------------------------+
# |  Full-DP synthesiser for sensitive X + sensitive y          |
# +------------------------------------------------------------+
# Wraps: DP X marginals + the mixed-margin DP copula + full-DP OLS y
# (beta, sigma).  Assumes gaussian identity link.  ZI and log-link deferred.
#
# Args:
#   d            data frame
#   y.nm         response column
#   x.spec       named list describing each X column:
#     list(age = list(type="continuous", breaks=seq(20, 85, length=11)),
#          female = list(type="categorical", levels=c("0", "1")),
#          race = list(type="categorical", levels=c("1","2","3")),
#          income = list(type="continuous", breaks=seq(0, 4000, length=11)))
#   eps.total    total privacy budget
#   delta        DP delta
#   B.y          public bound on |y| (legacy path only)
#   sens.method  "subagg" (default, sound mu-GDP) or "local" (deprecated, see below)
#   clip.lo/hi   public per-coefficient bounds.  Required for "subagg".
#   disp.lo/hi   public bounds on sigma^2.  Required for "subagg".
#   m            subagg block count
#   mu.split     fractions of mu^2 for the marginal, correlation and
#                coefficient channels (must sum to 1)
#   back.map     categorical back-map for the mixed copula, "probit" (default,
#                recommended) or "logit" (see assumptions A2)
#   seed         RNG seed
# Returns: list with d.syn, beta.dp, sigma2.dp, the mu budgets, etc.
#
# Mixed-margin copula, assumptions in application.
# With use.copula = TRUE the categorical covariates join the same Gaussian copula
# as the continuous block (each encoded by K-1 bounded normal-scored level
# indicators), so X_syn preserves their categorical-to-continuous joint
# structure.  This is descriptive fidelity only: beta_syn = beta_DP exactly
# regardless, so inferential fidelity is untouched.  The recovered joint carries
# the following assumptions, documented identically in the manuscript (S6) and
# the package (the synthexpln package, R/dp-full.R):
#   A1  Gaussian-copula fit: only the pairwise, latent-monotone association is
#       recovered.  Higher-order or non-elliptical dependence is not.
#   A2  Back-map: "probit" (argmax of the latents) is deterministic and transfers
#       the full released correlation.  "logit" (softmax then sample) injects
#       independent choice noise, so it attenuates the very association being
#       preserved.  Use it only when a stochastic sampler is specifically wanted.
#   A3  Sensitivity s.M^2 * p_total / n is a sound upper bound, attained only in
#       the continuous sub-block.  The categorical columns sit strictly inside the
#       clamp, so their channel is over-noised in the safe direction.
#   A4  p_total grows by K-1 per nominal, raising the correlation-channel noise
#       (acceptable: that channel is budget-starved, the coefficient dominates).
#   A5  Reference-level encoding: recovered associations are relative to level 1.
#   A6  De-attenuation (T25b): the categorical-continuous cross entries of the
#       released correlation are divided by the public biserial factor
#       c_k = dnorm(qnorm(p_k)) / sqrt(p_k (1 - p_k)) before sampling
#       (.mixed.deattenuate), correcting the point-biserial fit to the latent
#       scale the back-map needs.  Public post-processing (c_k reads the
#       Released proportions only), so the release and its sensitivity are
#       untouched.  Cross-categorical blocks are not corrected.
gen.syn.dp.full <- function(d, y.nm, x.spec, eps.total = 10,
                              delta = 1e-6, B.y = NULL,
                              sens.method = c("subagg", "local"),
                              clip.lo = NULL, clip.hi = NULL,
                              disp.lo = NULL, disp.hi = NULL, m = NULL,
                              mu.split = c(marg = 0.05, corr = 0.05, coef = 0.90),
                              eps.split = c(X = 0.1, beta = 0.8, sigma = 0.1),
                              use.copula = TRUE, eps.corr.frac = 0.5,
                              score.method = "adaptive",
                              back.map = c("probit", "logit"), seed = 6489) {
    sens.method <- match.arg(sens.method)
    score.method <- match.arg(score.method, c("box", "normal", "adaptive"))
    back.map <- match.arg(back.map)

    v.x.nms <- names(x.spec)
    if(!all(v.x.nms %in% names(d))) stop("x.spec names must match d")

    # Copula eligibility: at least two continuous columns
    v.is.cont <- vapply(v.x.nms,
                        function(nm) x.spec[[nm]][["type"]] == "continuous",
                        logical(1))
    v.cont.nms <- v.x.nms[v.is.cont]
    do.copula  <- isTRUE(use.copula) && length(v.cont.nms) >= 2L

    # ---- budget allocation --------------------------------------------------
    # subagg: every channel is a Gaussian mechanism, so the release composes in
    # One mu-GDP quadrature, mu^2 = mu_marg^2 + mu_corr^2 + mu_coef^2, converted
    # to (eps, delta) by the exact dual.  mu.split holds fractions of mu^2.  The
    # coefficient channel dominates by default: beta_syn = beta_DP exactly via
    # the projection, so X_syn's joint structure does not touch inferential
    # fidelity, and the marginal channels are starved in its favour.  Override
    # towards marg/corr for releases that also need descriptive fidelity of X.
    #
    # When the copula is off, the correlation share is folded back into the
    # marginals in quadrature, so the total mu is preserved either way.
    # The k marginal releases are k functions of the same row, so they too
    # compose in quadrature: each column gets mu_marg / sqrt(k).
    if(sens.method == "subagg") {
        if(is.null(clip.lo) || is.null(clip.hi))
            stop("sens.method = \"subagg\" requires the public per-coefficient ",
                 "bounds clip.lo and clip.hi")
        if(is.null(disp.lo) || is.null(disp.hi))
            stop("sens.method = \"subagg\" requires the public sigma^2 bounds ",
                 "disp.lo and disp.hi")
        if(abs(sum(mu.split) - 1) > 1e-8) stop("mu.split must sum to 1")
        s.mu.total <- gdp.mu.from.eps.delta(eps.total, delta)
        s.mu.coef  <- s.mu.total * sqrt(mu.split[["coef"]])
        s.mu.corr  <- if(do.copula) s.mu.total * sqrt(mu.split[["corr"]]) else 0
        s.mu.marg  <- if(do.copula) s.mu.total * sqrt(mu.split[["marg"]])
                      else gdp.compose(s.mu.total *
                                       sqrt(c(mu.split[["marg"]], mu.split[["corr"]])))
        s.mu.X.each <- s.mu.marg / sqrt(length(x.spec))
        s.eps.X.each <- NULL; s.eps.corr <- NULL
    } else {
        warning("sens.method = \"local\" calibrates beta to the DFBETA LOCAL ",
                "sensitivity, which is a function of the realised data. The ",
                "release therefore carries NO formal (eps, delta)-DP guarantee ",
                "Retained only for reproducibility. ",
                "Use sens.method = \"subagg\".")
        if(is.null(B.y)) stop("B.y required (public bound on |y|)")
        if(abs(sum(eps.split) - 1) > 1e-8) stop("eps.split must sum to 1")
        s.eps.X     <- eps.total * eps.split[["X"]]
        s.eps.beta  <- eps.total * eps.split[["beta"]]
        s.eps.sigma <- eps.total * eps.split[["sigma"]]
        s.eps.corr  <- if(do.copula) s.eps.X * eps.corr.frac else 0
        # Laplace(2/eps): the replace-one (Hamming-1) count sensitivity is 2
        # (review C1, remediation 1b).  Pure eps-DP, so it does not enter a
        # mu-GDP quadrature, which is one reason this path is deprecated.
        s.eps.X.each <- (s.eps.X - s.eps.corr) / length(x.spec)
        s.mu.total <- NA_real_; s.mu.coef <- NA_real_
        s.mu.corr  <- NA_real_; s.mu.marg <- NA_real_; s.mu.X.each <- NULL
    }

    # channel dispatch: the marginal mechanisms take exactly one of eps or mu.
    dp.marg.cont <- function(x, breaks)
        if(sens.method == "subagg") dp.hist.continuous(x, breaks, mu = s.mu.X.each)
        else dp.hist.continuous(x, breaks, eps = s.eps.X.each)
    dp.marg.cat <- function(x, levels)
        if(sens.method == "subagg") dp.cat(x, levels, mu = s.mu.X.each)
        else dp.cat(x, levels, eps = s.eps.X.each)

    # Fit real-data OLS to get beta_hat, sigma^2_hat
    fml <- reformulate(v.x.nms, response = y.nm)
    mo <- lm(fml, d)
    v.beta <- coef(mo)
    s.sigma2 <- summary(mo)[["sigma"]]^2
    m.X <- model.matrix(fml, d)
    n <- nrow(m.X); p <- ncol(m.X)

    # Stage 1: DP X marginals (+ optional mixed-margin Gaussian copula on X)
    set.seed(seed)
    l.dp.X <- list()
    if(do.copula) {
        v.cat.nms <- setdiff(v.x.nms, v.cont.nms)
        # DP marginals: histograms (continuous) and proportions (categorical),
        # all released.  The categoricals join the copula rather than
        # being sampled independently.
        for(nm in v.cont.nms)
            l.dp.X[[nm]] <- dp.marg.cont(d[[nm]], x.spec[[nm]][["breaks"]])
        for(nm in v.cat.nms)
            l.dp.X[[nm]] <- dp.marg.cat(d[[nm]], x.spec[[nm]][["levels"]])
        # Extended score matrix: continuous PIT scores, then K-1 categorical
        # indicator scores per categorical.  Every column is clamped to
        # [-s.M, s.M].  "normal"/"adaptive" winsorise at tau = 0.005.  "box"
        # clamps at the extreme quantile tau = 1/2n.  "adaptive" reads each
        # column's released DP histogram.
        s.tau <- if(score.method == "box") 1 / (2 * n) else 0.005
        s.M <- qnorm(1 - s.tau)
        m.W.cont <- vapply(v.cont.nms, function(nm) {
            v.bx <- range(x.spec[[nm]][["breaks"]])
            dp.copula.scores(d[[nm]], v.bx[1], v.bx[2], s.tau,
                             method = score.method, dp.hist = l.dp.X[[nm]])
        }, numeric(n))
        m.W <- cbind(m.W.cont,
                     do.call(cbind, lapply(v.cat.nms, function(nm)
                         cat.indicator.scores(d[[nm]], l.dp.X[[nm]], s.M))))
        # DP release of the extended correlation (p_total = continuous columns +
        # sum_c (K_c - 1) indicator columns), then one mixed copula sample.
        # dp.copula.corr reads p = ncol(m.W), so the released vech(S) sensitivity
        # is s.M^2 * p_total / n unchanged (the vech bound holds for any
        # clamped score matrix).  Attained only in the continuous sub-block.  See
        # the assumptions note on this function.
        m.R.dp <- if(sens.method == "subagg")
                      dp.copula.corr(m.W, NULL, delta, s.M, mu = s.mu.corr)
                  else dp.copula.corr(m.W, s.eps.corr, delta, s.M)
        d.syn <- sample.from.dp.mixed.copula(m.R.dp, l.dp.X, v.cont.nms,
                                             v.cat.nms, back.map, n)
        d.syn <- d.syn[v.x.nms]
    } else {
        m.R.dp <- NULL
        d.syn <- data.frame(.row = seq_len(n))
        for(nm in v.x.nms) {
            spec <- x.spec[[nm]]
            if(spec[["type"]] == "continuous") {
                l.dp.X[[nm]] <- dp.marg.cont(d[[nm]], spec[["breaks"]])
                d.syn[[nm]] <- sample.from.dp.hist(l.dp.X[[nm]], n)
            } else {
                l.dp.X[[nm]] <- dp.marg.cat(d[[nm]], spec[["levels"]])
                d.syn[[nm]] <- sample.from.dp.cat(l.dp.X[[nm]], n)
            }
        }
        d.syn[[".row"]] <- NULL
    }

    # ---- Stages 2 and 3: coefficient and dispersion --------------------------
    if(sens.method == "subagg") {
        # Subsample-and-aggregate on (beta, sigma^2) jointly.  Sound under
        # Full-row replacement, which is the relation that matters once X is
        # sensitive: the block partition is data-independent, and every
        # coordinate is clipped into a public box, so a one-row change moves the
        # mean by at most s_j = (hi_j - lo_j)/m whether it altered x, y or both.
        # The scaled aggregate has L2 sensitivity sqrt(p + 1), so the release is
        # mu_coef-GDP and composes in quadrature with the two X channels.
        #
        # This replaces both legacy channels, the DFBETA local-sensitivity beta
        # release (data-dependent, with no formal guarantee) and the Laplace
        # residual-norm sigma^2 release (valid only under response-row
        # replacement holding X fixed, which does not cover the sensitive-X case).
        l.sa <- subsample.aggregate.beta.loglink(
            d, y.nm, v.x.nms, "gaussian", "identity",
            eps = NA_real_, delta = delta, m = m,
            clip.lo = clip.lo, clip.hi = clip.hi,
            seed = seed + 10L, mu = s.mu.coef,
            disp = TRUE, disp.lo = disp.lo, disp.hi = disp.hi)
        v.beta.dp   <- l.sa[["beta.dp"]]
        s.sigma2.dp <- l.sa[["sigma2.dp"]]
    } else {
        l.sa <- NULL
        # Deprecated DFBETA local sensitivity.  Data-dependent, so it delivers
        # no formal (eps, delta)-DP.  Reachable only through
        # sens.method = "local", which warns.  Kept for reproducibility.
        v.sd <- apply(m.X, 2, stats::sd)
        v.s  <- ifelse(v.sd > 1e-8, v.sd, 1)
        m.dfb <- dfbeta(mo)
        s.LS  <- 2 * max(sqrt(rowSums(sweep(m.dfb, 2, v.s, "*")^2)))
        s.sigma.beta <- s.LS * dp.gm.multiplier(s.eps.beta, delta)
        set.seed(seed + 1L)
        v.eta.std <- rnorm(p, 0, s.sigma.beta)
        # 3-sigma truncation (post-processing)
        s.norm <- sqrt(sum(v.eta.std^2))
        s.max  <- 3 * sqrt(p) * s.sigma.beta
        if(s.norm > s.max) v.eta.std <- v.eta.std * (s.max / s.norm)
        v.beta.dp <- `names<-`((v.beta * v.s + v.eta.std) / v.s, names(v.beta))

        # Legacy sigma^2: Laplace on the residual norm g = ||(I - H)y||_2, whose
        # one-row y sensitivity is 2 B.y for any fixed design.  sigma^2 is the
        # debiased square, E[(g + L)^2] = g^2 + 2 b^2, so subtract 2 b^2.
        s.b <- 2 * B.y / s.eps.sigma
        set.seed(seed + 2L)
        s.g.dp <- sqrt(sum(residuals(mo)^2)) + rlaplace(1, 0, s.b)
        s.sigma2.dp <- max((s.g.dp^2 - 2 * s.b^2) / (n - p), 1e-6)
    }

    # Stage 4: construct y_syn using X_syn (not X_orig).
    # Use RHS-only formula: d.syn does not have y yet.
    fml.rhs <- reformulate(v.x.nms)
    m.X.syn <- model.matrix(fml.rhs, d.syn)
    # Some factor levels may be missing: align columns.
    v.common <- intersect(names(v.beta.dp), colnames(m.X.syn))
    v.beta.aligned <- setNames(rep(0, ncol(m.X.syn)), colnames(m.X.syn))
    v.beta.aligned[v.common] <- v.beta.dp[v.common]
    # Fresh residuals projected onto (I - H_syn) of X_syn
    m.H.syn <- m.X.syn %*% solve(crossprod(m.X.syn) + 1e-10 * diag(ncol(m.X.syn))) %*% t(m.X.syn)
    set.seed(seed + 3L)
    v.e.iid <- rnorm(n, 0, sqrt(s.sigma2.dp))
    v.e.proj <- as.numeric((diag(n) - m.H.syn) %*% v.e.iid)
    s.rss <- sum(v.e.proj^2)
    if(s.rss > 0)
        v.e.proj <- v.e.proj * sqrt((n - p) * s.sigma2.dp / s.rss)
    v.y.syn <- as.numeric(m.X.syn %*% v.beta.aligned) + v.e.proj

    d.syn[[y.nm]] <- v.y.syn

    # Not private, evaluation only: beta.hat and sigma2.hat are the full-data
    # fits, and subagg carries its own non-private intermediates (beta.agg,
    # n.failed).  The DP release is d.syn.  Do not compose the rest downstream.
    list(d.syn       = d.syn,
         beta.hat    = v.beta,
         beta.dp     = v.beta.dp,
         sigma2.hat  = s.sigma2,
         sigma2.dp   = s.sigma2.dp,
         sens.method = sens.method,
         eps.total   = eps.total,
         delta       = delta,
         # mu-GDP accounting (subagg).  mu.check re-composes the three channels
         # in quadrature and must reproduce mu.total.
         mu.total    = s.mu.total,
         mu.channels = c(marg = s.mu.marg, corr = s.mu.corr, coef = s.mu.coef),
         mu.check    = if(sens.method == "subagg")
                           gdp.compose(c(s.mu.marg, s.mu.corr, s.mu.coef))
                       else NA_real_,
         subagg      = l.sa,
         # legacy (sens.method = "local") diagnostics
         resid.norm.dp  = if(sens.method == "local") s.g.dp else NA_real_,
         eps.allocation = if(sens.method == "local")
                              c(X = s.eps.X, beta = s.eps.beta,
                                sigma = s.eps.sigma, corr = s.eps.corr)
                          else NULL,
         dp.X.marginals = l.dp.X,
         dp.corr        = m.R.dp,
         # Under subagg the whole release is calibrated to one (eps, delta) via
         # mu, so delta does not accumulate across channels as it did before.
         delta.total    = if(sens.method == "subagg") delta
                          else if(do.copula) 2 * delta else delta,
         use.copula     = do.copula)
}

# +------------------------------------------------------------+
# |  Monte Carlo driver                                         |
# +------------------------------------------------------------+
# This file defines a mechanism.  The Monte Carlo demonstration lives in
# 36-mc-all-variants.r, the single DP-variant driver, which takes every bound
# from dp-public-bounds.r and refuses to run on an invented one.
