#!/usr/bin/env Rscript
# +------------------------------------------------------------------------+
# | 65-dp-suffstat-survey.r -- survey of the regularised statistic route    |
# +------------------------------------------------------------------------+
#
# What this supports.  The survey, which must settle the
# build order for the coefficient-channel gap-closer.  This script is decision
# Analysis, not a mechanism.  It releases nothing, adds no generator, and its
# numbers feed the survey record, not the paper.  If the route it assesses is
# built, the build gets its own scripts and verify tasks.
#
# The question.  20-dp-precheck.r rejected the Gaussian mechanism on X'y with
# an unregularised refit through (X'X)^{-1}, on the raw and unit-box scales,
# with X public.  The survey asks whether that rejection generalises to the joint
# regularised release
#
#     release  vech(X'X) + E   and   X'y + eta        (X private), then solve
#     beta_DP  = (X'X + E + lambda I)^{-1} (X'y + eta)
#
# on a better-conditioned public scale.  The precheck's own log shows the unit
# box is part of the failure it measures: income's public box [0, 4000] is 8x
# its published sd (500), so the unit-box column variance collapses and the
# smallest eigenvalue of X'X lands at 21.89 (n = 2000).  A public scale that
# standardises by the DGP's published parameters, or whitens by the population
# second-moment matrix, buys conditioning at the price of a larger B_x.  Which
# side wins is a computation, and this script is that computation.
#
# The three public scales.  Each is a fixed affine map whose constants are
# functionals of the published DGP, not statistics of the analysis sample.
#
#   unitbox   x_c -> (x_c - lo) / (hi - lo) by dp.public.xspec(), and y by its
#             public clip range.  B_x = sqrt(5) (dp.public.B.x.normalised).
#   std       x_c -> (x_c - m) / s by the DGP's published marginal parameters,
#             age (48, 13) and income (2000, 500), and y by (50, 15).  Bounds are
#             the xspec box mapped through the same constants.
#   white     x -> W x with W = Sigma_pop^{-1/2}, Sigma_pop = E[x x'] over an
#             independent draw of N = 1e6 at the public seed 20260713, the
#             64-dp-design-constants.r construction.  X'X/n is then near the
#             identity by design.  B_x is the maximum of ||W x|| over the
#             vertices of the public support box (a convex quadratic attains
#             its box maximum at a vertex).
#
# The sensitivities, one-row replacement throughout.
#
#   X'y     Delta_y = 2 B_x B_y, as in 20-dp-precheck.r.
#   vech    Delta_S = 2 B_x^2, the Frobenius bound on the move of X'X.  It
#           showed this over-noises the released vech by up to sqrt(2), so it
#           is sound and conservative here.  A route that wins under the
#           conservative constant wins a fortiori.  The exact vech constant
#           and its attainment proof are build-phase work.
#
# The lambda.  The ridge penalty is swept over a public grid, and the survey
# reports every cell plus the per-(scale, split, eps) minimum.  Choosing the
# minimising lambda from the population error curve is itself public under
# this project's premise (the DGP is published), the same argument that makes
# the subagg box public.  Real-data deployments without a published DGP would
# instead spend budget on a released eigenvalue, as AdaSSP does (Wang 2018),
# which is a build-phase concern and out of the survey's scope.
#
# Budget.  mu-GDP throughout, mu = 1 / dp.gm.multiplier(eps, delta), the same
# conversion as 62-mc-dp-gdp.r (eps {1, 2, 5, 10} -> mu {0.237, 0.448, 1.020,
# 1.848} at delta = 1e-6).  The private-X setting splits mu^2 between the two
# channels as mu_S = sqrt(f) mu and mu_y = sqrt(1 - f) mu.  The public-X
# setting spends everything on X'y.
#
# The metric.  Mean over noise draws of the max |% deviation| of the refitted
# coefficients from beta_hat, back-transformed exactly to the analyst's raw
# scale, over all 6 coefficients and again over the substantive filter
# |beta_hat| >= 0.01 that 62-mc-dp-gdp.r applies.  Comparators, quoted from
# the record rather than recomputed, the deprecated local path 12.3% (section 14,
# B = 25, no guarantee), subagg kappa = 1.5 / 3 at 29.2% / 21.5% (section 25,
# B = 200, under the block-count rule, and they read 34.0% / 66.7% at the
# superseded m = 20), Route 1 public-X 9059.5% (section 25), the unit-box
# precheck 484.5% analytic (results/dp-precheck.csv).
#
# Output: results/dp-suffstat-survey.csv, and the run log.
# +------------------------------------------------------------------------+

source("lib/00-shared.r")
source("lib/dp-gaussian-mechanism.r")   # dp.gm.multiplier, the ONE calibration
source("lib/dp-public-bounds.r")        # the ONLY permitted source of a DP bound

DELTA    <- 1e-6
V.EPS    <- c(1, 2, 5, 10)
V.FSPLIT <- c(0.3, 0.5, 0.7)                       # mu^2 fraction on the vech channel
V.LREL   <- c(0, 1e-4, 3e-4, 1e-3, 3e-3, 0.01, 0.03, 0.1, 0.3, 1)  # lambda / n
N.MC     <- 2000L                                  # noise draws per private-X cell
N.MC.PUB <- 10000L                                 # noise draws per public-X cell
MC.SEED  <- 20260716L                              # noise MC only, not the analysis seed
N.POP    <- 1e6L                                   # public-draw size, as 64- uses
PUB.SEED <- 20260713L                              # the PUBLIC seed of 64-

# +------------------------------------------------------------------------+
# | The analysis sample, publicly clipped (the estimand premise of every    |
# | route here, as in 20-dp-precheck.r)                                     |
# +------------------------------------------------------------------------+
load("results/d.sim.scenarios.RData")

prep <- function(d) {
    d[["female"]] <- factor(as.character(d[["female"]]), levels = c("0", "1"))
    d[["race"]]   <- factor(as.character(d[["race"]]),   levels = c("1", "2", "3"))
    d
}

d1 <- dp.public.clip.x(dp.public.clip.y(prep(l.scenarios[[1L]][["data"]]),
                                        "response", 1L))
mo    <- lm(response ~ age + female + race + income, d1)
v.beta <- coef(mo)
m.X   <- model.matrix(response ~ age + female + race + income, d1)
v.y   <- d1[["response"]]
s.n   <- nrow(m.X)
s.p   <- length(v.beta)
v.substantive <- abs(v.beta) >= 0.01               # the 62- filter

# +------------------------------------------------------------------------+
# | The public whitening matrix, from the 64- construction                  |
# +------------------------------------------------------------------------+
# Lift sim.scenario() out of 00-simulate.r without executing its main block,
# exactly as 64-dp-design-constants.r does.  The population draw does not see
# the seed-6489 analysis sample.
v.exprs <- parse("simulate/00-simulate.r")
b.simdef <- Map(function(e) {
    is.call(e) && identical(as.character(e[[1L]]), "<-") &&
        identical(as.character(e[[2L]]), "sim.scenario")
}, v.exprs) |> unlist()
invisible(Map(function(e) eval(e, envir = globalenv()), v.exprs[b.simdef]))
stopifnot(is.function(sim.scenario))

d.pop <- dp.public.clip.x(prep(sim.scenario(1L, n = N.POP, seed = PUB.SEED)[["data"]]))
m.X.pop <- model.matrix(~ age + female + race + income, d.pop)
m.Sigma.pop <- crossprod(m.X.pop) / N.POP
l.eig <- eigen(m.Sigma.pop, symmetric = TRUE)
m.W <- l.eig[["vectors"]] %*% diag(1 / sqrt(l.eig[["values"]])) %*% t(l.eig[["vectors"]])
rm(d.pop, m.X.pop)

# +------------------------------------------------------------------------+
# | The scale maps                                                          |
# +------------------------------------------------------------------------+
# Each scale is x_s = M x on the augmented coefficient basis (intercept, age,
# female1, race2, race3, income) plus a public affine map of y.  The refit on
# the mapped data returns beta_s = b^{-1} (M')^{-1} beta, so the exact back
# map is beta = b * t(M) %*% beta_s, asserted below on the noiseless refit.
#
# B_x is the max of ||M x|| over the public support box, which for unitbox and
# std is available in closed form and for white is a vertex enumeration (the
# box maximum of a convex quadratic sits at a vertex.  The categoricals
# contribute their finite dummy patterns).
v.xspec.age    <- range(dp.public.xspec()[["age"]][["breaks"]])      # [20, 90]
v.xspec.income <- range(dp.public.xspec()[["income"]][["breaks"]])   # [0, 4000]
l.By1 <- dp.public.B.y(1L)                                           # clip [-40, 140]

mk.M.diag <- function(v.shift, v.scale) {
    # x_s = (x - shift) / scale coordinate-wise.  The shift is absorbed through
    # the constant column of x, so it lands in the first column of M: row j
    # reads x_s[j] = x[j] / scale_j - (shift_j / scale_j) * x[1] with x[1] = 1.
    m.M <- diag(1 / v.scale)
    m.M[, 1L] <- -v.shift / v.scale
    m.M[1L, 1L] <- 1
    dimnames(m.M) <- list(names(v.beta), names(v.beta))
    m.M
}

# Vertices of the public support box on the augmented basis.
m.box.vertices <- local({
    l.g <- expand.grid(age    = v.xspec.age,
                       female = c(0, 1),
                       race   = c("1", "2", "3"),
                       income = v.xspec.income)
    m.V <- cbind(1, l.g[["age"]], l.g[["female"]],
                 as.integer(l.g[["race"]] == "2"),
                 as.integer(l.g[["race"]] == "3"),
                 l.g[["income"]])
    colnames(m.V) <- names(v.beta)
    m.V
})

B.x.of <- function(m.M, source)
    list(B.x = sqrt(max(rowSums((m.box.vertices %*% t(m.M))^2))), source = source)

l.scales <- list(
    unitbox = local({
        v.shift <- c(0, v.xspec.age[[1L]], 0, 0, 0, v.xspec.income[[1L]])
        v.scale <- c(1, diff(v.xspec.age), 1, 1, 1, diff(v.xspec.income))
        list(m.M = mk.M.diag(v.shift, v.scale),
             y.a = l.By1[["y.clip"]][[1L]], y.b = diff(l.By1[["y.clip"]]),
             l.Bx = dp.public.B.x.normalised(2L, 2L),
             l.By = list(B.y = 1, source = paste0(
                 "1 by the PUBLIC clip of dp.public.B.y(1): y in [-40, 140] ",
                 "rescaled to [0, 1]. Not by min-max scaling on the sample.")))
    }),
    std = local({
        v.shift <- c(0, 48, 0, 0, 0, 2000)
        v.scale <- c(1, 13, 1, 1, 1, 500)
        m.M <- mk.M.diag(v.shift, v.scale)
        list(m.M = m.M, y.a = 50, y.b = 15,
             l.Bx = B.x.of(m.M, paste0(
                 "max ||x||_2 over the public xspec box standardised by the ",
                 "DGP's published marginal parameters, age (48, 13) and income ",
                 "(2000, 500). A functional of the published DGP, not of the ",
                 "analysis sample.")),
             l.By = list(B.y = 6, source = paste0(
                 "the public clip [-40, 140] standardised by the DGP marginal ",
                 "(50, 15), so |y_s| <= 6.")))
    }),
    white = list(m.M = m.W, y.a = 50, y.b = 15,
        l.Bx = B.x.of(m.W, paste0(
            "max ||Sigma_pop^{-1/2} x||_2 over the vertices of the public ",
            "xspec box, with Sigma_pop = E[x x'] under the published DGP at ",
            "N = 1e6, public seed 20260713 (the 64-dp-design-constants.r ",
            "construction). A functional of the published DGP.")),
        l.By = list(B.y = 6, source = paste0(
            "the public clip [-40, 140] standardised by the DGP marginal ",
            "(50, 15), so |y_s| <= 6."))))

# +------------------------------------------------------------------------+
# | The error engine                                                        |
# +------------------------------------------------------------------------+
pct.dev <- function(m.beta.dp) {
    # rows = draws, cols = coefficients, on the raw analyst scale
    m.rel <- 100 * abs(sweep(m.beta.dp, 2L, v.beta, `-`)) /
        rep(pmax(abs(v.beta), 1e-10), each = nrow(m.beta.dp))
    c(e.max.pct      = mean(apply(m.rel, 1L, max)),
      e.max.pct.subst = mean(apply(m.rel[, v.substantive, drop = FALSE], 1L, max)),
      e.mean.pct     = mean(m.rel),
      p.gt50         = mean(apply(m.rel, 1L, max) > 50))
}

l.rows <- list()
emit <- function(setting, scale, f.S, eps, mu, lambda.rel, l.Bx, l.By, v.stat,
                 n.fail = 0L)
    l.rows[[length(l.rows) + 1L]] <<- data.frame(
        setting = setting, scale = scale, f.S = f.S, eps = eps, mu = mu,
        lambda.rel = lambda.rel, B.x = l.Bx[["B.x"]], B.y = l.By[["B.y"]],
        n = s.n, p = s.p, delta = DELTA,
        e.max.pct = v.stat[["e.max.pct"]],
        e.max.pct.subst = v.stat[["e.max.pct.subst"]],
        e.mean.pct = v.stat[["e.mean.pct"]], p.gt50 = v.stat[["p.gt50"]],
        n.fail = n.fail, stringsAsFactors = FALSE)

cat(strrep("=", 70), "\n")
cat("SURVEY -- regularised sufficient-statistic release, scenario 1\n")
cat(strrep("=", 70), "\n")
cat(sprintf("n = %d, p = %d, substantive filter keeps %d of %d coefficients\n",
            s.n, s.p, sum(v.substantive), s.p))

invisible(Map(function(s.scale) {
    l.s <- l.scales[[s.scale]]
    m.M <- l.s[["m.M"]]

    m.Xs <- m.X %*% t(m.M)
    v.ys <- (v.y - l.s[["y.a"]]) / l.s[["y.b"]]
    m.A  <- crossprod(m.Xs)
    v.b  <- crossprod(m.Xs, v.ys)

    # The exact back map, asserted before it is used on anything noisy.  The
    # y shift is absorbed by the intercept, whose column of X is the constant 1.
    back <- function(m.beta.s) {
        m.out <- m.beta.s %*% m.M * l.s[["y.b"]]
        m.out[, 1L] <- m.out[, 1L] + l.s[["y.a"]]
        m.out
    }
    stopifnot(max(abs(back(t(solve(m.A, v.b))) - rbind(v.beta))) < 1e-6)

    v.eig <- eigen(m.A, symmetric = TRUE, only.values = TRUE)[["values"]]
    cat(sprintf("\n--- scale %s ---\n", s.scale))
    dp.bounds.announce(B.x = l.s[["l.Bx"]], B.y = l.s[["l.By"]])
    cat(sprintf("  eigenvalues of X'X: min %.4g, max %.4g, condition %.4g\n",
                min(v.eig), max(v.eig), max(v.eig) / min(v.eig)))

    s.Bx <- l.s[["l.Bx"]][["B.x"]]
    s.By <- l.s[["l.By"]][["B.y"]]
    s.Delta.y <- 2 * s.Bx * s.By
    s.Delta.S <- 2 * s.Bx^2

    # -------- public X: noise on X'y only, ridge grid ---------------------
    invisible(Map(function(s.eps) {
        s.mu <- 1 / dp.gm.multiplier(s.eps, DELTA)
        s.sigma <- s.Delta.y / s.mu
        set.seed(MC.SEED)
        m.eta <- matrix(rnorm(N.MC.PUB * s.p, sd = s.sigma), ncol = s.p)
        invisible(Map(function(s.lrel) {
            m.solve <- solve(m.A + s.lrel * s.n * diag(s.p))
            # beta_s(lambda) = (A + lambda n I)^{-1} (b + eta): the ridge bias
            # sits in m.solve %*% v.b and the noise in m.solve %*% eta.  The
            # length-p vector recycles down the columns of the p x N matrix,
            # which is the correct orientation.
            m.beta.s <- t(as.numeric(m.solve %*% v.b) + m.solve %*% t(m.eta))
            emit("public-x", s.scale, NA_real_, s.eps, s.mu, s.lrel,
                 l.s[["l.Bx"]], l.s[["l.By"]], pct.dev(back(m.beta.s)))
        }, V.LREL))
    }, V.EPS))

    # -------- private X: joint vech(X'X) and X'y release, ridge grid ------
    invisible(Map(function(s.eps) {
        s.mu <- 1 / dp.gm.multiplier(s.eps, DELTA)
        invisible(Map(function(s.f) {
            s.sigma.S <- s.Delta.S / (sqrt(s.f) * s.mu)
            s.sigma.y <- s.Delta.y / (sqrt(1 - s.f) * s.mu)
            set.seed(MC.SEED)
            # One draw set per (scale, split, eps), reused across the lambda
            # grid through a single eigendecomposition per draw.
            l.acc <- Map(function(s.l) matrix(NA_real_, N.MC, s.p), V.LREL)
            v.fail <- integer(length(V.LREL))
            for(i in seq_len(N.MC)) {
                m.E <- matrix(0, s.p, s.p)
                m.E[upper.tri(m.E, diag = TRUE)] <- rnorm(s.p * (s.p + 1) / 2,
                                                          sd = s.sigma.S)
                m.E <- m.E + t(m.E) - diag(diag(m.E))
                v.bn <- as.numeric(v.b) + rnorm(s.p, sd = s.sigma.y)
                l.e <- eigen(m.A + m.E, symmetric = TRUE)
                v.proj <- as.numeric(crossprod(l.e[["vectors"]], v.bn))
                for(k in seq_along(V.LREL)) {
                    v.d <- l.e[["values"]] + V.LREL[[k]] * s.n
                    if(min(v.d) <= 1e-8) { v.fail[[k]] <- v.fail[[k]] + 1L; next }
                    l.acc[[k]][i, ] <- as.numeric(l.e[["vectors"]] %*%
                                                  (v.proj / v.d))
                }
            }
            invisible(Map(function(s.lrel, m.acc, s.nf) {
                m.ok <- m.acc[stats::complete.cases(m.acc), , drop = FALSE]
                if(nrow(m.ok) == 0L) return(invisible(NULL))
                emit("private-x", s.scale, s.f, s.eps, s.mu, s.lrel,
                     l.s[["l.Bx"]], l.s[["l.By"]], pct.dev(back(m.ok)), s.nf)
            }, V.LREL, l.acc, as.list(v.fail)))
        }, V.FSPLIT))
    }, V.EPS))
}, names(l.scales)))

d.out <- do.call(rbind, l.rows)
write.csv(d.out, "results/dp-suffstat-survey.csv", row.names = FALSE)

# +------------------------------------------------------------------------+
# | The headline: best lambda per (setting, scale, split, eps)              |
# +------------------------------------------------------------------------+
cat("\n", strrep("=", 70), "\n", sep = "")
cat("BEST CELL PER (setting, scale, split) AT eps = 10, over the lambda grid\n")
cat("(e.max.pct over all 6 coefficients; subst applies the 62- filter)\n\n")
d.10 <- d.out[d.out[["eps"]] == 10, ]
v.key <- paste(d.10[["setting"]], d.10[["scale"]],
               ifelse(is.na(d.10[["f.S"]]), "-", d.10[["f.S"]]))
l.key <- split(d.10, v.key)
invisible(Map(function(d.g) {
    d.b <- d.g[which.min(d.g[["e.max.pct"]]), ]
    cat(sprintf("  %-9s %-8s f.S=%-4s lambda/n=%-6g  max%% %8.1f  subst %8.1f  mean%% %6.1f  fail %d\n",
                d.b[["setting"]], d.b[["scale"]],
                format(d.b[["f.S"]]), d.b[["lambda.rel"]],
                d.b[["e.max.pct"]], d.b[["e.max.pct.subst"]],
                d.b[["e.mean.pct"]], d.b[["n.fail"]]))
}, l.key))

cat("\nComparators, quoted from the recorded values (not recomputed here):\n")
cat("  deprecated local path   12.3%  (s14, B = 25, no guarantee)\n")
cat("  subagg kappa = 1.5      29.2%  (s25, B = 200, sound)\n")
cat("  subagg kappa = 3        21.5%  (s25, B = 200, sound)\n")
cat("  Route 1 public-X      9059.5%  (s25, B = 200, sound)\n")
cat("  precheck unitbox X'y   484.5%  (dp-precheck.csv, analytic)\n")
cat("\nWritten to dp-suffstat-survey.csv\n")
cat("Copy every reported number FROM THAT FILE, with this script named.\n")
