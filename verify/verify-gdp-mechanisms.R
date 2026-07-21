#!/usr/bin/env r
# Assertion-based verification of the mu-GDP helpers and sound sensitivities
# of the DP mechanisms: 28-dp-projection.r and 33-dp-projection-public-x.r
# (Tasks 1-5), the ridge variance channel and exact ridge recovery of
# 34-dp-ridge-projection.r (Task 6), the legacy sensitive-X residual-norm
# release of 35-dp-full.r (Task 7, pinned so its numbers stay reproducible),
# and the sound mu-GDP sensitive-X releases, namely V4's three-channel
# quadrature (Task 8),
# the full-row replacement sensitivity of the joint (beta, sigma^2)
# subsample-and-aggregate (Task 9), and V5's logistic release (Task 10).
.dp.suppress.demo <- TRUE
.dp.public.no.main <- TRUE
source("lib/28-dp-projection.r")
source("lib/33-dp-projection-public-x.r")

# --- Task 1: mu-GDP helpers ----------------------------------------------
# Composition: two mu=1 mechanisms compose to sqrt(2).
stopifnot(abs(gdp.compose(c(1, 1)) - sqrt(2)) < 1e-12)

# Dual conversion is monotone increasing in mu at fixed eps.
stopifnot(gdp.to.delta(0.5, 1) < gdp.to.delta(1.0, 1))

# Round-trip: mu -> (eps, delta) -> mu recovers mu.
mu0 <- 1.3
d0 <- gdp.to.delta(mu0, 2)
stopifnot(abs(gdp.mu.from.eps.delta(2, d0) - mu0) < 1e-6)

# mu-GDP sigma matches the analytic multiplier for eps >= 1.
for (eps in c(1, 2, 5, 10)) {
    s.gdp <- gdp.sigma(1, eps, 1e-6)
    s.ana <- dp.gm.multiplier(eps, 1e-6)
    stopifnot(abs(s.gdp - s.ana) / s.ana < 1e-3)
}
cat("verify-gdp-mechanisms: mu-GDP helpers OK\n")

# --- Task 2: Route 1 global sensitivities --------------------------------
# The global sensitivity must upper-bound the actual worst-case beta change
# over every single-row y replacement in [-B.y, B.y].
set.seed(1)
m.X <- cbind(1, matrix(rnorm(30 * 2), 30, 2))
B.y <- 4
v.y <- pmin(pmax(as.numeric(m.X %*% c(1, 0.5, -0.3)) + rnorm(30), -B.y), B.y)
m.inv <- solve(crossprod(m.X))
b0 <- as.numeric(m.inv %*% crossprod(m.X, v.y))
gs <- gs.beta.ols.public(m.X, B.y)
worst <- 0
for (i in seq_len(30)) for (yr in c(-B.y, B.y)) {
    v.y2 <- v.y; v.y2[i] <- yr
    b2 <- as.numeric(m.inv %*% crossprod(m.X, v.y2))
    worst <- max(worst, sqrt(sum((b2 - b0)^2)))
}
stopifnot(worst <= gs + 1e-8)                 # sensitivity is a valid bound
stopifnot(gs > 0)

# The variance sensitivity must bound the worst-case sigma_hat^2 change over
# the whole bounded universe, not just the observed y.  The adversarial sign
# pattern y_j = -B.y sign(h_ij), y_i flipped between +/- B.y, attains
# 4 B.y^2 max_i max{1 - h_ii, s_i} / (n - p) whenever max s_i >= max(1 - h).
m.H <- m.X %*% m.inv %*% t(m.X)
v.h <- diag(m.H)
v.s <- rowSums(abs(m.H)) - abs(v.h)
gs.v <- gs.sigma2(m.X, B.y)
s2.of <- function(v) sum((v - m.H %*% v)^2) / (nrow(m.X) - ncol(m.X))
worst.v <- 0
for (i in seq_len(nrow(m.X))) {
    v.adv <- -B.y * sign(m.H[i, ]); v.adv[m.H[i, ] == 0] <- B.y
    v.y1 <- v.adv; v.y1[i] <- B.y
    v.y2 <- v.adv; v.y2[i] <- -B.y
    worst.v <- max(worst.v, abs(s2.of(v.y1) - s2.of(v.y2)))
}
stopifnot(worst.v <= gs.v + 1e-10)            # K_X bound is valid
if (max(v.s) >= max(1 - v.h))
    stopifnot(worst.v >= gs.v - 1e-8)         # and attained (exact constant)
# The naive 4 B.y^2 / (n - p) under-states the sensitivity whenever the
# hat-row constant K_X exceeds 1, which the assertion below covers.
if (max(pmax(1 - v.h, v.s)) > 1 + 1e-8)
    stopifnot(worst.v > 4 * B.y^2 / (nrow(m.X) - ncol(m.X)))
cat("verify-gdp-mechanisms: Route 1 sensitivities OK\n")

# --- Task 3: Route 1 Gaussianised release, exact recovery, mu accounting ---
set.seed(2)
n3 <- 400; p3 <- 3
X3 <- cbind(1, matrix(rnorm(n3 * p3), n3, p3))
y3 <- as.numeric(X3 %*% c(1, 0.5, -0.3, 0.2) + rnorm(n3))
d3 <- data.frame(y = y3, x1 = X3[, 2], x2 = X3[, 3], x3 = X3[, 4])
# B.y is a design constant of this test's own DGP, not a statistic of its sample.
#
# Here y = 1 + 0.5 x1 - 0.3 x2 + 0.2 x3 + N(0, 1) with x ~ N(0, 1) i.i.d., so
# y ~ N(1, sqrt(0.38 + 1)) = N(1, 1.175) and |y| <= 12 at 9.4 sd.  Public,
# derived from the DGP above, and independent of the realised draw.
B.y3 <- 12
res3 <- gen.syn.dp.ols.public(d3, y.nm = "y", x.nms = c("x1", "x2", "x3"),
                              eps = 10, delta = 1e-6, B.y = B.y3,
                              seed = 42)
stopifnot(!is.null(res3[["mu.total"]]), res3[["mu.total"]] > 0)
# mu.total equals the composition of the two Gaussian releases.
stopifnot(abs(res3[["mu.total"]] -
              gdp.compose(c(res3[["mu.beta"]], res3[["mu.v"]]))) < 1e-8)
# mu.total is calibrated to (eps = 10, delta = 1e-6).
stopifnot(abs(gdp.to.delta(res3[["mu.total"]], 10) - 1e-6) < 1e-8)
# Analyst OLS refit on synthetic data recovers beta.dp exactly.
fit3 <- lm(y ~ x1 + x2 + x3, res3[["d.syn"]])
stopifnot(max(abs(coef(fit3) - res3[["beta.dp"]])) < 1e-6)
cat("verify-gdp-mechanisms: Route 1 exact recovery + GDP accounting OK\n")

# --- Task 4: subsample-and-aggregate under mu-GDP ------------------------
# The scaled aggregate has L2 sensitivity sqrt(p), so sigma.Q equals
# sqrt(p) * dp.gm.multiplier for eps >= 1, and the release is mu-calibrated.
p4 <- 4
stopifnot(abs(gdp.sigma(sqrt(p4), 10, 1e-6) -
              sqrt(p4) * dp.gm.multiplier(10, 1e-6)) /
          gdp.sigma(sqrt(p4), 10, 1e-6) < 1e-3)
# The box below is a public constant, fixed in the source and independent of the
# sample.  See verify-dp-soundness.R Task 1, which asserts that property.
v.lo4 <- c(`(Intercept)` = -6, x1 = -6, x2 = -6, x3 = -6)
v.hi4 <- c(`(Intercept)` =  6, x1 =  6, x2 =  6, x3 =  6)
res4 <- gen.syn.dp.projected.subagg(
    d3, y.nm = "y", family = "gaussian", link = "identity",
    x.nms = c("x1", "x2", "x3"), eps = 10, delta = 1e-6, m = 20,
    clip.lo = v.lo4, clip.hi = v.hi4,
    dispersion = dp.dispersion.public(1, source = "design constant, unit error variance"),
    seed = 42)
stopifnot(!is.null(res4[["mu.total"]]), res4[["mu.total"]] > 0,
          is.finite(res4[["mu.total"]]))
stopifnot(abs(gdp.to.delta(res4[["mu.total"]], 10) - 1e-6) < 1e-8)
# The release must reproduce (beta_DP, sigma2) exactly: that is the projection
# theorem's exact-fidelity property, and it must survive drawing the response
# from the public family rather than warm-starting it from the real y.
fit4 <- lm(y ~ x1 + x2 + x3, res4[["d.syn"]])
stopifnot(max(abs(coef(fit4) - res4[["beta.dp"]])) < 1e-8)
stopifnot(abs(sigma(fit4)^2 - 1) < 1e-8)
cat("verify-gdp-mechanisms: Route 2 mu-GDP calibration OK\n")

# --- Task 5: subsample-and-aggregate composition for the DP-ZI generator ---
# The zero-inflated generator privatises a logistic zero-model at mu_Z and a
# log-link positive component at mu_P by subsample-and-aggregate, composing
# sqrt(mu_Z^2 + mu_P^2)-GDP.  Verify the extended mechanism (binomial family,
# supplied mu) is mu-calibrated and that the split composes to the dual.
set.seed(3)
n5 <- 800
X5 <- cbind(1, matrix(rnorm(n5 * 2), n5, 2))
z5 <- rbinom(n5, 1, plogis(as.numeric(X5 %*% c(-0.2, 0.6, -0.4))))
d.z5 <- data.frame(.Z = z5, x1 = X5[, 2], x2 = X5[, 3])
mu.dual <- gdp.mu.from.eps.delta(10, 1e-6)
mu.Z <- mu.dual * sqrt(0.5); mu.P <- mu.dual * sqrt(0.5)
# Again a public constant box, fixed in the source and independent of the sample.
v.lo5 <- c(`(Intercept)` = -4, x1 = -4, x2 = -4)
v.hi5 <- c(`(Intercept)` =  4, x1 =  4, x2 =  4)
l.z5 <- subsample.aggregate.beta.loglink(
    d.z5, ".Z", c("x1", "x2"), "binomial", "logit", 10, 1e-6, m = 20,
    clip.lo = v.lo5, clip.hi = v.hi5, seed = 42, mu = mu.Z)
stopifnot(abs(l.z5[["mu.total"]] - mu.Z) < 1e-8)       # binomial subagg at mu_Z
stopifnot(length(l.z5[["beta.dp"]]) == 3)
stopifnot(abs(gdp.compose(c(mu.Z, mu.P)) - mu.dual) < 1e-10)     # quadrature
stopifnot(abs(gdp.to.delta(gdp.compose(c(mu.Z, mu.P)), 10) - 1e-6) < 1e-8)
cat("verify-gdp-mechanisms: DP-ZI subagg composition (binomial + mu) OK\n")

# --- Task 6: ridge variance sensitivity (public X) ------------------------
# The ridge RSS is y'My with M = (I - H_l)'(I - H_l), H_l the ridge hat
# matrix X (X'X + lambda I)^{-1} X', so the Task-2 one-row identity holds
# with M in place of I - H and the worst case over the bounded universe is
# |Delta RSS| = 4 B.y^2 max_i max{M_ii, sum_{j != i} |M_ij|}.
# gs.sigma2.ridge must bound (and, at the adversarial sign pattern, attain)
# the brute-force worst case computed by actual ridge refits, and must
# reduce to gs.sigma2's K_X constant as lambda -> 0.
.dp.ridge.no.main <- TRUE
source("lib/34-dp-ridge-projection.r")
set.seed(4)
n6 <- 30L; lambda6 <- 1; B.y6 <- 4
m.X6 <- cbind(1, matrix(rnorm(n6 * 2), n6, 2))
m.X6[1, 2] <- 8                                   # leverage row pushes K > 1
p6 <- ncol(m.X6)
m.Ainv6 <- solve(crossprod(m.X6) + lambda6 * diag(p6))
m.M6 <- crossprod(diag(n6) - m.X6 %*% m.Ainv6 %*% t(m.X6))   # definitional M
v.m6 <- diag(m.M6)
v.s6 <- rowSums(abs(m.M6)) - abs(v.m6)
stopifnot(max(pmax(v.m6, v.s6)) > 1)              # naive constant must fail here
gs.v6 <- gs.sigma2.ridge(m.X6, lambda6, B.y6)
s2r.of <- function(v) {
    v.b <- as.numeric(m.Ainv6 %*% crossprod(m.X6, v))
    sum((v - m.X6 %*% v.b)^2) / (n6 - p6)
}
worst.v6 <- 0
for (i in seq_len(n6)) {
    v.adv <- -B.y6 * sign(m.M6[i, ]); v.adv[m.M6[i, ] == 0] <- B.y6
    v.y1 <- v.adv; v.y1[i] <- B.y6
    v.y2 <- v.adv; v.y2[i] <- -B.y6
    worst.v6 <- max(worst.v6, abs(s2r.of(v.y1) - s2r.of(v.y2)))
}
stopifnot(worst.v6 <= gs.v6 + 1e-10)              # bound is valid
stopifnot(worst.v6 >= gs.v6 - 1e-8)               # attained: max s_i > 1 >= max M_ii
stopifnot(worst.v6 > 4 * B.y6^2 / (n6 - p6))      # naive constant under-states
stopifnot(abs(gs.sigma2.ridge(m.X6, 1e-10, B.y6) - gs.sigma2(m.X6, B.y6)) <
          1e-6 * gs.sigma2(m.X6, B.y6))           # lambda -> 0 recovers K_X

# --- The ridge coefficient channel's global sensitivity ---------------
# The block above brute-forces the ridge variance sensitivity with real rigour
# and not once looks at the coefficient channel, which this block covers.  A
# local sensitivity, such as the ridge DFBETA read off the realised residuals,
# gives neighbouring datasets different noise scales, so no valid bound exists
# and there is no theorem to instantiate.
#
# Delta_beta = 2 B_y max_i ||A^{-1} x_i||, A = X'X + lambda I.  Brute-force it.
gs.b6 <- gs.beta.ridge.public(m.X6, lambda6, B.y6)
b.of <- function(v) as.numeric(m.Ainv6 %*% crossprod(m.X6, v))
set.seed(63)
v.y6.base <- pmin(pmax(as.numeric(m.X6 %*% c(1, 0.5, -0.3)) + rnorm(n6), -B.y6), B.y6)
worst.b6 <- max(vapply(seq_len(n6), function(i)
    max(vapply(c(-B.y6, B.y6), function(yr) {
        v2 <- v.y6.base; v2[i] <- yr
        sqrt(sum((b.of(v2) - b.of(v.y6.base))^2))
    }, numeric(1))), numeric(1)))
stopifnot(worst.b6 <= gs.b6 + 1e-9)               # the bound HOLDS
stopifnot(abs(gs.beta.ridge.public(m.X6, 1e-10, B.y6) -
              gs.beta.ols.public(m.X6, B.y6)) < 1e-6)   # lambda -> 0 is the OLS one

# A local sensitivity moves between neighbours, giving the two neighbours two
# noise scales and therefore no bound.  The assertion below confirms it moves.
v.sd6 <- apply(m.X6, 2, stats::sd); v.ss6 <- ifelse(v.sd6 > 1e-8, v.sd6, 1)
ls.of <- function(v) sens.beta.ridge(sweep(m.X6, 2, v.ss6, "/"), v, lambda6)[["LS"]]
# Flip the leverage row to the far corner.  Note it must be -B.y6 and not +B.y6:
# row 1 has x = 8, so its mean is about 5, and the clip above already pinned its
# response at +B.y6.  Writing +B.y6 back is a no-op, the "neighbour" is the same
# dataset, and the assertion passes vacuously.  (It did, on the first attempt.)
v.y6.nb <- `[<-`(v.y6.base, 1L, -B.y6)
stopifnot(!identical(v.y6.base, v.y6.nb))                 # a REAL neighbour
stopifnot(!isTRUE(all.equal(ls.of(v.y6.base), ls.of(v.y6.nb))))
# Exact ridge recovery: the analyst's ridge(lambda) refit on the synthetic
# data must return beta_DP to machine precision, which requires the
# lambda-corrected fitted component X (I + lambda (X'X)^{-1}) beta_DP.  A
# plain X beta_DP fitted part refits to (I - lambda A^{-1}) beta_DP only.
set.seed(9)
d6 <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
d6[["response"]] <- 1 + 0.5 * d6[["x1"]] - 0.3 * d6[["x2"]] + rnorm(80)
# Design constant again, not the sample maximum (see Task 3).  Here
# y = 1 + 0.5 x1 - 0.3 x2 + N(0, 1) with x ~ N(0, 1), so y ~ N(1, 1.158) and
# |y| <= 10 at 7.8 sd.
B.y6.pub <- 10
res6 <- gen.syn.dp.ridge.public(d6, "response", c("x1", "x2"), lambda = 1,
                                eps = 5, delta = 1e-6,
                                B.y = B.y6.pub,
                                seed = 11)
v.b6 <- fit.ridge(res6[["d.syn"]], "response", c("x1", "x2"), lambda = 1)
stopifnot(max(abs(v.b6 - res6[["beta.dp"]])) < 1e-8)
cat("verify-gdp-mechanisms: ridge variance sensitivity + exact recovery OK\n")

# --- Task 7: sensitive-X variance channel, residual-norm release -----------
# gen.syn.dp.full cannot use K_X (the design is private), so its sigma^2
# channel releases the residual norm g = ||(I - H)y||_2 at sensitivity
# 2 B.y: a one-row y change moves the norm by at most ||dy||_2 <= 2 B.y for
# any fixed design, the projection having operator norm 1, so no design
# constant appears.  sigma^2 is recovered by the debiased square
# (g_dp^2 - 2 b^2) / (n - p) as post-processing.  Checks: (i) the 2 B.y
# Lipschitz bound holds by brute force on actual projections, (ii) the
# release is seed-exact at the corrected scale with the debias, (iii) the
# released norm is returned for audit.
.dp.full.no.main <- TRUE
source("lib/35-dp-full.r")
set.seed(7)
n7 <- 200L
d7 <- data.frame(age = pmin(pmax(rnorm(n7, 50, 10), 1), 99))
d7[["y"]] <- 1 + 0.1 * d7[["age"]] + rnorm(n7)
x.spec7 <- list(age = list(type = "continuous",
                           breaks = seq(0, 100, length = 11)))
B.y7 <- 20
mo7 <- lm(y ~ age, d7)
m.X7 <- model.matrix(mo7)
m.H7 <- m.X7 %*% solve(crossprod(m.X7)) %*% t(m.X7)
g.of <- function(v) sqrt(sum(((diag(n7) - m.H7) %*% v)^2))
v.y7 <- d7[["y"]]
worst.g <- 0
for (i in seq_len(20)) {                      # flip row i between +/- B.y
    v.a <- v.y7; v.a[i] <- B.y7
    v.b <- v.y7; v.b[i] <- -B.y7
    worst.g <- max(worst.g, abs(g.of(v.a) - g.of(v.b)))
}
stopifnot(worst.g <= 2 * B.y7 + 1e-10)        # design-free Lipschitz bound
# sens.method = "local" is explicit here, pinning the legacy channel so its
# numbers stay reproducible.  It warns by design.
res7 <- suppressWarnings(
    gen.syn.dp.full(d7, "y", x.spec7, eps.total = 10, delta = 1e-6,
                    B.y = B.y7, sens.method = "local",
                    eps.split = c(X = 0.1, beta = 0.1, sigma = 0.8),
                    seed = 42))
s.g7 <- sqrt(sum(residuals(mo7)^2))
s.b7 <- 2 * B.y7 / 8                          # eps.sigma = 10 * 0.8
set.seed(42 + 2)
s.gdp7 <- s.g7 + rlaplace(1, 0, s.b7)
stopifnot(is.numeric(res7[["resid.norm.dp"]]),
          abs(res7[["resid.norm.dp"]] - s.gdp7) < 1e-12)
stopifnot(abs(res7[["sigma2.dp"]] -
              max((s.gdp7^2 - 2 * s.b7^2) / (n7 - 2L), 1e-6)) < 1e-12)
cat("verify-gdp-mechanisms: sensitive-X residual-norm release OK\n")

# --- Task 8: V4 sound sensitive-X release (mu-GDP, subagg) ----------------
# gen.syn.dp.full composes three Gaussian channels (X marginals, X
# correlations, and the joint (beta, sigma^2) subagg) in one mu-GDP quadrature.
source("lib/53-dp-logistic.r")
set.seed(6489)
n8 <- 2000
d8 <- data.frame(
    age    = runif(n8, 20, 90),
    income = rnorm(n8, 2000, 500),
    race   = factor(sample(c("1", "2", "3"), n8, TRUE), levels = c("1", "2", "3")))
d8[["y"]] <- 1 + 0.05 * d8[["age"]] + 0.001 * d8[["income"]] + rnorm(n8, 0, 2)
x.spec8 <- list(
    age    = list(type = "continuous", breaks = seq(20, 90, length = 11)),
    income = list(type = "continuous", breaks = seq(0, 4000, length = 11)),
    race   = list(type = "categorical", levels = c("1", "2", "3")))
v.cf8 <- c("(Intercept)", "age", "income", "race2", "race3")
lo8 <- `names<-`(c(-20, -2, -2, -10, -10), v.cf8)
hi8 <- `names<-`(c( 20,  2,  2,  10,  10), v.cf8)

res8 <- gen.syn.dp.full(d8, "y", x.spec8, eps.total = 10, delta = 1e-6,
                        clip.lo = lo8, clip.hi = hi8,
                        disp.lo = 0, disp.hi = 30, m = 20, seed = 6489)
stopifnot(res8[["sens.method"]] == "subagg")             # sound path is default
# the three channels re-compose in quadrature to mu.total ...
stopifnot(abs(res8[["mu.check"]] - res8[["mu.total"]]) < 1e-9)
# ... and mu.total is exactly the multiplier the (eps, delta) dual admits.
stopifnot(abs(gdp.to.delta(res8[["mu.total"]], 10) - 1e-6) < 1e-9)
# delta does not accumulate across channels any more (it was 2 * delta).
stopifnot(abs(res8[["delta.total"]] - 1e-6) < 1e-12)
# the projection invariant survives the new coefficient channel: the analyst
# refitting the synthetic data still recovers beta_DP exactly.
fit8 <- lm(reformulate(names(x.spec8), response = "y"), res8[["d.syn"]])
stopifnot(max(abs(coef(fit8) - res8[["beta.dp"]][names(coef(fit8))])) < 1e-6)
cat("verify-gdp-mechanisms: V4 sound sensitive-X release (mu quadrature) OK\n")

# --- Task 9: Full-row replacement sensitivity of the (beta, sigma^2) subagg ---
# A residual-norm sigma^2 channel is valid only under response-row replacement
# holding X fixed, which is the wrong relation once X is sensitive.
# The joint subagg is sound under full-row replacement: a data-independent
# partition plus a public clip box means one row moves the aggregate by at most
# s_j = (hi_j - lo_j)/m, whether it altered x, y, or both.  Same seed => same
# Gaussian draw, so the released differences are the aggregate differences.
sa.args9 <- list(y.nm = "y", x.nms = names(x.spec8), family = "gaussian",
                 link = "identity", eps = NA_real_, delta = 1e-6, m = 20,
                 clip.lo = lo8, clip.hi = hi8, seed = 6489, mu = 1,
                 disp = TRUE, disp.lo = 0, disp.hi = 30)
sa9a <- do.call(subsample.aggregate.beta.loglink, c(list(d = d8), sa.args9))
d9 <- d8                                     # one adversarial row: x AND y
d9[1, "age"]    <- 90
d9[1, "income"] <- 4000
d9[1, "race"]   <- factor("3", levels = c("1", "2", "3"))
d9[1, "y"]      <- 1e4                       # wildly outside the data range
sa9b <- do.call(subsample.aggregate.beta.loglink, c(list(d = d9), sa.args9))
v.s9  <- sa9a[["s.j"]]
v.mv9 <- abs(sa9b[["beta.dp"]] - sa9a[["beta.dp"]])
stopifnot(all(v.mv9 <= v.s9[names(v.mv9)] + 1e-9))         # every beta coord
stopifnot(abs(sa9b[["sigma2.dp"]] - sa9a[["sigma2.dp"]]) <=
          v.s9[[".sigma2"]] + 1e-9)                        # and the dispersion
cat("verify-gdp-mechanisms: full-row replacement sensitivity bound holds OK\n")

# --- Task 10: V5 sound sensitive-X logistic (S6 third extension) ----------
set.seed(6489)
n10 <- 4000
d10 <- data.frame(
    age    = runif(n10, 20, 90),
    income = rnorm(n10, 2000, 500),
    race   = factor(sample(c("1", "2", "3"), n10, TRUE), levels = c("1", "2", "3")))
d10[["y"]] <- rbinom(n10, 1,
    plogis(-2 + 0.03 * d10[["age"]] + 0.0005 * d10[["income"]]))
res10 <- gen.syn.dp.logistic(d10, "y", names(x.spec8), eps = 10, delta = 1e-6,
                             x.public = FALSE, x.spec = x.spec8, use.copula = TRUE,
                             clip.lo = `names<-`(c(-8, -1, -1, -5, -5), v.cf8),
                             clip.hi = `names<-`(c( 8,  1,  1,  5,  5), v.cf8),
                             m = 20, seed = 42)
stopifnot(res10[["sens.method"]] == "subagg")
stopifnot(abs(res10[["mu.check"]] - res10[["mu.total"]]) < 1e-9)
stopifnot(abs(gdp.to.delta(res10[["mu.total"]], 10) - 1e-6) < 1e-9)
stopifnot(abs(res10[["delta.total"]] - 1e-6) < 1e-12)      # no delta doubling
stopifnot(all(res10[["d.syn"]][["y"]] %in% c(0L, 1L)))     # Bernoulli tail
stopifnot(!isTRUE(all.equal(res10[["d.syn"]][["age"]], d10[["age"]])))  # X is synthetic
cat("verify-gdp-mechanisms: V5 sound sensitive-X logistic OK\n")

# --- Task 11: Route 3, the whitened sufficient-statistic release -----
# Mechanism checks only.  The sensitivity guard, the attainment, and the
# neighbour-invariance live in verify-dp-soundness.R Task 8.
source("lib/dp-public-bounds.r")
source("lib/66-dp-suffstat.r")
set.seed(6489)
n11 <- 500L
d11 <- data.frame(
    age    = runif(n11, 25, 80),
    female = factor(sample(c("0", "1"), n11, TRUE), levels = c("0", "1")),
    race   = factor(sample(c("1", "2", "3"), n11, TRUE), levels = c("1", "2", "3")),
    income = rnorm(n11, 2000, 400))
d11[["response"]] <- 50 + 0.5 * d11[["age"]] + 0.005 * d11[["income"]] +
    rnorm(n11, 0, 10)
l11.w <- dp.public.whitening(1L)
l11.y <- dp.public.y.std(1L)
l11.s <- dp.public.suffstat.sens(l11.w, l11.y)

res11 <- gen.syn.dp.suffstats(d11, "response", sc = 1L, eps = 10, delta = 1e-6,
                              whitening = l11.w, y.std = l11.y, sens = l11.s,
                              x.public = TRUE, seed = 42)
# The mu accounting: channels compose in quadrature to mu.total, and mu.total
# maps back to the requested (eps, delta) exactly.
stopifnot(abs(sqrt(sum(res11[["mu.channels"]]^2)) - res11[["mu.total"]]) < 1e-9)
stopifnot(abs(gdp.to.delta(res11[["mu.total"]], 10) - 1e-6) < 1e-9)
# Exact refit recovery under the public-X premise (the projection property).
v11.refit <- coef(lm(response ~ age + female + race + income, res11[["d.syn"]]))
stopifnot(max(abs(v11.refit - res11[["beta.dp"]])) < 1e-8)
# The inference pair.  beta.infer is
# the lambda-0 solve of the released statistics, bias.fold is the exact
# de-shrink beta.infer - beta.dp, and the recentred covariance is symmetric
# with a positive diagonal.  At eps = 10 the eigenvalue floor cannot bind,
# so the plain solve must reproduce beta.infer exactly.
stopifnot(max(abs(res11[["beta.infer"]] - res11[["beta.dp"]] -
                  res11[["bias.fold"]])) < 1e-9)
v11.b0 <- as.numeric(l11.w[["W"]] %*% solve(res11[["S.dp"]], res11[["b.dp"]])) *
    l11.y[["scale"]]
v11.b0[[1L]] <- v11.b0[[1L]] + l11.y[["centre"]]
stopifnot(max(abs(v11.b0 - res11[["beta.infer"]])) < 1e-8,
          max(abs(res11[["Cov.beta.infer"]] -
                  t(res11[["Cov.beta.infer"]]))) < 1e-8,
          all(diag(res11[["Cov.beta.infer"]]) > 0))
# Deterministic under the seed.
res11b <- gen.syn.dp.suffstats(d11, "response", sc = 1L, eps = 10, delta = 1e-6,
                               whitening = l11.w, y.std = l11.y, sens = l11.s,
                               x.public = TRUE, seed = 42)
stopifnot(identical(res11[["beta.dp"]], res11b[["beta.dp"]]))
# The public-dispersion mode drops the y'y channel and passes the declared
# value through, on the same quadrature.
l11.disp <- dp.dispersion.public(100, source = "verify task constant")
res11c <- gen.syn.dp.suffstats(d11, "response", sc = 1L, eps = 10, delta = 1e-6,
                               whitening = l11.w, y.std = l11.y, sens = l11.s,
                               dispersion = l11.disp, seed = 42)
stopifnot(res11c[["mu.channels"]][["yy"]] == 0,
          identical(res11c[["sigma2.dp"]], 100),
          is.na(res11c[["q.dp"]]),
          abs(sqrt(sum(res11c[["mu.channels"]]^2)) - res11c[["mu.total"]]) < 1e-9)
# PD repair: at a small budget the matrix noise dwarfs S, and the solve must
# still return finite coefficients through the raised effective lambda.
res11d <- gen.syn.dp.suffstats(d11, "response", sc = 1L, eps = 0.1, delta = 1e-6,
                               whitening = l11.w, y.std = l11.y, sens = l11.s,
                               seed = 7)
stopifnot(all(is.finite(res11d[["beta.dp"]])),
          res11d[["lambda.eff"]] >= res11d[["lambda"]])
# The released-eigenmin rule spends its share: the eig channel enters the
# quadrature and the total still resolves to (eps, delta).
res11e <- gen.syn.dp.suffstats(d11, "response", sc = 1L, eps = 10, delta = 1e-6,
                               whitening = l11.w, y.std = l11.y, sens = l11.s,
                               lambda.method = "released-eigenmin", seed = 42)
stopifnot(res11e[["mu.channels"]][["eig"]] > 0,
          abs(sqrt(sum(res11e[["mu.channels"]]^2)) - res11e[["mu.total"]]) < 1e-9)
# A generator does not invent its own bounds (by design).
stopifnot(inherits(tryCatch(gen.syn.dp.suffstats(d11, "response", sc = 1L),
                            error = function(e) e), "error"))
cat("verify-gdp-mechanisms: Route 3 suffstat release OK\n")

# --- Task 12: the public block-count rule ----------------------------
# m = floor(n_eff / (5 p)): the largest block count leaving five rows per
# released coefficient per block.  n_eff = n, or n * (1 - pi0) on a two-part
# positive channel with a published zero rate carried by a dp.rate provenance
# object.  A realised rate (mean(y == 0)) must not reach the rule.
stopifnot(dp.subagg.blocks(2000, 6) == 66L)
stopifnot(dp.subagg.blocks(2000, 6, dp.rate.public(0.5, "test rate")) == 33L)
stopifnot(inherits(tryCatch(dp.subagg.blocks(2000, 6, rate = 0.5),
                            error = function(e) e), "error"))
stopifnot(inherits(tryCatch(dp.rate.public(1.0, "test"),
                            error = function(e) e), "error"))
stopifnot(inherits(tryCatch(dp.rate.public(0.5, ""),
                            error = function(e) e), "error"))
stopifnot(inherits(tryCatch(dp.subagg.blocks(50, 6),
                            error = function(e) e), "error"))
# Engine self-default and the partition guard.  A caller that supplies a block
# assignment partitioned a fixed superset of the rows, so nrow(d) inside the
# engine is a restricted, sensitive count the rule must not read.
set.seed(1)
d.t12 <- data.frame(x1 = rnorm(300), x2 = rnorm(300))
d.t12[["y"]] <- 1 + d.t12[["x1"]] - d.t12[["x2"]] + rnorm(300)
l.t12 <- subsample.aggregate.beta.loglink(d.t12, "y", c("x1", "x2"),
    "gaussian", "identity", eps = 5, delta = 1e-6, m = NULL,
    clip.lo = -5, clip.hi = 5, seed = 1)
stopifnot(l.t12[["m"]] == 20)                    # floor(300 / (5 * 3)) = 20
v.blk <- sample(rep_len(1:10, 300))
stopifnot(inherits(tryCatch(
    subsample.aggregate.beta.loglink(d.t12, "y", c("x1", "x2"),
        "gaussian", "identity", eps = 5, delta = 1e-6, m = NULL,
        clip.lo = -5, clip.hi = 5, seed = 1, block = v.blk),
    error = function(e) e), "error"))
l.t12b <- subsample.aggregate.beta.loglink(d.t12, "y", c("x1", "x2"),
    "gaussian", "identity", eps = 5, delta = 1e-6, m = 10,
    clip.lo = -5, clip.hi = 5, seed = 1, block = v.blk)
stopifnot(l.t12b[["m"]] == 10)                   # explicit m always overrides
# The published-rate accessor serves a provenance object.
stopifnot(inherits(dp.public.zero.rate(4), "dp.rate"),
          inherits(tryCatch(dp.public.zero.rate(1),
                            error = function(e) e), "error"))
# ZI pair: m = NULL requires the published rate as a provenance object and
# yields per-channel block counts, floor(n / (5 p)) for the zero-model on all
# rows and floor(n (1 - pi0) / (5 p)) for the positive fit.
set.seed(2)
n.t12 <- 600L
d.zi <- data.frame(x1 = rnorm(n.t12), x2 = rnorm(n.t12))
v.pos.t12 <- rbinom(n.t12, 1L, 0.5) == 1L
d.zi[["response"]] <- ifelse(v.pos.t12, rgamma(n.t12, 2, 1), 0)
v.box.t12 <- c(`(Intercept)` = 3, x1 = 3, x2 = 3)
stopifnot(inherits(tryCatch(
    gen.syn.dp.projected.zi(d.zi, "response", family = "quasi-gamma",
        x.nms = c("x1", "x2"), eps = 5, m = NULL,
        clip.lo = -v.box.t12, clip.hi = v.box.t12, clip.z.lo = -v.box.t12,
        clip.z.hi = v.box.t12, dispersion = dp.dispersion.public(1, "test")),
    error = function(e) e), "error"))          # no pos.rate: must stop
l.zi.t12 <- gen.syn.dp.projected.zi(d.zi, "response", family = "quasi-gamma",
    x.nms = c("x1", "x2"), eps = 5, m = NULL,
    pos.rate = dp.rate.public(0.5, "test rate"),
    clip.lo = -v.box.t12, clip.hi = v.box.t12, clip.z.lo = -v.box.t12,
    clip.z.hi = v.box.t12, dispersion = dp.dispersion.public(1, "test"))
stopifnot(identical(l.zi.t12[["m"]], c(z = 40L, p = 20L)))
# floor(600 / 15) = 40 (zero-model), floor(600 * 0.5 / 15) = 20 (positive).
cat("verify-gdp-mechanisms: block-count rule OK\n")
