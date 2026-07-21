#!/usr/bin/env Rscript
# +------------------------------------------------------------------------+
# |  verify-dp-soundness.R -- sensitivity invariants for the DP releases    |
# +------------------------------------------------------------------------+
#
# What this file checks, and why it is separate from verify-gdp-mechanisms.R.
#
# verify-gdp-mechanisms.R checks the mu-GDP arithmetic, meaning composition,
# the exact dual, and calibration.  Correct arithmetic over an invalid
# sensitivity still gives an invalid release, so this file asserts the property
# that the arithmetic assumes:
#
#   The noise is calibrated to the sensitivity of the function the code
#   Actually releases.
#
# That property has two requirements.  Every bound must be a functional of
# public quantities rather than of the sensitive sample, and the bound must be
# taken over the object the mechanism releases rather than over a related one.
#
# The tasks below assert both, from several directions.
#
# Run:  Rscript verify-dp-soundness.R    (expects "OK" from every task, exit 0)
#
# +------------------------------------------------------------------------+
.sobol.main <- TRUE
source("lib/00-shared.r")
source("lib/gen-syn-family.r")
source("lib/07-calibrated-generator.r")
source("lib/09-symmetric-generator.r")
source("lib/11-projection-generator.r")
source("lib/28-dp-projection.r")
.dp.public.no.main <- TRUE
source("lib/33-dp-projection-public-x.r")
.dp.full.no.main <- TRUE
source("lib/35-dp-full.r")                                  # V4's channels, for Task 5
rm(.sobol.main, .dp.public.no.main, .dp.full.no.main)

ok <- function(msg) cat(sprintf("verify-dp-soundness: %s OK\n", msg))

skew <- function(v) mean((v - mean(v))^3) / sd(v)^3
kurt <- function(v) mean((v - mean(v))^4) / sd(v)^4 - 3

# A public box.  Fixed constants, chosen from the design, not from a sample.
V.LO <- c(`(Intercept)` = -8, x1 = -8, x2 = -8)
V.HI <- c(`(Intercept)` =  8, x1 =  8, x2 =  8)

# +------------------------------------------------------------------------+
# |  Task 1 -- the noise calibration must not move between neighbours       |
# +------------------------------------------------------------------------+
# Catches a noise calibration that varies with the sample.
#
# This is the executable form of the sensitivity hypothesis, which binds the
# clip box once and applies the same box to both neighbouring datasets.  A box
# derived from the data gives the two neighbours different boxes, so the
# theorem cannot be
# instantiated: there is no valid sensitivity bound, hence no valid noise
# calibration, hence no theorem -- not a weaker epsilon, and not any epsilon.
#
# Operationally: s_j and sigma_j are the mechanism's calibration.  Build D and D'
# one row apart (adversarially, to the corners).  If the calibration moves at
# all, it is a function of the private data.  Tolerance is exactly zero.
set.seed(11L)
n1 <- 600L
x1a <- rnorm(n1); x1b <- rnorm(n1)
d1  <- data.frame(response = rgamma(n1, shape = 2,
                    rate = 2 / exp(1.2 + 0.5 * x1a - 0.3 * x1b)),
                  x1 = x1a, x2 = x1b)
d1.nb <- d1
d1.nb[1L, "x1"]       <- max(d1[["x1"]]) * 3
d1.nb[1L, "x2"]       <- min(d1[["x2"]]) * 3
d1.nb[1L, "response"] <- max(d1[["response"]]) * 5

sa.at <- function(dd) subsample.aggregate.beta.loglink(dd, "response",
    c("x1", "x2"), "gamma", "log", eps = 10, delta = 1e-6, m = 10,
    clip.lo = V.LO, clip.hi = V.HI, seed = 7L)
a1 <- sa.at(d1); a2 <- sa.at(d1.nb)

stopifnot(identical(a1[["s.j"]],     a2[["s.j"]]))       # sensitivity
stopifnot(identical(a1[["sigma"]],   a2[["sigma"]]))     # noise sd
stopifnot(identical(a1[["clip.lo"]], a2[["clip.lo"]]))   # both box endpoints,
stopifnot(identical(a1[["clip.hi"]], a2[["clip.hi"]]))   #   as the bound requires

# The unsound call must also be refused rather than silently defaulted.
stopifnot(inherits(tryCatch(subsample.aggregate.beta.loglink(d1, "response",
    c("x1", "x2"), "gamma", "log", eps = 10, delta = 1e-6, m = 10, seed = 7L),
    error = function(e) e), "error"))
ok("Task 1: noise calibration is invariant between neighbouring datasets")

# +------------------------------------------------------------------------+
# |  Task 2 -- the release must not carry the empirical residual law        |
# +------------------------------------------------------------------------+
# Catches a release that carries the empirical residual law.
#
# The subsample-and-aggregate proposition guarantees "any release depending on
# the data only through beta_DP".  inject.virtual.y set y.virt = y + mu_DP -
# mu_hat, which retains the real y, so the real residuals reached the release
# and the hypothesis was not met.  A DP sigma^2 channel alone does not fix this:
# the three datasets below have residual sd exactly 2, so no dispersion channel
# can tell them apart, yet their shapes differ (Gaussian, heavy-tailed, and
# skewed).
#
# Hold X, beta and the residual SD fixed.  Vary only the residual shape.  The
# release's residual skewness and kurtosis must not track the input's.
set.seed(31L)
n2   <- 800L
x2a  <- rnorm(n2); x2b <- rnorm(n2)
eta2 <- 2 + 1.5 * x2a - 0.8 * x2b
std  <- function(v) 2 * (v - mean(v)) / sd(v)          # residual sd == 2, exactly

shape.of <- function(v.e) {
    d <- data.frame(response = eta2 + v.e, x1 = x2a, x2 = x2b)
    res <- gen.syn.dp.projected.subagg(d, y.nm = "response", family = "gaussian",
             link = "identity", x.nms = c("x1", "x2"), eps = 10, delta = 1e-6,
             m = 20, clip.lo = V.LO, clip.hi = V.HI,
             dispersion = dp.dispersion.public(4,
                 source = "design constant, residual variance 4"),
             seed = 99L)
    r <- residuals(lm(response ~ x1 + x2, res[["d.syn"]]))
    c(skew = skew(r), kurt = kurt(r))
}
m2 <- rbind(gauss = shape.of(std(rnorm(n2))),
            heavy = shape.of(std(rt(n2, df = 3))),
            skewd = shape.of(std(rexp(n2, rate = 1))))

# The three inputs have skew {~0, ~0, ~2} and kurt {~0, ~7, ~5}.  If the release
# echoed them, the rows below would differ by that much.  They must not: every
# row is a draw from the same public family at the same public dispersion.
stopifnot(max(abs(m2[, "skew"] - mean(m2[, "skew"]))) < 0.4)
stopifnot(max(abs(m2[, "kurt"] - mean(m2[, "kurt"]))) < 1.0)
ok("Task 2: release does not echo the empirical residual distribution")

# +------------------------------------------------------------------------+
# |  Task 3 -- a "public" constant must not be a sample statistic           |
# +------------------------------------------------------------------------+
# Catches a sample statistic presented as a public constant, which Task 1
# cannot catch.
#
# Note carefully why Task 1 misses it.  If someone reads the coefficients off
# Supplement S1 and hardcodes them, the box is a constant in the source, so it is
# trivially invariant between neighbours and Task 1 passes.  But S1's column is
# headed "Original" under the note "Seed 6489, n = 2000": it is coef(lm(...)) on
# the analysis sample.  The provenance reads "Supplement S1", which sounds
# impeccable, and the box is still a function of the sensitive data.
#
# The tell is statistical.  A population constant is a functional of the DGP: it
# is reproducible from any public seed and does not coincide with the fit on any
# particular sample.  A sample statistic matches its own sample exactly and no
# other.  So: the shipped design constants must be (a) stable across public seeds
# and (b) not equal to the analysis-sample fit.
if(file.exists("results/dp-design-constants.rds") &&
   file.exists("results/d.sim.scenarios.RData")) {
    l.k <- readRDS("results/dp-design-constants.rds")
    load("results/d.sim.scenarios.RData")
    d3 <- l.scenarios[[1L]][["data"]]
    d3 <- `[[<-`(d3, "female", value = factor(as.character(d3[["female"]]),
                                              levels = c("0", "1")))
    d3 <- `[[<-`(d3, "race", value = factor(as.character(d3[["race"]]),
                                            levels = c("1", "2", "3")))
    v.samp <- coef(lm(response ~ age + female + race + income, d3))
    v.pop  <- l.k[["sc1"]][["beta"]]

    # (a) they must not be the same object.  Equality here means someone sourced
    #     the "public" box from the sample, which is the S1 trap.
    stopifnot(!isTRUE(all.equal(unname(v.pop), unname(v.samp), tolerance = 1e-8)))

    # (b) but they must agree to within sampling error -- if they disagreed by
    #     more than the 50% asserted below, the "population" value would be
    #     measuring something else and the box would be wrong rather than
    #     merely leaky.
    s.rel <- max(abs(v.pop - v.samp) / abs(v.pop))
    stopifnot(s.rel > 1e-4, s.rel < 0.5)

    # (c) the provenance field must exist and must name its source.
    stopifnot(nzchar(l.k[["sc1"]][["source"]]),
              grepl("population", l.k[["sc1"]][["source"]], fixed = TRUE))
    ok(sprintf("Task 3: design constants are population values (%.1f%% from the sample fit, not 0%%)",
               100 * s.rel))
} else {
    stop("Task 3 needs dp-design-constants.rds; run 64-dp-design-constants.r")
}

# +------------------------------------------------------------------------+
# |  Task 4 -- epsilon must affect the release                              |
# +------------------------------------------------------------------------+
# The general check, which catches faults nobody has anticipated.
#
# If a channel reaches the release without passing through the noised mechanism,
# then turning epsilon down does not protect it.  An unprotected channel
# therefore shows an attack whose success is flat in epsilon.  A retained
# residual channel behaves exactly that way, scoring AUC 1.000 at eps = 10 and
# at eps = 1 alike, since the residuals do not pass through the noise.
#
# So: run a distinguishing game on a neighbouring pair at a large and a small
# epsilon.  The advantage at the small epsilon must be materially smaller.  A
# flat curve means something is escaping the accounting.
set.seed(41L)
n4  <- 500L
x4a <- rnorm(n4 - 1L); x4b <- rnorm(n4 - 1L)
y4  <- 2 + 1.5 * x4a - 0.8 * x4b + rnorm(n4 - 1L, 0, 2)
# The two worlds differ in one row: the probe row sits on, or far off, the line.
w4 <- function(off) data.frame(response = c(y4, 2 + off), x1 = c(x4a, 0),
                               x2 = c(x4b, 0))
d4.0 <- w4(0); d4.1 <- w4(20 * 2)

adv.at <- function(s.eps) {
    stat <- function(dd, s) {
        res <- tryCatch(gen.syn.dp.projected.subagg(dd, y.nm = "response",
                 family = "gaussian", link = "identity", x.nms = c("x1", "x2"),
                 eps = s.eps, delta = 1e-6, m = 20,
                 clip.lo = V.LO, clip.hi = V.HI,
                 dispersion = dp.dispersion.private(0, 40,
                     source = "design constant, sigma^2 in [0, 40]"),
                 seed = s), error = function(e) NULL)
        if(is.null(res)) return(NA_real_)
        sigma(lm(response ~ x1 + x2, res[["d.syn"]]))
    }
    v0 <- vapply(seq_len(40L), function(b) stat(d4.0, 1000L + b), numeric(1))
    v1 <- vapply(seq_len(40L), function(b) stat(d4.1, 5000L + b), numeric(1))
    v0 <- v0[is.finite(v0)]; v1 <- v1[is.finite(v1)]
    s.r <- rank(c(v0, v1))
    s.auc <- (sum(s.r[seq_along(v1) + length(v0)]) -
                  length(v1) * (length(v1) + 1) / 2) / (length(v0) * length(v1))
    abs(s.auc - 0.5) * 2                                  # advantage over chance
}
s.loose <- adv.at(10)
s.tight <- adv.at(0.5)

# The mu-GDP advantage bound falls by ~5x between eps 10 and eps 0.5, so a sound
# mechanism's empirical advantage must fall too.  A flat curve is the tell.
stopifnot(s.tight < s.loose - 0.15 || s.loose < 0.2)
ok(sprintf("Task 4: epsilon bites (advantage %.2f at eps 10 -> %.2f at eps 0.5)",
           s.loose, s.tight))

# +------------------------------------------------------------------------+
# |  Task 5 -- every channel'S claimed sensitivity must be its TRUE one     |
# +------------------------------------------------------------------------+
# Catches a bound taken over an object other than the one released, which
# Tasks 1-4 all miss.
#
# Note why they miss it.  The V4/V5 marginal channel's bound, sqrt(2), is a
# genuine public constant, correctly derived, and hardcoded in the source.  It is
# therefore invariant between neighbours (Task 1 passes), it is not a sample
# statistic (Task 3 passes), and it noises a channel that did not touch the
# response (Tasks 2 and 4 pass).  Every previous test is satisfied.  The bound is
# the sensitivity of the wrong function.
#
#   dp.hist.continuous built its counts with table() on the bin indices, and
#   table() drops the bins no row occupies.  The counts of the occupied bins were
#   left-packed and the zeros padded onto the end, so the released vector was not
#   the histogram: with an empty interior bin, every count above the gap slid down
#   one slot.  sqrt(2) is a histogram's L2 sensitivity (a replace-one swap
#   decreases one bin's count and increases another's).  The left-packed vector
#   has no such bound.  One row moving into an empty interior bin drags every
#   count above it, so its L2 sensitivity is O(n), not O(1).
#
# The test therefore cannot be "is the constant public".  It has to be: Take the shipped
# Mechanism, hand it the worst neighbouring pair you can build, and measure how
# far its released vector moves.  Re-deriving the count vector inside the test is
# how this survived: the arithmetic was checked against the intent, and
# not against the code.
V5.BREAKS <- seq(0, 5, length.out = 6L)                 # B = 5 public bins
S5.MU     <- 1e8                    # noise sd sqrt(2)/mu ~ 1e-8, negligible vs counts
S5.K      <- 100L
s.tol     <- 1e-6

# Recover the released count vector through the shipped mechanism.  probs is
# counts / sum(counts), and at this mu the noise is eight orders below the counts,
# so probs * n is what the mechanism actually binned.
counts.of <- function(v.x)
    dp.hist.continuous(v.x, V5.BREAKS, mu = S5.MU)[["probs"]] * length(v.x)

# Neighbouring pairs, one row apart, chosen the way an adversary would: to create
# or destroy an empty interior bin, which is what left-packing mishandles.
l.pairs <- list(
    list(nm = "open an interior gap",
         a = c(rep(0.5, S5.K), rep(2.5, S5.K)),
         b = c(rep(0.5, S5.K), rep(2.5, S5.K - 1L), 1.5)),
    list(nm = "close an interior gap",
         a = c(rep(0.5, S5.K), rep(2.5, S5.K - 1L), 1.5),
         b = c(rep(0.5, S5.K), rep(2.5, S5.K))),
    list(nm = "row crosses two gaps",
         a = c(rep(0.5, S5.K), rep(4.5, S5.K)),
         b = c(rep(0.5, S5.K), rep(4.5, S5.K - 1L), 2.5)),
    # Control: every bin occupied, so no packing can happen.  Ordinary data
    # does not exercise this path, so a test on realistic data alone
    # would find nothing.
    list(nm = "dense, no gap (control)",
         a = c(rep(0.5, S5.K), rep(1.5, S5.K), rep(2.5, S5.K), rep(3.5, S5.K),
               rep(4.5, S5.K)),
         b = c(rep(0.5, S5.K), rep(1.5, S5.K), rep(2.5, S5.K), rep(3.5, S5.K),
               rep(4.5, S5.K - 1L), 1.5)))

v.sens <- vapply(l.pairs, function(l.p)
    sqrt(sum((counts.of(l.p[["a"]]) - counts.of(l.p[["b"]]))^2)), numeric(1))
v.sens <- `names<-`(v.sens, vapply(l.pairs, `[[`, character(1), "nm"))

# sqrt(2) is the bound dp.hist.continuous calibrates to (noise sd sqrt(2)/mu).
# Anything above it is under-noising by exactly that ratio.
if(any(v.sens > sqrt(2) + s.tol))
    stop(sprintf(paste0("histogram channel under-noised by up to %.0fx: ",
                        "worst pair \"%s\" moves the released vector by %.1f, ",
                        "against the sqrt(2) = 1.414 it calibrates to."),
                 max(v.sens) / sqrt(2), names(which.max(v.sens)), max(v.sens)))

# The categorical channel makes the same sqrt(2) claim.  It builds its counts with
# factor(levels = ), so it was always right -- assert it so it stays that way.
V5.LEV <- c("a", "b", "c", "d")
cat.counts.of <- function(v.x) dp.cat(v.x, V5.LEV, mu = S5.MU)[["probs"]] * length(v.x)
s.cat <- sqrt(sum((cat.counts.of(c(rep("a", S5.K), rep("b", S5.K), rep("d", S5.K))) -
                   cat.counts.of(c(rep("a", S5.K), rep("b", S5.K),
                                   rep("d", S5.K - 1L), "c")))^2))
stopifnot(s.cat <= sqrt(2) + s.tol)

# The correlation channel.  What the mechanism releases is vech(S): it draws
# p(p+1)/2 iid normals, lays them on upper.tri(m.S, diag = TRUE) and mirrors
# (35-dp-full.r:222-226).  The mirroring is post-processing.  A Gaussian mechanism
# is mu-GDP at mu = Delta / sigma with Delta the L2 sensitivity of the released
# Vector, so Delta must be the sensitivity of vech(S) and of nothing else.
#
# ||vech(D)||_2 <= ||D||_F always, so a Frobenius measurement is a sufficient check
# that can only ever be over-satisfied.  It cannot see slack, and slack is what was
# there.  Worse, the check was `measured <= Delta`, which any large enough Delta
# passes.  A bound nothing attains is not a sensitivity.
#
# Two changes, and the second governs the result:
#   1. measure the released vector, vech(S), not the full matrix;
#   2. assert the bound is attained by the worst pair, not merely respected.
#
# The score map is the covariate's own released DP histogram, which is public by
# the time this channel runs.  Its sensitivity is therefore taken conditional on
# that map: one fixed map applied to both neighbours, which is how the code
# composes the two channels.
s.tau5 <- 0.005
s.M5   <- qnorm(1 - s.tau5)
s.n5   <- 400L
s.p5   <- 2L
set.seed(53L)
v.bg.u <- runif(s.n5 - 1L, 0.2, 4.8)
v.bg.v <- runif(s.n5 - 1L, 0.2, 4.8)
l.map5 <- Map(function(v.bg) dp.hist.continuous(v.bg, V5.BREAKS, mu = S5.MU),
              list(v.bg.u, v.bg.v))

# Both neighbours share every background row and differ only in the probe row, which
# is the replace-one adjacency both propositions use.
scores5 <- function(v.last) mapply(function(v.bg, s.x, l.h)
    dp.copula.scores(c(v.bg, s.x), V5.BREAKS[[1L]], V5.BREAKS[[6L]], s.tau5,
                     method = "adaptive", dp.hist = l.h),
    list(v.bg.u, v.bg.v), v.last, l.map5)

# The L2 norm of the half-vectorisation: the p(p+1)/2 entries the mechanism
# actually noises.  This is the only norm entitled to set the noise scale.
vech.l2 <- function(m.D) sqrt(sum(m.D[upper.tri(m.D, diag = TRUE)]^2))

move.of <- function(v.a, v.b) {
    m.A <- scores5(v.a)
    m.B <- scores5(v.b)
    # the clamp the bound rests on, checked on the shipped score map
    stopifnot(max(abs(m.A)) <= s.M5 + s.tol, max(abs(m.B)) <= s.M5 + s.tol)
    vech.l2(crossprod(m.A) / s.n5 - crossprod(m.B) / s.n5)
}

# Read Delta from the shipped mechanism, not re-derive it here.  A test that
# recomputes the constant it is guarding tests only its own arithmetic, and the two
# copies then drift.  dp.copula.corr calls this same function.
s.Delta5 <- dp.copula.corr.delta(s.n5, s.p5, s.M5)

# The maximiser, and it is realisable by data rather than by an adversary picking
# scores directly: a record outside the public box clamps to exactly +-M, so both
# probe rows sit on the clamp, and opposite signs in the second column make them
# Orthogonal.  Orthogonality eliminates the cross term and equal clamped lengths
# eliminate the diagonal of the move, which is exactly where the bound is attained.
s.sharp <- move.of(c(+1e6, +1e6), c(+1e6, -1e6))        # w = (M, M) against v = (M, -M)

# Soundness on hand-picked corners, including the zero-row case (a centre record
# scores ~0).
l.corr <- list(
    list(nm = "clamp corner, orthogonal (maximiser)", a = c(+1e6, +1e6), b = c(+1e6, -1e6)),
    list(nm = "centre row against clamp corner",      a = c(2.5, 2.5),   b = c(-1e6, +1e6)),
    list(nm = "clamp corner, sign-flipped both cols", a = c(+1e6, +1e6), b = c(-1e6, -1e6)),
    list(nm = "clamp corner, orthogonal (other way)", a = c(+1e6, +1e6), b = c(-1e6, +1e6)))
v.corr <- vapply(l.corr, function(l.p) move.of(l.p[["a"]], l.p[["b"]]), numeric(1))
v.corr <- `names<-`(v.corr, vapply(l.corr, `[[`, character(1), "nm"))

# Soundness under random search: No neighbouring pair may exceed Delta.  The draws
# span the public box and well outside it, so the clamp is exercised throughout.
set.seed(431L)
v.rand <- replicate(500L, move.of(runif(2L, -8, 13), runif(2L, -8, 13)))

if(any(c(v.corr, v.rand) > s.Delta5 + s.tol))
    stop(sprintf(paste0("correlation channel UNDER-noised: a neighbouring pair moves the ",
                        "released vector by %.6f against the Delta = %.6f it calibrates ",
                        "to. The mu-GDP guarantee is VOID."),
                 max(c(v.corr, v.rand)), s.Delta5))

# Sharpness.  Delta must be the sensitivity, not merely an upper bound on it.  This
# is the assertion a Frobenius measurement is structurally unable to make, and it
# fails on a constant that is only an upper bound, reporting it
# as a 2.00x
# over-noising.  Over-noising costs utility and not privacy, so this is not a
# soundness failure -- but a Delta that no dataset attains is a Delta nobody
# derived from the released object, and that is how the wrong object gets noised.
if(abs(s.sharp - s.Delta5) > s.tol)
    stop(sprintf(paste0("correlation channel's Delta is not its sensitivity: the worst ",
                        "neighbouring pair moves the released vector by %.6f, but the ",
                        "channel calibrates to Delta = %.6f (%.2fx). Re-derive Delta for ",
                        "vech(S), which is what the code releases."),
                 s.sharp, s.Delta5, s.Delta5 / s.sharp))

ok(sprintf(paste0("Task 5: channel sensitivities are the ones claimed (hist %.3f, cat %.3f ",
                  "vs sqrt(2) = 1.414; corr vech move %.6f == Delta %.6f, ATTAINED; ",
                  "%d random pairs, max %.6f)"),
           max(v.sens), s.cat, s.sharp, s.Delta5, length(v.rand), max(v.rand)))

# +------------------------------------------------------------------------+
# |  Task 6 -- the DP-ZI positive block must stay mu_P-GDP under a toggle   |
# +------------------------------------------------------------------------+
# Catches a partition that depends on the sensitive row count, which Tasks 1-5
# all miss for the same reason, namely that they reason about a fixed row set.
#
# gen.syn.dp.projected.zi privatises the positive component by
#   subsample.aggregate.beta.loglink(d[y > 0], ..., mu = mu_P),
# claiming mu_P-GDP under one-row replacement of the full data.  That release is
# (sqrt(p)/sigma.Q)-GDP if and only if the scaled aggregate Q = a/s has L2 sensitivity
# <= sqrt(p).  The subsample-and-aggregate proposition proves this for a
# partition of a fixed row set: one-row replacement perturbs exactly one block.
# The positive block, however, is fit on d[y > 0], and the mechanism partitioned it
# with sample(rep_len(1:m, n)) at n = |P| = sum(y > 0).  |P| is sensitive -- the
# zero-model channel exists to protect it -- so a replacement that toggles one
# row's positivity changes |P| by one, sample() then draws over a
# different-length vector, and the whole partition reshuffles.  Many blocks
# change, the "one block" premise fails, and the true sensitivity of a exceeds
# (hi_j - lo_j)/m: attains about 3x sqrt(p) below.
#
# The fix (28-dp-projection.r): partition the fixed full N rows and restrict that
# assignment to the positives, passed as `block`.  A toggle then changes exactly
# one block's positive subset, so a moves by at most s_j and the bound holds.
#
# The test hands the shipped channel adversarial neighbours built to clip: small
# |P|, and log-mean slopes beyond the +-8 box so each thin block clips to a corner
# by its composition.  It asserts:
#   (a) sound (block over the fixed N): every toggle pair moves Q by <= sqrt(p);
#   (b) discrimination (the partition over |P|): some pair attains > sqrt(p),
#       so (a) is a real constraint rather than a vacuous one.
#   (c) control (value-only, |P| fixed): <= sqrt(p), the harness self-check.
S6.M     <- 20L
S6.SQRTP <- sqrt(length(V.LO))                          # p = 3, the same public box
S6.SEED  <- 909L                                        # a public partition seed

# adversarial ZI data: n.pos corner-clipping positives + n.zero exact zeros.
make.zi6 <- function(n.pos, n.zero, seed) {
    set.seed(seed)
    typ <- rep(c(1, -1), length.out = n.pos)
    x1p <- runif(n.pos, -0.6, 0.6); x2p <- runif(n.pos, -0.6, 0.6)
    v.eta <- pmin(pmax(1 + typ * 12 * x1p + typ * 10 * x2p, -5), 6)
    rbind(data.frame(response = exp(v.eta), x1 = x1p, x2 = x2p),
          data.frame(response = rep(0, n.zero),
                     x1 = runif(n.zero, -0.6, 0.6),
                     x2 = runif(n.zero, -0.6, 0.6)))
}

# scaled aggregate Q = a / s from the shipped channel on the positive rows of dd,
# with the caller partition `blk` (NULL reproduces the |P| partition).
s6.Q <- function(dd, blk) tryCatch({
    d.pos <- dd[dd[["response"]] > 0, , drop = FALSE]
    l <- subsample.aggregate.beta.loglink(d.pos, "response", c("x1", "x2"),
             "quasi-gamma", "log", eps = 10, delta = 1e-6, m = S6.M,
             clip.lo = V.LO, clip.hi = V.HI, seed = 6489L, mu = 1, block = blk)
    a <- l[["beta.agg"]]
    a / l[["s.j"]][names(a)]
}, error = function(e) NULL)
s6.l2 <- function(u, v) if(is.null(u) || is.null(v)) NA_real_ else sqrt(sum((u - v)^2))

# one adversarial base -> (sound toggle move, |P| toggle move, control move).
s6.moves <- function(rep.i) {
    n.pos <- 90L + (rep.i %% 5L) * 10L
    dd <- make.zi6(n.pos, n.zero = 40L, seed = 3000L + rep.i)
    v.pos <- dd[["response"]] > 0
    set.seed(S6.SEED); v.blk <- sample(rep_len(seq_len(S6.M), nrow(dd)))
    k <- which.max(ifelse(v.pos, abs(dd[["x1"]]) + abs(dd[["x2"]]), -Inf))

    dd.t <- `[<-`(dd, k, "response", 0)                 # toggle probe row to zero
    dd.c <- `[<-`(dd, k, "response", dd[k, "response"] * 1.7 + 0.5)  # value only

    Q.base <- s6.Q(dd,   v.blk[v.pos])
    c(sound   = s6.l2(Q.base, s6.Q(dd.t, v.blk[dd.t[["response"]] > 0])),
      prefix  = s6.l2(s6.Q(dd, NULL), s6.Q(dd.t, NULL)),
      control = s6.l2(Q.base, s6.Q(dd.c, v.blk[dd.c[["response"]] > 0])))
}

m.s6 <- vapply(seq_len(12L), s6.moves, numeric(3))
s6.sound   <- max(m.s6["sound", ],   na.rm = TRUE)
s6.prefix  <- max(m.s6["prefix", ],  na.rm = TRUE)
s6.control <- max(m.s6["control", ], na.rm = TRUE)

stopifnot(s6.control <= S6.SQRTP + 1e-6)                # (c) harness sanity
if(s6.prefix <= S6.SQRTP)                               # (b) discrimination
    stop(sprintf(paste0("Task 6 is no longer discriminating: the |P| partition moved ",
                        "the aggregate by only %.3f, not past sqrt(p) = %.3f. ",
                        "Strengthen the adversarial data before trusting (a)."),
                 s6.prefix, S6.SQRTP))
if(s6.sound > S6.SQRTP + 1e-6)                          # (a) THE GUARD
    stop(sprintf(paste0("DP-ZI positive block UNDER-noised: a positivity toggle ",
                        "moves the scaled aggregate by %.4f against the sqrt(p) = ",
                        "%.4f it calibrates to (%.2fx). The mu_P-GDP guarantee is ",
                        "VOID. Partition the fixed full N, not |P|."),
                 s6.sound, S6.SQRTP, s6.sound / S6.SQRTP))

ok(sprintf(paste0("Task 6: DP-ZI positive block is toggle-sound (fixed-N ",
                  "partition max move %.4f <= sqrt(p) = %.4f, and the |P| ",
                  "partition attains %.4f = %.2fx; control %.4f)"),
           s6.sound, S6.SQRTP, s6.prefix, s6.prefix / S6.SQRTP, s6.control))

# +------------------------------------------------------------------------+
# |  Task 7 -- the mixed-margin copula's released vech(S) must stay bounded |
# +------------------------------------------------------------------------+
# Categorical covariates enter the correlation channel, and each is
# encoded by K-1 bounded normal-scored level indicators (cat.indicator.scores),
# appended to the continuous scores.  dp.copula.corr reads p = ncol(m.W), so the
# released vech(S) sensitivity is M^2 * p_total / n with p_total counting the
# indicator columns.  This hands the shipped score construction adversarial
# neighbours that drive the continuous rows onto the clamp and flip a
# categorical level, and asserts:
#   (a) sound: every such pair moves vech(S) by <= Delta(p_total).  The guard.
#   (b) discrimination: with a rare level, a level-flip pair moves it past Delta(p_c), the
#       continuous-only bound, so the indicator columns genuinely enter the
#       release and p_total, not p_c, is the width of the sensitivity.
#   (c) A3: the max move is strictly below Delta(p_total).  Unlike the all-
#       continuous channel (Task 5, attained), a categorical column takes one of
#       K fixed winsorised values strictly inside the clamp, so the mixed bound
#       is a sound upper bound that over-noises the categorical channel in the
#       safe direction.  The exact heterogeneous bound remains follow-on work.
S7.TAU <- 0.005; S7.M <- qnorm(1 - S7.TAU); S7.N <- 400L
S7.LEV <- c("a", "b", "c")                              # "c" is the rare level (discriminating)
set.seed(57L)
v7.u <- runif(S7.N - 1L, 0.2, 4.8); v7.v <- runif(S7.N - 1L, 0.2, 4.8)
v7.g <- c(rep("a", round((S7.N - 1L) * 0.60)), rep("b", round((S7.N - 1L) * 0.38)))
v7.g <- c(v7.g, rep("c", (S7.N - 1L) - length(v7.g)))   # ~2% prevalence of "c"
l7.hu  <- dp.hist.continuous(v7.u, V5.BREAKS, mu = S5.MU)
l7.hv  <- dp.hist.continuous(v7.v, V5.BREAKS, mu = S5.MU)
l7.cat <- dp.cat(v7.g, S7.LEV, mu = S5.MU)              # released proportions (public map)
scores7 <- function(s.u, s.v, s.g) cbind(
    dp.copula.scores(c(v7.u, s.u), V5.BREAKS[[1L]], V5.BREAKS[[6L]], S7.TAU, "adaptive", l7.hu),
    dp.copula.scores(c(v7.v, s.v), V5.BREAKS[[1L]], V5.BREAKS[[6L]], S7.TAU, "adaptive", l7.hv),
    cat.indicator.scores(factor(c(v7.g, s.g), levels = S7.LEV), l7.cat, S7.M))
move7 <- function(a, b) {
    m.A <- scores7(a[["u"]], a[["v"]], a[["g"]]); m.B <- scores7(b[["u"]], b[["v"]], b[["g"]])
    stopifnot(max(abs(m.A)) <= S7.M + s.tol, max(abs(m.B)) <= S7.M + s.tol)
    vech.l2(crossprod(m.A) / S7.N - crossprod(m.B) / S7.N)
}
s7.ptot <- ncol(scores7(1, 1, "a"))                     # 2 continuous + K-1 = 2 indicators = 4
s7.Dtot <- dp.copula.corr.delta(S7.N, s7.ptot, S7.M)    # read from the shipped constant
s7.Dpc  <- dp.copula.corr.delta(S7.N, 2L,      S7.M)    # the continuous-only bound
l7 <- list(
    list(nm = "clamp corners + flip to rare c", a = list(u=+1e6,v=+1e6,g="a"), b = list(u=+1e6,v=-1e6,g="c")),
    list(nm = "clamp corners + flip c to b",    a = list(u=+1e6,v=+1e6,g="c"), b = list(u=-1e6,v=+1e6,g="b")),
    list(nm = "flip only, continuous centred",  a = list(u=2.5, v=2.5, g="a"), b = list(u=2.5, v=2.5, g="c")))
v7c <- vapply(l7, function(p) move7(p[["a"]], p[["b"]]), numeric(1))
set.seed(571L)
v7r <- replicate(500L, { g2 <- sample(S7.LEV, 2L, replace = TRUE)
    move7(list(u=runif(1,-8,13), v=runif(1,-8,13), g=g2[1L]),
          list(u=runif(1,-8,13), v=runif(1,-8,13), g=g2[2L])) })
s7.max <- max(c(v7c, v7r))
if(s7.max > s7.Dtot + s.tol)                            # (a) THE GUARD
    stop(sprintf(paste0("mixed copula UNDER-noised: a neighbouring pair moves vech(S) by ",
                        "%.6f against Delta(p_total) = %.6f. The mu-GDP guarantee is VOID."),
                 s7.max, s7.Dtot))
if(max(v7c) <= s7.Dpc)                                  # (b) discrimination
    stop(sprintf(paste0("Task 7 is no longer discriminating: no level-flip pair cleared the continuous-only ",
                        "Delta(p_c) = %.6f (max %.6f), so p_total is not shown to bind. ",
                        "Sharpen the rare level."),
                 s7.Dpc, max(v7c)))
if(s7.max >= s7.Dtot - s.tol)                           # (c) A3: sound upper bound, not attained
    stop(sprintf(paste0("mixed Delta(p_total) = %.6f is being ATTAINED (max %.6f); A3 requires ",
                        "the categorical columns to sit off the clamp, so it must be a STRICT ",
                        "upper bound. If a mixed pair attains it, re-derive the sensitivity."),
                 s7.Dtot, s7.max))
ok(sprintf(paste0("Task 7: mixed copula respects Delta(p_total) = %.6f (max %.6f, %.2fx), clears ",
                  "Delta(p_c) = %.6f so p_total binds, and stays strictly below it per A3 ",
                  "(categorical channel over-noised, attained bound is follow-on work)"),
           s7.Dtot, s7.max, s7.max / s7.Dtot, s7.Dpc))

# +------------------------------------------------------------------------+
# |  Task 8 -- the suffstat release's three channels must be bounded, and   |
# |  two of them attained                                                   |
# +------------------------------------------------------------------------+
# Route 3 (66-dp-suffstat.r) releases vech(X_w'X_w), X_w'y_s, and
# y_s'y_s at the constants of dp.public.suffstat.sens().  This drives
# adversarial one-row replacements through the shipped map
# (dp.suffstat.stats) and asserts, per channel:
#   (a) sound: every realisable pair moves the release by <= its Delta.
#   (b) attained: X'y reaches 2 B_x B_y (same x at the max vertex, y flipped
#       across the clip) and y'y reaches B_y^2 (y from the centre to the
#       clip), so those two constants are exact sensitivities, not merely
#       bounds, as required.
#   (c) discrimination: the vech max exceeds B_x^2, so a constant that drops the
#       sqrt(2) (the Frobenius-vs-vech confusion in the other direction) is
#       rejected by a realisable pair.
#   (d) A3: the vech max stays strictly below sqrt(2) B_x^2.  The bound's
#       equality needs two one-hot max-norm whitened rows, which the box
#       does not contain, so the channel over-noises in the safe
#       direction.  The exact constant is recorded follow-on work.
# Plus (e): the noise calibration and the public-curve lambda must be
# Identical between neighbouring datasets (the Task 1 property, for the new
# mechanism).
source("lib/dp-public-bounds.r")                            # the constants Task 8 drives
source("lib/66-dp-suffstat.r")
l8.w <- dp.public.whitening(1L)
l8.y <- dp.public.y.std(1L)
l8.s <- dp.public.suffstat.sens(l8.w, l8.y)
m8.V <- .dp.box.vertices() %*% l8.w[["W"]]              # whitened box vertices
v8.ys <- c(-l8.y[["B.y"]], 0, l8.y[["B.y"]])            # standardised y extremes

# The row-level move formulas, tied to the shipped map below.
s8.vech.max <- 0; s8.xty.max <- 0; s8.yy.max <- 0
for(a in seq_len(nrow(m8.V))) for(b in seq_len(nrow(m8.V)))
    for(ya in v8.ys) for(yb in v8.ys) {
        u <- m8.V[b, ]; v <- m8.V[a, ]
        s8.vech.max <- max(s8.vech.max, vech.l2(tcrossprod(u) - tcrossprod(v)))
        s8.xty.max  <- max(s8.xty.max, sqrt(sum((u * yb - v * ya)^2)))
        s8.yy.max   <- max(s8.yy.max, abs(yb^2 - ya^2))
    }

# Tie the sweep to the shipped map: one replacement driven through
# dp.suffstat.stats must move the statistics by exactly the row formula.
d8 <- data.frame(age = runif(200L, 25, 80),
                 female = factor(sample(c("0", "1"), 200L, TRUE),
                                 levels = c("0", "1")),
                 race = factor(sample(c("1", "2", "3"), 200L, TRUE),
                               levels = c("1", "2", "3")),
                 income = runif(200L, 500, 3500),
                 response = runif(200L, 10, 90))
d8b <- d8
d8b[1L, c("age", "income", "response")] <- c(90, 4000, 140)
d8b[1L, "female"] <- "1"; d8b[1L, "race"] <- "3"
l8.a <- dp.suffstat.stats(d8, "response", l8.w, l8.y)
l8.b <- dp.suffstat.stats(d8b, "response", l8.w, l8.y)
u8 <- (model.matrix(~ age + female + race + income, dp.public.clip.x(d8b)) %*%
       l8.w[["W"]])[1L, ]
v8 <- (model.matrix(~ age + female + race + income, dp.public.clip.x(d8))  %*%
       l8.w[["W"]])[1L, ]
y8a <- (d8[1L, "response"]  - l8.y[["centre"]]) / l8.y[["scale"]]
y8b <- (pmin(d8b[1L, "response"], l8.y[["y.clip"]][[2L]]) -
        l8.y[["centre"]]) / l8.y[["scale"]]
if(abs(vech.l2(l8.b[["S"]] - l8.a[["S"]]) -
       vech.l2(tcrossprod(u8) - tcrossprod(v8))) > s.tol ||
   abs(sqrt(sum((l8.b[["b"]] - l8.a[["b"]])^2)) -
       sqrt(sum((u8 * y8b - v8 * y8a)^2))) > s.tol)
    stop("Task 8: the row-level move formula does not match the shipped map, ",
         "so the vertex sweep is not measuring the release the code makes")

if(s8.vech.max > l8.s[["Delta.vech"]] + s.tol ||
   s8.xty.max  > l8.s[["Delta.xty"]]  + s.tol ||
   s8.yy.max   > l8.s[["Delta.yy"]]   + s.tol)          # (a) THE GUARD
    stop(sprintf(paste0("suffstat release UNDER-noised: moves (%.4f, %.4f, %.4f) against ",
                        "Deltas (%.4f, %.4f, %.4f). The mu-GDP guarantee is VOID."),
                 s8.vech.max, s8.xty.max, s8.yy.max, l8.s[["Delta.vech"]],
                 l8.s[["Delta.xty"]], l8.s[["Delta.yy"]]))
if(s8.xty.max < l8.s[["Delta.xty"]] - s.tol ||
   s8.yy.max  < l8.s[["Delta.yy"]]  - s.tol)            # (b) ATTAINED
    stop(sprintf(paste0("Task 8 is no longer discriminating: X'y or y'y fails to attain its constant ",
                        "(%.4f of %.4f, %.4f of %.4f); a bound nothing attains was derived ",
                        "from a different object."),
                 s8.xty.max, l8.s[["Delta.xty"]], s8.yy.max, l8.s[["Delta.yy"]]))
if(s8.vech.max <= l8.w[["B.x"]]^2)                      # (c) discrimination
    stop(sprintf(paste0("Task 8 is no longer discriminating: the vech sweep max %.4f fails to clear ",
                        "B_x^2 = %.4f, so a sqrt(2)-dropping constant would not be rejected."),
                 s8.vech.max, l8.w[["B.x"]]^2))
if(s8.vech.max >= l8.s[["Delta.vech"]] - s.tol)         # (d) A3: strict
    stop(sprintf(paste0("the vech bound sqrt(2) B_x^2 = %.4f is being ATTAINED (%.4f); its ",
                        "equality needs one-hot max-norm rows the box excludes, so if a pair ",
                        "attains it, re-derive the sensitivity."),
                 l8.s[["Delta.vech"]], s8.vech.max))

# (e) the calibration is invariant between neighbours, and so is the lambda.
l8.ma <- gen.syn.dp.suffstats(d8, "response", sc = 1L, eps = 2, delta = 1e-6,
                              whitening = l8.w, y.std = l8.y, sens = l8.s,
                              seed = 8L)
l8.mb <- gen.syn.dp.suffstats(d8b, "response", sc = 1L, eps = 2, delta = 1e-6,
                              whitening = l8.w, y.std = l8.y, sens = l8.s,
                              seed = 8L)
if(!identical(l8.ma[["sigma.channels"]], l8.mb[["sigma.channels"]]) ||
   !identical(l8.ma[["lambda"]], l8.mb[["lambda"]]))
    stop("Task 8: the noise calibration or the public-curve lambda moved between ",
         "neighbouring datasets, so it is a function of the sample")

ok(sprintf(paste0("Task 8: suffstat channels bounded (vech %.4f of %.4f, %.2fx, strict per ",
                  "A3 and clears B_x^2 so the sqrt(2) binds), X'y and y'y attained, ",
                  "calibration and lambda neighbour-invariant"),
           s8.vech.max, l8.s[["Delta.vech"]], s8.vech.max / l8.s[["Delta.vech"]]))

cat("\nverify-dp-soundness: all tasks passed.\n")
