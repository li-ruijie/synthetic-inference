#!/usr/bin/env r
# +------------------------------------------------------------+
# |  Unified DP-Projection comparison: all 4 variants           |
# +------------------------------------------------------------+
#
# Single-run comparison of:
#   [V1] gen.syn.dp.projected     -- analysis-level DP (any X)
#   [V2] gen.syn.dp.ols.public    -- full DP, public X, OLS target
#   [V3] gen.syn.dp.ridge.public  -- full DP, public X, ridge(1) target
#   [V4] gen.syn.dp.full          -- full DP, sensitive X, OLS target
#
# Scenario 1 (gaussian identity, n=2000, p=6) at matched eps, matched seeds and
# matched B.  Output: results/mc-all-variants.csv.
#
# +------------------------------------------------------------+
# |  Where the bounds come from                                |
# +------------------------------------------------------------+
# A bound such as B.y <- 1.5 * max(abs(response)) is the sample maximum, so a
# "public" bound taken that way would be a statistic of the very data being
# protected.  B_y defines Route 1's sensitivity
# (Delta_beta = 2 B_y max_i ||(X'X)^-1 x_i||), and a noise scale that is a
# function of the private data instantiates no DP theorem at all.
#
# gen.syn.dp.full defaults to
# sens.method = "subagg", under which B.y and eps.split are dead (subagg takes its
# sensitivity from clipping the per-block coefficients into a public box, so it
# needs no bound on |y|) and clip.lo/clip.hi/disp.lo/disp.hi are required.  This
# Every bound comes from dp-public-bounds.r and from nowhere else.  The four
# variants compared are unchanged: this table is what the paper's DP-variant
# comparison rests on, so nothing was added to it and nothing was taken away.
#
# +------------------------------------------------------------+
# |  Which variants need A public bound on |y|, and why         |
# +------------------------------------------------------------+
# V2 and V3 (Route 1) need one as a theorem: OLS and ridge beta_hat have
# Unbounded global sensitivity when y is unbounded, since one row at |y| ->
# infinity moves beta_hat arbitrarily.  Scenario 1's bound comes from the DGP's
# own N(50, 15) marginal, clipped at mean +/- 6 sd (dp.public.B.y), which is
# public, costs no budget, and makes the bounded-response premise TRUE rather
# than assumed.
#
# V1 is calibrated to DFBETA local sensitivity, which ignores any bound on |y|.
# The driver therefore passes none, and passing one would leave the calibration
# unchanged while misrepresenting the guarantee this route provides.
#
# Route 1's estimand is therefore the clipped-OLS (and clipped-ridge) one, so V2
# and V3 are generated from, and scored against, the publicly clipped response.
# At 6 sd the clip binds with probability ~1e-9 per row, so on this sample it is
# the identity map and the clipped and unclipped baselines coincide numerically.
# The estimand it defines is still the clipped one, and the public-covariate
# exact release proposition must say so.
#
# On a scenario with no finite B_y (the gamma, Poisson and binomial families)
# Routes V1, V2 and V3 are genuinely unavailable.  They are then skipped and
# recorded as unavailable in the output, not run on an invented bound.  V4
# (Route 2, subagg) needs no bound on |y| and runs on any scenario.
#
# V1 is deprecated.  Its DFBETA local sensitivity is a function of
# the realised data, so it makes no formal (eps, delta)-DP claim, whatever bound
# it is handed.  It is kept in the comparison as the analysis-level baseline it
# has always been, and is otherwise untouched.
#
# +------------------------------------------------------------+
# |  V3 uses global sensitivities on both channels              |
# +------------------------------------------------------------+
# The sigma^2 channel uses gs.sigma2.ridge, a global sensitivity in the public X
# and the public B_y.  A beta channel calibrated to the ridge DFBETA would read
# the realised residuals, which is a local sensitivity and unsound.
#
# gen.syn.dp.ridge.public calibrates beta to
# gs.beta.ridge.public(m.X, lambda, B.y) = 2 B_y max_i ||A^-1 x_i||, with
# A = X'X + lambda I: a global sensitivity in the public X, the public lambda and
# the public B_y, reducing to gs.beta.ols.public exactly as lambda -> 0
# (34-dp-ridge-projection.r:55, used at :187).  Both of V3's channels are Gaussian
# and compose in one mu quadrature, so V3 now carries the mu-GDP it claims.
#
# The V3 column below therefore no longer reads small.  It sits beside V2 (9782%
# against 9983% at eps = 10), which is what theory predicts for the same route on
# the same X with the same B_y, differing only in the ridge shrinkage.  A large
# gap between them would be an error rather than a finding.
#
# +------------------------------------------------------------+
# |  kappa                                                     |
# +------------------------------------------------------------+
# V4's public box is +/- kappa * |beta_pop|, zero-centred.  kappa is swept
# over c(1.5, 3, 5) and carried into the output as a column, so utility is
# reported as A function of the box width: the correct answer to the obvious
# referee question of whether the public box smuggles in the answer.  V1, V2 and
# V3 carry no coefficient box, so kappa does not enter them and their rows record
# kappa = NA.  V4's mu.split is the generator default, c(marg = 0.05,
# corr = 0.05, coef = 0.90), the faithful translation of its earlier eps.split =
# c(X = 0.1, beta = 0.8, sigma = 0.1): under subagg the dispersion rides along in
# the coefficient aggregate as a (p + 1)-th coordinate, so beta + sigma = 0.9
# becomes coef = 0.90 and X = 0.1 becomes marg + corr = 0.10.
#
# Reading maxpct.  It is scale-dependent, and under the sound worst-case noise
# the income coefficient (beta ~ 0.013 for a predictor of order 1e3) dominates it
# on every Gaussian scenario.  The percentages are comparable across variants and
# epsilons, which is what this driver is for, and are not a claim about absolute
# utility.  62-mc-dp-gdp.r records the same caveat and reports decision agreement
# (sign, significance, coverage) alongside.
#
# +------------------------------------------------------------+

.sobol.main <- TRUE
.dp.public.no.main <- TRUE
.dp.ridge.no.main <- TRUE
.dp.full.no.main <- TRUE
source("lib/00-shared.r")
source("lib/gen-syn-family.r")
source("lib/07-calibrated-generator.r")
source("lib/09-symmetric-generator.r")
source("lib/11-projection-generator.r")
source("lib/28-dp-projection.r")
source("lib/33-dp-projection-public-x.r")
source("lib/34-dp-ridge-projection.r")
source("lib/35-dp-full.r")
source("lib/dp-public-bounds.r")            # the only permitted source of a bound
rm(.sobol.main, .dp.public.no.main, .dp.ridge.no.main, .dp.full.no.main)

load("results/d.sim.scenarios.RData")
S.SC <- 1L                              # this script is Scenario 1 only

d <- l.scenarios[[S.SC]][["data"]]
d <- `[[<-`(d, "female", value = factor(as.character(d[["female"]]),
                                        levels = c("0", "1")))
d <- `[[<-`(d, "race", value = factor(as.character(d[["race"]]),
                                      levels = c("1", "2", "3")))

s.B      <- 15L
v.eps    <- c(1, 2, 5, 10)
v.kappa  <- c(1.5, 3, 5)                # public box half-width, in |beta_pop|
s.delta  <- 1e-6
s.lambda <- 1                           # V3's ridge parameter
v.x.nms  <- c("age", "female", "race", "income")
fml      <- response ~ age + female + race + income
sup      <- c(response = "real", age = "real", income = "real")
x.spec   <- dp.public.xspec()           # PUBLIC breaks and levels

v.beta.o <- coef(lm(fml, d))            # V1's and V4's target: OLS on d

# +------------------------------------------------------------+
# |  The public response bound, and Route 1's clipping step     |
# +------------------------------------------------------------+
# dp.public.B.y errors where no finite bound on |y| exists, and that error is the
# point of the module. It says Route 1 does not exist for the scenario, not
# that the driver should go and find a bound elsewhere.  Catching it here is not a
# workaround.  It records V1, V2 and V3 as unavailable and skips them, which is
# the only sound treatment of a mechanism whose sensitivity is undefined.
l.By     <- tryCatch(dp.public.B.y(S.SC), error = function(e) NULL)
b.route1 <- !is.null(l.By)
d.r1     <- if(b.route1) dp.public.clip.y(d, "response", S.SC) else NULL
# y.clip is NULL where the DGP is structurally bounded (scenario 9), and the
# symmetric public range is then the bound itself.
# Route 1's own targets, on the response its bound is a bound for.
v.beta.r1  <- if(b.route1) coef(lm(fml, d.r1)) else NULL
v.ridge.r1 <- if(b.route1) fit.ridge(d.r1, "response", v.x.nms, s.lambda) else NULL

# +------------------------------------------------------------+
# |  One replication of one variant                            |
# +------------------------------------------------------------+
# Returns the vector of per-coefficient percentage deviations, or NULL if the
# replication failed.  Each variant is scored against its own target: V1 and V4
# against OLS on d, V2 against clipped OLS, V3 against clipped ridge(lambda).
one.run <- function(s.var, s.eps, s.kappa, s.seed) {
    res <- tryCatch(switch(s.var,
        V1 = gen.syn.dp.projected(d, y.nm = "response", family = "gaussian",
                link = "identity", supports = sup, eps = s.eps, delta = s.delta,
                seed = s.seed),
        V2 = gen.syn.dp.ols.public(d.r1, y.nm = "response", x.nms = v.x.nms,
                eps = s.eps, delta = s.delta, B.y = l.By[["B.y"]], seed = s.seed),
        V3 = gen.syn.dp.ridge.public(d.r1, y.nm = "response", x.nms = v.x.nms,
                lambda = s.lambda, eps = s.eps, delta = s.delta,
                B.y = l.By[["B.y"]], seed = s.seed),
        # Subagg: the public coefficient box and the public sigma^2 box.  B.y and
        # eps.split are dead here and are not passed.
        V4 = {
            l.box  <- dp.public.box(S.SC, s.kappa, names(v.beta.o))
            l.disp <- dp.public.disp.box(S.SC, s.kappa)
            gen.syn.dp.full(d, "response", x.spec, eps.total = s.eps,
                delta = s.delta,
                clip.lo = l.box[["lo"]],  clip.hi = l.box[["hi"]],
                disp.lo = l.disp[["lo"]], disp.hi = l.disp[["hi"]], seed = s.seed)
        }), error = function(e) NULL)
    if(is.null(res)) return(NULL)

    # fit.ridge returns an unnamed vector, so V3's comparison is positional.
    if(s.var == "V3") {
        v.b <- tryCatch(fit.ridge(res[["d.syn"]], "response", v.x.nms, s.lambda),
                        error = function(e) NULL)
        if(is.null(v.b)) return(NULL)
        return(abs(100 * (v.b - v.ridge.r1) / abs(v.ridge.r1)))
    }
    v.b <- tryCatch(coef(lm(fml, res[["d.syn"]])), error = function(e) NULL)
    if(is.null(v.b)) return(NULL)
    v.o      <- if(s.var == "V2") v.beta.r1 else v.beta.o
    v.common <- intersect(names(v.o), names(v.b))
    if(length(v.common) < 2L) return(NULL)
    abs(100 * (v.b[v.common] - v.o[v.common]) / abs(v.o[v.common]))
}

# One cell: B replications at a fixed (variant, eps, kappa), on the same seed
# block as every other cell, so the four variants are paired.  n.success is
# Reported, not assumed: silent replication loss is how this fleet came to write
# empty result files while every driver's tryCatch reported success.
cell <- function(s.var, s.kappa, s.eps) {
    l.p    <- Filter(Negate(is.null),
                     Map(function(s.seed) one.run(s.var, s.eps, s.kappa, s.seed),
                         6489L + seq_len(s.B)))
    s.n    <- length(l.p)
    v.max  <- vapply(l.p, max,  numeric(1))
    v.mean <- vapply(l.p, mean, numeric(1))
    data.frame(eps = s.eps, kappa = s.kappa, variant = s.var,
               maxpct.mean = if(s.n > 0L) mean(v.max)  else NA_real_,
               maxpct.sd   = if(s.n > 1L) sd(v.max)    else NA_real_,
               meanpct     = if(s.n > 0L) mean(v.mean) else NA_real_,
               n.success   = s.n, n.total = s.B, status = "ok")
}

# A scenario with no public B_y has no Route 1.  Record that, rather than
# inventing a bound or dropping the row.
cell.no.route1 <- function(s.var, s.eps)
    data.frame(eps = s.eps, kappa = NA_real_, variant = s.var,
               maxpct.mean = NA_real_, maxpct.sd = NA_real_, meanpct = NA_real_,
               n.success = 0L, n.total = 0L,
               status = "unavailable: no public bound on |y|")

# +------------------------------------------------------------+
# |  Provenance, printed                                       |
# +------------------------------------------------------------+
cat("UNIFIED DP-PROJECTION COMPARISON (Scenario 1, n=2000, p=6)\n")
cat(strrep("=", 80), "\n")
cat(sprintf("B = %d reps, delta = %g, kappa swept over %s, subagg m = %d (block-count rule)\n\n",
            s.B, s.delta, paste(v.kappa, collapse = ", "),
            dp.subagg.blocks(nrow(d), length(v.beta.o))))
dp.bounds.announce(
    coef = dp.public.box(S.SC, v.kappa[[2L]], names(v.beta.o)),
    disp = dp.public.disp.box(S.SC, v.kappa[[2L]]),
    B.y  = if(b.route1) l.By
           else list(source = "NONE EXISTS. V1, V2 and V3 are unavailable and are skipped."))

# +------------------------------------------------------------+
# |  Run                                                       |
# +------------------------------------------------------------+
# kappa enters V4 only, so V1, V2 and V3 are run once per eps rather than three
# times over identical numbers.
d.grid <- rbind(
    expand.grid(kappa = NA_real_, variant = c("V1", "V2", "V3"), eps = v.eps,
                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE),
    expand.grid(kappa = v.kappa, variant = "V4", eps = v.eps,
                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))

l.cells <- Map(function(s.var, s.kappa, s.eps) {
    if(!b.route1 && s.var %in% c("V1", "V2", "V3"))
        return(cell.no.route1(s.var, s.eps))
    cat(sprintf("  %s  eps %-3g %s\n", s.var, s.eps,
                if(is.na(s.kappa)) "" else sprintf(" kappa %-4g", s.kappa)))
    flush(stdout())
    cell(s.var, s.kappa, s.eps)
}, d.grid[["variant"]], d.grid[["kappa"]], d.grid[["eps"]])

d.out <- Reduce(rbind, l.cells)
d.out <- d.out[order(d.out[["eps"]], d.out[["variant"]], d.out[["kappa"]]), ]

# Refuse to write an all-zero result file.  A driver that exits 0 having dropped
# every replication is the failure this repair is for.
if(all(d.out[["n.success"]] == 0L))
    stop("36-: EVERY replication failed. The generators are erroring and the ",
         "tryCatch in one.run is swallowing it. Do not write an empty result file.")

# +------------------------------------------------------------+
# |  Report                                                    |
# +------------------------------------------------------------+
fmt <- function(s.var, s.eps, s.kappa) {
    v.sel <- d.out[["eps"]] == s.eps & d.out[["variant"]] == s.var &
             `if`(is.na(s.kappa), is.na(d.out[["kappa"]]),
                  !is.na(d.out[["kappa"]]) & d.out[["kappa"]] == s.kappa)
    r <- d.out[v.sel, ][1L, ]
    if(r[["status"]] != "ok")      sprintf("%11s", "unavail.")
    else if(r[["n.success"]] == 0) sprintf("%11s", "DROPPED")
    else sprintf("%9.1f%%", r[["maxpct.mean"]])
}

v.cols <- c(sprintf("%-11s", c("V1", "V2", "V3")),
            vapply(v.kappa, function(k) sprintf("%-11s", sprintf("V4 k=%g", k)),
                   character(1)))
cat("\n"); cat(strrep("=", 80), "\n")
cat("mean max|%| coefficient deviation\n\n")
cat(sprintf("%-6s %s\n", "eps", paste(v.cols, collapse = " ")))
cat(strrep("-", 80), "\n")
invisible(Map(function(s.eps) {
    v.f <- c(vapply(c("V1", "V2", "V3"),
                    function(v) sprintf("%-11s", fmt(v, s.eps, NA_real_)), character(1)),
             vapply(v.kappa,
                    function(k) sprintf("%-11s", fmt("V4", s.eps, k)), character(1)))
    cat(sprintf("%-6g %s\n", s.eps, paste(v.f, collapse = " ")))
}, v.eps))

# output is written to the repository root
write.csv(d.out, "results/mc-all-variants.csv", row.names = FALSE)

cat("\n"); print(d.out, row.names = FALSE)
cat("\nResults written to results/mc-all-variants.csv\n")
cat(strrep("=", 80), "\n")
cat("Legend:\n")
cat("  V1 = gen.syn.dp.projected      (analysis-level, deprecated, no formal DP)\n")
cat("  V2 = gen.syn.dp.ols.public     (Route 1, full DP, public X, OLS)\n")
cat("  V3 = gen.syn.dp.ridge.public   (Route 1, full DP, public X, ridge lambda=1; BOTH\n")
cat("                                  channels are global-sensitivity, so V3\n")
cat("                                  carries the mu-GDP it claims. See header.)\n")
cat("  V4 = gen.syn.dp.full           (Route 2, full DP, sensitive X, subagg OLS)\n")
