#!/usr/bin/env r
# +------------------------------------------------------------+
# |  Monte Carlo: mixed-margin DP copula, categorical fidelity  |
# +------------------------------------------------------------+
#
# The headline.  A companion instrument measures the continuous-block
# correlation fidelity of the DP Gaussian copula.  This one
# measures the categorical-to-continuous association the mixed-margin extension
# (35-dp-full.r, cat.indicator.scores + the probit/logit back-map) is meant to
# preserve.  Same Scenario 1 (gaussian identity, n = 2000), same public bounds,
# same descriptive-fidelity split mu.split = c(marg = 0.25, corr = 0.25,
# coef = 0.50) as 58-, so the two tables read side by side.
#
# Metric: the correlation ratio eta between each categorical (female binary, race
# nominal 3-level) and each continuous covariate (age, income).  eta^2 =
# SS_between / SS_total from the one-way decomposition of the continuous column
# across the categorical's levels, so eta lies in [0, 1] and is scale-free.  For
# the binary female it equals the absolute point-biserial correlation.  It is the
# categorical analogue of the continuous 0.606 -> 0.604 correlation recovery, so
# the story is the same: recovered toward the real value under the mixed copula,
# collapsed to the ~0 noise floor under independent margins.
#
# Three methods, not two:
#   mixed-probit  use.copula = TRUE,  back.map = "probit" (the recommended default)
#   mixed-logit   use.copula = TRUE,  back.map = "logit"  (substantiates A2: the
#                 softmax step injects independent choice noise, attenuating the
#                 recovered association)
#   independent   use.copula = FALSE (categoricals drawn from independent DP
#                 marginals.  The ~0 baseline)
#
# kappa is fixed at 3, not swept.  The recovered association lives entirely in the
# correlation channel (dp.copula.corr releases the extended vech(S)).  The public
# coefficient box +/- kappa |beta_pop| sets the noise on the separate subagg beta
# channel and touches none of the X columns.  eta is therefore exactly
# kappa-invariant, and 58- (section 17) already confirmed the continuous-block
# X-fidelity is kappa-invariant empirically for the very same channel.
#
# Controls carried alongside the headline:
#   frob.cont  Frobenius of cor(X_syn continuous) - cor(X continuous), i.e. the
#              58- metric, to confirm the continuous block still recovers once the
#              categorical columns join the copula.
#   beta.err   max|beta_syn - beta_dp|, the inferential invariant.  Must stay at
#              machine zero: the mixed copula is descriptive-only, beta_syn =
#              beta_dp exactly whatever joint structure X_syn carries.


# Single-thread BLAS: the n x n projection-matrix ops segfault under multi-worker
# PSOCK if each worker also multi-threads its BLAS.  This run is sequential, but
# the setting is kept for parity with 58- and costs nothing.
Sys.setenv(OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OMP_NUM_THREADS = "1")

.dp.full.no.main <- TRUE
source("lib/35-dp-full.r")
source("lib/dp-public-bounds.r")

load("results/d.sim.scenarios.RData")
d <- l.scenarios[[1]][["data"]]
d[["female"]] <- factor(as.character(d[["female"]]), levels = c("0", "1"))
d[["race"]]   <- factor(as.character(d[["race"]]),   levels = c("1", "2", "3"))

S.SC   <- 1L                                    # this script is Scenario 1 only
x.spec <- dp.public.xspec()                     # public breaks and levels

v.cont     <- c("age", "income")
m.cor.orig <- cor(d[v.cont])
v.beta.o   <- coef(lm(response ~ age + female + race + income, d))
l.box      <- dp.public.box(S.SC, 3, names(v.beta.o))
l.disp     <- dp.public.disp.box(S.SC, 3)

s.B        <- 100L
v.eps      <- c(1, 2, 5, 10)
s.kappa    <- 3                                 # FIXED (see header)
v.mu.split <- c(marg = 0.25, corr = 0.25, coef = 0.50)   # descriptive split (58-)
v.methods  <- c("mixed-probit", "mixed-logit", "independent")

# +------------------------------------------------------------+
# |  Association metrics                                        |
# +------------------------------------------------------------+
# Correlation ratio eta of a continuous v.num across the levels of a categorical
# v.cat.  Empty levels contribute 0 to the between-group sum of squares and are
# filtered so an unsampled synthetic level does not turn the estimate into NaN.
eta.ratio <- function(v.cat, v.num) {
    v.ok  <- !is.na(v.cat) & !is.na(v.num)
    v.cat <- v.cat[v.ok]; v.num <- v.num[v.ok]
    if(length(v.num) < 2L) return(NA_real_)
    s.gm  <- mean(v.num)
    l.g   <- Filter(function(g) length(g) > 0L, split(v.num, v.cat))
    ss.b  <- sum(vapply(l.g, function(g) length(g) * (mean(g) - s.gm)^2, numeric(1)))
    ss.t  <- sum((v.num - s.gm)^2)
    if(ss.t <= 0) return(NA_real_)
    sqrt(ss.b / ss.t)
}

# Spread of the group-mean income across race levels (income units), the metric
# the validation reported.  Kept as a scale-dependent secondary alongside eta.
spread.gm <- function(v.cat, v.num) {
    v.m <- vapply(split(v.num, v.cat),
                  function(g) if(length(g) > 0L) mean(g) else NA_real_, numeric(1))
    diff(range(v.m, na.rm = TRUE))
}

# The full metric vector on one data frame (real or synthetic).
metrics <- function(dd) {
    c(eta.female.age     = eta.ratio(dd[["female"]], dd[["age"]]),
      eta.female.income  = eta.ratio(dd[["female"]], dd[["income"]]),
      eta.race.age       = eta.ratio(dd[["race"]],   dd[["age"]]),
      eta.race.income    = eta.ratio(dd[["race"]],   dd[["income"]]),
      spread.race.income = spread.gm(dd[["race"]],   dd[["income"]]))
}
v.metric.nms <- names(metrics(d))

# +------------------------------------------------------------+
# |  One release                                               |
# +------------------------------------------------------------+
# Generate at (eps, method) on seed s.seed, then measure the association metrics
# on X_syn plus the two controls.  Returns a named vector, or NULL if the release
# or the refit failed (not silently counted as success).
mc.one <- function(s.seed, s.eps, s.method) {
    b.copula <- s.method != "independent"
    s.back   <- if(s.method == "mixed-logit") "logit" else "probit"
    res <- tryCatch(
        gen.syn.dp.full(d, "response", x.spec, eps.total = s.eps, delta = 1e-6,
                        clip.lo = l.box[["lo"]],  clip.hi = l.box[["hi"]],
                        disp.lo = l.disp[["lo"]], disp.hi = l.disp[["hi"]],
                        mu.split = v.mu.split, use.copula = b.copula,
                        back.map = s.back, seed = s.seed),
        error = function(e) NULL)
    if(is.null(res)) return(NULL)
    d.syn <- res[["d.syn"]]
    b.syn <- tryCatch(coef(lm(response ~ age + female + race + income, d.syn)),
                      error = function(e) NULL)
    if(is.null(b.syn)) return(NULL)
    v.common <- intersect(names(v.beta.o), names(b.syn))
    c(metrics(d.syn),
      frob.cont = sqrt(sum((cor(d.syn[v.cont]) - m.cor.orig)^2)),
      beta.err  = max(abs(b.syn[v.common] - res[["beta.dp"]][v.common])))
}

# +------------------------------------------------------------+
# |  Provenance, printed                                       |
# +------------------------------------------------------------+
dp.bounds.announce(
    coef = l.box,
    disp = l.disp)
cat(sprintf("kappa FIXED at %g; mu.split = (marg %.2f, corr %.2f, coef %.2f)\n",
            s.kappa, v.mu.split[["marg"]], v.mu.split[["corr"]], v.mu.split[["coef"]]))
cat(sprintf("B = %d releases per (eps, method) cell, delta = 1e-6\n\n", s.B))

v.real <- metrics(d)
cat("REAL-DATA ASSOCIATION TARGETS (Scenario 1, n = 2000):\n")
invisible(Map(function(nm, x) cat(sprintf("  %-20s %.4f\n", nm, x)),
              v.metric.nms, v.real[v.metric.nms]))
cat("\n")

# +------------------------------------------------------------+
# |  Run                                                       |
# +------------------------------------------------------------+
# One independent seed per (eps, method, rep) row, drawn once from the master seed
# 6489, exactly as 58- draws its grid seeds.  Each release is fully determined by
# its own seed, so the run is reproducible whatever the evaluation order.
l.grid <- expand.grid(eps = v.eps, method = v.methods, rep = seq_len(s.B),
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
set.seed(6489L)
v.seeds <- sample.int(.Machine[["integer.max"]], nrow(l.grid))
s.n     <- nrow(l.grid)

l.raw <- Map(function(i, s.seed, s.eps, s.method) {
    if(i %% 50L == 0L) { cat(sprintf("  %d/%d\n", i, s.n)); flush(stdout()) }
    tryCatch(mc.one(s.seed, s.eps, s.method), error = function(e) NULL)
}, seq_len(s.n), v.seeds, l.grid[["eps"]], l.grid[["method"]])
cat("MAP DONE\n"); flush(stdout())

l.grid[["ok"]] <- vapply(l.raw, Negate(is.null), logical(1))

# +------------------------------------------------------------+
# |  Aggregate to cells                                        |
# +------------------------------------------------------------+
# For each (eps, method) cell, the mean and sd of every eta over the successful
# releases, the mean continuous-block Frobenius, the worst beta invariant error,
# and n.success (reported, not assumed).
v.eta.nms <- c("eta.female.age", "eta.female.income",
               "eta.race.age", "eta.race.income")

d.cells <- expand.grid(eps = v.eps, method = v.methods,
                       KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
cell <- function(s.eps, s.method) {
    v.sel  <- which(l.grid[["eps"]] == s.eps & l.grid[["method"]] == s.method &
                    l.grid[["ok"]])
    m.m    <- if(length(v.sel) > 0L)
                  do.call(rbind, l.raw[v.sel]) else matrix(numeric(0), 0L, 0L)
    col    <- function(nm, f) if(nrow(m.m) > 0L) f(m.m[, nm]) else NA_real_
    v.eta.mean <- vapply(v.eta.nms, function(nm) col(nm, mean), numeric(1))
    v.eta.sd   <- vapply(v.eta.nms, function(nm)
                         if(nrow(m.m) > 1L) sd(m.m[, nm]) else NA_real_, numeric(1))
    data.frame(scenario = S.SC, eps = s.eps, method = s.method,
               t(`names<-`(v.eta.mean, paste0(v.eta.nms, ".mean"))),
               t(`names<-`(v.eta.sd,   paste0(v.eta.nms, ".sd"))),
               spread.race.income = col("spread.race.income", mean),
               frob.cont          = col("frob.cont", mean),
               beta.err.max       = col("beta.err", max),
               n.success = length(v.sel), n.total = s.B,
               check.names = FALSE)
}

d.body <- Reduce(rbind, Map(cell, d.cells[["eps"]], d.cells[["method"]]))

# The real-data targets as a self-contained header block: method = "real",
# eps = NA, so the manuscript reads the target and the recovered value from one
# file.  The controls are NA for the real data (there is no synthetic release to
# compare against itself).
d.real <- data.frame(
    scenario = S.SC, eps = NA_real_, method = "real",
    t(`names<-`(v.real[v.eta.nms],   paste0(v.eta.nms, ".mean"))),
    t(`names<-`(rep(NA_real_, 4L),   paste0(v.eta.nms, ".sd"))),
    spread.race.income = v.real[["spread.race.income"]],
    frob.cont = NA_real_, beta.err.max = NA_real_,
    n.success = 1L, n.total = 1L, check.names = FALSE)

d.out <- rbind(d.real, d.body)

if(all(d.body[["n.success"]] == 0L))
    stop("58e-: EVERY release failed. The generator is erroring and the tryCatch ",
         "in mc.one is swallowing it. Do not write an empty result file.")

write.csv(d.out, "results/mc-mixed-copula.csv", row.names = FALSE)

# +------------------------------------------------------------+
# |  Report                                                    |
# +------------------------------------------------------------+
cat(strrep("=", 78), "\n")
cat("eta(race, income): real target vs recovered, by method and eps\n\n")
cat(sprintf("%-14s %-6s %-18s %-12s %-10s\n",
            "method", "eps", "eta.race.income", "frob.cont", "beta.err"))
cat(strrep("-", 78), "\n")
cat(sprintf("%-14s %-6s %-18.4f %-12s %-10s\n",
            "real", "-", v.real[["eta.race.income"]], "-", "-"))
invisible(Map(function(s.method, s.eps) {
    r <- d.body[d.body[["method"]] == s.method & d.body[["eps"]] == s.eps, ]
    cat(sprintf("%-14s %-6g %-18s %-12.4f %-10.1e\n",
                s.method, s.eps,
                sprintf("%.4f (%.4f)", r[["eta.race.income.mean"]],
                        r[["eta.race.income.sd"]]),
                r[["frob.cont"]], r[["beta.err.max"]]))
}, d.cells[["method"]], d.cells[["eps"]]))

cat("\n"); print(d.out, row.names = FALSE)
cat("\nWrote mc-mixed-copula.csv\n")
