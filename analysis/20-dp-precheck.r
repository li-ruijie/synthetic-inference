#!/usr/bin/env Rscript
# +------------------------------------------------------------------+
# | 20-dp-precheck.r -- DP analytical pre-check (Gaussian on X'y)     |
# +------------------------------------------------------------------+
#
# What this supports.  Supplementary Material S6, "Gaussian mechanism on the
# raw sufficient statistic": the paper's argument that output perturbation on
# X'y is not a viable route to private regression, which is the reason the
# coefficient-space mechanisms exist at all.
#
# The mechanism.  Release X'y + eta, eta ~ N(0, sigma^2 I_p), and refit
#
#     beta_DP = (X'X)^-1 (X'y + eta) = beta_hat + (X'X)^-1 eta.
#
# Unbiased, but the noise is amplified through (X'X)^-1: a coefficient on a rare
# categorical level sits on a small eigenvalue of X'X, so its entry of (X'X)^-1
# is large and its deviation is orders of magnitude worse than a well-represented
# predictor's.  That amplification is the finding, and it is why the paper puts
# the noise on the coefficient instead.
#
# +------------------------------------------------------------------+
# | Every bound here comes from dp-public-bounds.r, and not        |
# | from the response.                                               |
# +------------------------------------------------------------------+
#
# Two constructions that must not be used:
#
#   Raw path:
#       l2.sensitivity <- function(m.X, v.y)
#           max(sqrt(rowSums(m.X^2))) * max(abs(v.y))     # <- the sample maximum
#
#   Normalised path:
#       v.y.n    <- (v.y - min(v.y)) / diff(range(v.y))   # <- the sample range
#       s.delta2 <- max(v.row.norms) * 1                  # "y in [0,1]", by construction
#
# Min-max scaling y by its own sample range makes B_y = 1 true for the observed
# sample by construction.  That is not a bounded universe.  It is the
# circularity dp-public-bounds.r exists to prevent, in its purest form.  Two
# neighbouring datasets get two different scalings, hence two different noise
# scales, hence no sensitivity bound holds and no theorem can be instantiated.
# The paper's remark on where the bounds come from says this.
#
# Both constructions also drop the factor of 2 for one-row replacement, the
# adjacency both propositions use.  Both errors shrink Delta_2, so the baseline
# would be noised less than a sound calibration allows.  A paper whose thesis is
# that data-derived bounds void the guarantee cannot calibrate its own baseline
# to a sample maximum.
#
# Every bound comes from dp-public-bounds.r and prints its `source`.
# Delta_2 = 2 B_x B_y under one-row replacement, which is what
# 42-baseline-comparison.r uses for the Monte Carlo of the same mechanism.  The
# two scripts must not disagree about the sensitivity of one mechanism.
#
# X is public here, and that is not A loophole.  The mechanism inverts X'X in the
# clear, so it cannot run unless X is public -- the same premise
# the public-covariate exact release proposition makes, so a bound read off X
# is sound and a bound read off y is unsound.  That asymmetry is the point.
#
# Where no public B_y exists, it says so rather than inventing one.
# Scenarios 3 (gamma) and 4 (zero-inflated quasi-gamma) have a response unbounded
# above, so X'y has no finite global sensitivity there and no epsilon calibrates
# the mechanism at all.  It is not merely too noisy, it is undefined.
# dp.public.B.y() already refuses those scenarios.  This script prints the refusal
# and lets it stand as the result.  That is a stronger statement than the "27%"
# and "67%" the supplement prints, and it is the same statement
# the remark on where the bounds come from makes for Route 1.  Route 2
# (subsample-and-aggregate) needs no bound on |y| whatever, so it is
# the paper's default.
#
# delta is 1e-6 throughout, the value the paper states everywhere else.  It used
# to be 1/n^2, which differs per dataset, so the PHI figure was reported at a
# different delta from the simulation figures and the two were not comparable.
#
# Output: dp-precheck.csv -- every number in S6.
# +------------------------------------------------------------------+

source("lib/00-shared.r")
source("lib/dp-gaussian-mechanism.r")  # dp.gm.multiplier, the one calibration
source("lib/dp-public-bounds.r")       # the ONLY permitted source of a DP bound

load("results/d.sim.scenarios.RData")

DELTA   <- 1e-6                        # the paper's delta, everywhere
V.EPS   <- c(0.1, 0.5, 1, 2, 5, 10)
N.MC    <- 10000L                      # draws for E[max |% deviation|]
MC.SEED <- 6489L

l.rows <- list()                       # accumulates the CSV

# +------------------------------------------------------------------+
# | Model fitting                                                     |
# +------------------------------------------------------------------+
# The two-part families are fitted on the positive subset, matching every other
# driver (see 64-dp-design-constants.r's fit.pop).
fit.model <- function(d, y.nm, family, link) {
    v.pred <- setdiff(colnames(d), y.nm)
    fml <- reformulate(v.pred, response = y.nm)
    d.fit <- `if`(link == "log" && any(d[[y.nm]] == 0), d[d[[y.nm]] > 0, ], d)
    mo <- `if`(family == "gaussian",
               lm(fml, d.fit),
               glm(fml, d.fit,
                   family = quasi(link = "log", variance = "mu^2"),
                   start = c(1, rep(0, ncol(model.matrix(fml, d.fit)) - 1L)),
                   control = glm.control(maxit = 200)))
    list(mo = mo, m.X = model.matrix(fml, d.fit), v.y = d.fit[[y.nm]])
}

# +------------------------------------------------------------------+
# | The public unit normalisation                                     |
# +------------------------------------------------------------------+
# Map each continuous covariate to [0, 1] by its public range and the response to
# [0, 1] by its public clip range.  Dummies are already in {0, 1}.  Nothing here
# reads a sample statistic: l.xrange and v.yclip are both supplied by the caller
# from dp-public-bounds.r.
#
# On that scale B_y = 1 by the public clip (not by construction on the sample),
# and B_x = sqrt(1 + n.cont + n.factor) in closed form.  This is the genuine best
# case: the smallest Delta_2 the mechanism can soundly be given.
unit.normalise <- function(d, y.nm, l.xrange, v.yclip) {
    for (nm in names(l.xrange)) {
        v.r <- l.xrange[[nm]]
        d[[nm]] <- (pmin(pmax(d[[nm]], v.r[[1L]]), v.r[[2L]]) - v.r[[1L]]) /
            (v.r[[2L]] - v.r[[1L]])
    }
    v.y <- pmin(pmax(d[[y.nm]], v.yclip[[1L]]), v.yclip[[2L]])
    d[[y.nm]] <- (v.y - v.yclip[[1L]]) / (v.yclip[[2L]] - v.yclip[[1L]])
    d
}

# +------------------------------------------------------------------+
# | The analysis                                                      |
# +------------------------------------------------------------------+
# Delta_2 = 2 B_x B_y, the L2 sensitivity of X'y under one-row replacement:
# replacing row i moves X'y by x_i' y_i' - x_i y_i, of norm at most
# ||x_i'|| |y_i'| + ||x_i|| |y_i| <= 2 B_x B_y.
analyse <- function(s.label, s.scale, d, y.nm, family, link, l.Bx, l.By) {
    l.f <- fit.model(d, y.nm, family, link)
    v.beta <- coef(l.f[["mo"]])
    m.X <- l.f[["m.X"]]
    s.p <- length(v.beta)
    s.n <- nrow(m.X)

    s.delta2 <- 2 * l.Bx[["B.x"]] * l.By[["B.y"]]

    m.XtX     <- crossprod(m.X)
    m.XtX.inv <- solve(m.XtX)
    v.var.diag <- diag(m.XtX.inv)
    v.eig <- sort(eigen(m.XtX, symmetric = TRUE, only.values = TRUE)[["values"]])

    cat(sprintf("\n%s -- %s  (n = %d, p = %d)\n", s.label, s.scale, s.n, s.p))
    dp.bounds.announce(B.x = l.Bx, B.y = l.By)
    cat(sprintf("  Delta_2 = 2 * B_x * B_y = 2 * %.4g * %.4g = %.4g   (delta = %.1e)\n",
                l.Bx[["B.x"]], l.By[["B.y"]], s.delta2, DELTA))
    cat(sprintf("  smallest eigenvalue of X'X = %.4g,  condition number = %.4g\n",
                min(v.eig), max(v.eig) / max(min(v.eig), 1e-10)))

    cat(sprintf("\n  %7s  %12s  %11s  %8s  %8s  %10s\n",
                "epsilon", "sigma_DP", "E|max%|", "P(>5%)", "P(>50%)", "ridge|max%|"))
    cat("  ", strrep("-", 68), "\n", sep = "")

    invisible(Map(function(s.eps) {
        s.sigma <- s.delta2 * dp.gm.multiplier(s.eps, DELTA)

        # sd(beta_DP_j - beta_hat_j) = sigma_DP * sqrt([(X'X)^-1]_jj)
        v.rel.sd <- s.sigma * sqrt(v.var.diag) / pmax(abs(v.beta), 1e-10)

        set.seed(MC.SEED)
        m.z <- matrix(rnorm(N.MC * s.p), nrow = N.MC)
        v.max.pct <- apply(abs(sweep(m.z, 2, v.rel.sd, `*`)), 1, max) * 100

        # Ridge shrinkage (the ridge-bias equation), for comparison: systematic
        # bias from drawing off a noisy conditional, distinct from the random
        # deviation.
        m.shrink <- solve(diag(s.p) + s.sigma^2 * m.XtX.inv)
        v.beta.r <- as.numeric(m.shrink %*% v.beta)
        s.max.r  <- 100 * max(abs(v.beta.r - v.beta) / pmax(abs(v.beta), 1e-10))

        cat(sprintf("  %7.1f  %12.4g  %10.1f%%  %7.1f%%  %7.1f%%  %9.1f%%\n",
                    s.eps, s.sigma, mean(v.max.pct),
                    100 * mean(v.max.pct > 5), 100 * mean(v.max.pct > 50), s.max.r))

        l.rows[[length(l.rows) + 1L]] <<- data.frame(
            dataset = s.label, scale = s.scale, available = TRUE,
            n = s.n, p = s.p, B.x = l.Bx[["B.x"]], B.y = l.By[["B.y"]],
            Delta2 = s.delta2, eps = s.eps, delta = DELTA, sigma.dp = s.sigma,
            e.max.pct = mean(v.max.pct),
            p.gt5 = mean(v.max.pct > 5), p.gt50 = mean(v.max.pct > 50),
            ridge.max.pct = s.max.r, stringsAsFactors = FALSE)
    }, V.EPS))
}

# The refusal.  Where the response is unbounded above, X'y has no finite global
# sensitivity, so there is no Delta_2 and no epsilon to sweep.  The module's own
# error message is the result: print it and let it stand.
refuse <- function(s.label, sc) {
    s.msg <- tryCatch({ dp.public.B.y(sc); NA_character_ },
                      error = function(e) conditionMessage(e))
    cat(sprintf("\n%s\n", s.label))
    cat("  NO PUBLIC B_y EXISTS, SO THE MECHANISM DOES NOT EXIST.\n")
    cat(sprintf("  %s\n", s.msg))
    cat("  Delta_2 = 2 B_x B_y is unbounded, so no (eps, delta) calibrates it at\n")
    cat("  any budget. This is not a bad utility trade-off, it is the absence of\n")
    cat("  a mechanism, and the supplement states this for Route 1.\n")

    l.rows[[length(l.rows) + 1L]] <<- data.frame(
        dataset = s.label, scale = NA_character_, available = FALSE,
        n = NA_integer_, p = NA_integer_, B.x = NA_real_, B.y = NA_real_,
        Delta2 = Inf, eps = NA_real_, delta = DELTA, sigma.dp = Inf,
        e.max.pct = NA_real_, p.gt5 = NA_real_, p.gt50 = NA_real_,
        ridge.max.pct = NA_real_, stringsAsFactors = FALSE)
}

# +------------------------------------------------------------------+
# | Run: simulation scenarios                                         |
# +------------------------------------------------------------------+
cat(strrep("=", 70), "\n")
cat("DP ANALYTICAL PRE-CHECK -- Gaussian mechanism on X'y\n")
cat(strrep("=", 70), "\n")
cat("Delta_2 = 2 B_x B_y (one-row replacement). Every bound is public and\n")
cat("announces its provenance below. No bound is a statistic of the sample.\n")

SIM.XRANGE <- list(age = c(20, 90), income = c(0, 4000))   # dp.public.xspec()
SIM.NCONT  <- 2L                                            # age, income
SIM.NFACT  <- 2L                                            # female, race

prep.sim <- function(sc) {
    d <- l.scenarios[[sc]][["data"]]
    d[["female"]] <- factor(as.character(d[["female"]]), levels = c("0", "1"))
    d[["race"]]   <- factor(as.character(d[["race"]]),   levels = c("1", "2", "3"))
    d
}

# --- Scenario 1: Gaussian identity. A public B_y exists via the 6-sd clip. ---
d1 <- dp.public.clip.x(dp.public.clip.y(prep.sim(1L), "response", 1L))
analyse("Scenario 1 gaussian(identity)", "raw", d1, "response",
        "gaussian", "identity", dp.public.B.x(), dp.public.B.y(1L))

analyse("Scenario 1 gaussian(identity)", "normalised",
        unit.normalise(prep.sim(1L), "response", SIM.XRANGE,
                       dp.public.B.y(1L)[["y.clip"]]),
        "response", "gaussian", "identity",
        dp.public.B.x.normalised(SIM.NCONT, SIM.NFACT),
        list(B.y = 1, source = paste0("1 by the PUBLIC clip of dp.public.B.y(1): ",
                                      "y in [-40, 140] rescaled to [0, 1]. Not by ",
                                      "min-max scaling on the sample.")))

# --- Scenarios 3 and 4: the response is unbounded above. No mechanism. ---
refuse("Scenario 3 gamma(log)", 3L)
refuse("Scenario 4 ZI quasi-gamma(log)", 4L)

# +------------------------------------------------------------------+
# | Run: PHI application                                              |
# +------------------------------------------------------------------+
# Waist circumference has a public bound -- a physiological range -- so unlike
# scenarios 3 and 4 the mechanism exists here, and it is catastrophic.
s.phi.path <- "phi_phase1_baseline.dta"
if (file.exists(s.phi.path)) {
    library(haven)
    d.raw <- read_dta(s.phi.path)
    d.app <- data.frame(
        waist   = as.numeric(d.raw[["waist_circum"]]),
        age     = 2016 - as.numeric(d.raw[["birthyear"]]),
        female  = factor(as.numeric(d.raw[["female"]])),
        race    = factor(as.numeric(d.raw[["race"]])),
        edu     = factor(findInterval(as.numeric(d.raw[["education"]]),
                                      c(5, 8)) + 1),
        housing = factor(as.numeric(d.raw[["housingtype"]])))
    d.app <- d.app[complete.cases(d.app), ]

    l.By.phi <- dp.public.B.y.phi()
    d.app[["waist"]] <- pmin(pmax(d.app[["waist"]], l.By.phi[["y.clip"]][[1L]]),
                             l.By.phi[["y.clip"]][[2L]])

    # PHI has no published support box, so B_x comes from the design matrix the
    # mechanism already treats as public (it inverts X'X in the clear).
    m.X.phi <- model.matrix(waist ~ age + female + race + edu + housing, d.app)
    l.Bx.phi <- dp.public.B.x.design(
        m.X.phi,
        source = paste0("max ||x||_2 over the REALISED design matrix, which this ",
                        "mechanism already treats as public: it inverts X'X in ",
                        "the clear, so it cannot run otherwise. A functional of ",
                        "public data. The response is the only sensitive object, ",
                        "and its bound is physiological."))

    analyse("PHI waist OLS", "raw", d.app, "waist", "gaussian", "identity",
            l.Bx.phi, l.By.phi)

    # age is the only continuous covariate.  Its range is public with X.
    PHI.XRANGE <- list(age = range(d.app[["age"]]))
    analyse("PHI waist OLS", "normalised",
            unit.normalise(d.app, "waist", PHI.XRANGE, l.By.phi[["y.clip"]]),
            "waist", "gaussian", "identity",
            dp.public.B.x.normalised(1L, 4L),   # age; female, race, edu, housing
            list(B.y = 1, source = paste0("1 by the PUBLIC physiological clip ",
                                          "[0, 250] cm rescaled to [0, 1].")))
} else {
    cat("\n\nPHI data not found at ", s.phi.path, " -- skipping.\n", sep = "")
}

# +------------------------------------------------------------------+
# | Emit                                                              |
# +------------------------------------------------------------------+
d.out <- do.call(rbind, l.rows)
write.csv(d.out, "results/dp-precheck.csv", row.names = FALSE)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("Written to dp-precheck.csv\n")
cat("Copy every reported number FROM THAT FILE. Not from this log, and not\n")
cat("by analogy with a sample-calibrated one.\n")
