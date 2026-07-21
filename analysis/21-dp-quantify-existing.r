#!/usr/bin/env Rscript
# +------------------------------------------------------------------+
# | Per-Person Influence via DFBETA                                   |
# +------------------------------------------------------------------+
#
# Supports the S6 run-in heading "Per-person influence (DFBETA)":
#   "the per-person disclosure can be quantified via the leave-one-out
#    influence [DFBETA; belsley-regressi-1980, cook-detectio-1977]"
#
# The projection generator releases exactly p+1 statistics from the
# original data (X'y and RSS). Its disclosure risk is equivalent to
# publishing regression output. We quantify this via:
#
# 1. DFBETA [cook-detectio-1977]: leave-one-out influence on each
#    coefficient.
#    dfbeta_ij = (X'X)^{-1} x_i e_i / (1 - h_ii)
#    max |dfbeta_ij| / |beta_j| = maximum per-person % influence
#    on any coefficient.
#
# 2. Cook's distance [cook-detectio-1977]: overall influence.
#    D_i = (e_i^2 h_ii) / (p * MSE * (1 - h_ii)^2)
#
# 3. Leverage: h_ii = x_i' (X'X)^{-1} x_i
#    High leverage means one person disproportionately determines
#    the fit. Singleton strata give h_ii -> 1.
#
# 4. Membership inference AUC (computed in 16-report-privacy.r).
#    AUC ~ 0.49 => no practical distinguishability.
#
# Key findings reported in the paper:
#   - Simulation (n=2000, p=6): max single-person influence is 1-7%.
#   - PHI application: singleton strata yield >1000% influence on
#     housing coefficients.
#   - "At least 3-5 colluding individuals would be needed to shift
#     any coefficient by 5%" (simulation scenarios).

# +------------------------------------------------------------------+
# | Setup                                                             |
# +------------------------------------------------------------------+
.shared.loaded <- TRUE
source("lib/00-shared.r")
load("results/d.sim.scenarios.RData")

# +------------------------------------------------------------------+
# | Core functions                                                    |
# +------------------------------------------------------------------+

## Fit model (same families as the simulation study).
fit.model <- function(d, y.nm, family, link) {
    v.pred <- setdiff(colnames(d), y.nm)
    fml <- reformulate(v.pred, response = y.nm)
    `if`(family == "gaussian",
         lm(fml, d),
         {
             d.fit <- `if`(link == "log" && any(d[[y.nm]] == 0),
                           d[d[[y.nm]] > 0, ], d)
             s.np <- ncol(model.matrix(fml, d.fit))
             `if`(family == "poisson",
                  glm(fml, d.fit, family = poisson()),
                  glm(fml, d.fit,
                      family = quasi(link = "log", variance = "mu^2"),
                      start = c(1, rep(0, s.np - 1)),
                      control = glm.control(maxit = 100)))
         })
}

## Quantify per-person influence on regression coefficients.
##
## For each person i and coefficient j, the DFBETA is:
##   dfbeta_ij = (X'X)^{-1} x_i e_i / (1 - h_ii)
##
## We express this as a percentage of |beta_j|:
##   inf_ij = |dfbeta_ij / beta_j| * 100
##
## The maximum over all (i, j) pairs is the worst-case single-person
## disclosure: how much one person can shift any published coefficient
## simply by being in or out of the dataset.
##
## We then compute the number of colluding individuals needed to shift
## a coefficient by 5%:
##   n(5%, optimistic) = ceil(5 / max_i inf_ij)   (aligned influences)
##   n(5%, typical)    = ceil(5 / median_i inf_ij) (random alignment)
quantify.influence <- function(s.label, d, y.nm, family, link) {
    mo <- fit.model(d, y.nm, family, link)
    v.beta <- coef(mo)
    s.p <- length(v.beta)

    v.pred <- setdiff(colnames(d), y.nm)
    fml <- reformulate(v.pred, response = y.nm)
    d.fit <- `if`(link == "log" && any(d[[y.nm]] == 0),
                  d[d[[y.nm]] > 0, ], d)
    s.n <- nrow(d.fit)

    # Compute influence measures: hat values, DFBETA, Cook's distance.
    l.inf <- influence(mo)
    m.dfbeta <- dfbetas(mo)  # standardised dfbetas
    m.dfb.raw <- `if`(family == "gaussian",
                      l.inf[["coefficients"]],  # lm: raw dfbeta in $coefficients
                      {
                          # GLM: compute raw dfbeta from hat values.
                          # dfbeta_i = (X'X)^{-1} x_i * w_i * r_i / (1 - h_ii)
                          # where w_i are working weights and r_i are working
                          # residuals from the IRLS fit.
                          m.X <- model.matrix(fml, d.fit)
                          v.e <- residuals(mo, type = "response")
                          v.h <- l.inf[["hat"]]
                          m.XtX.inv <- solve(crossprod(m.X))
                          # dfbeta_i = (X'X)^{-1} x_i * e_i / (1 - h_ii)
                          # for GLM with weights: use working residuals
                          v.wr <- residuals(mo, type = "working")
                          v.w <- weights(mo, type = "working")
                          v.we <- `if`(is.null(v.w) || length(v.w) == 0,
                                       v.wr, v.w * v.wr)
                          t(mapply(function(i) {
                              as.numeric(m.XtX.inv %*% m.X[i, ] *
                                         v.we[i] / (1 - v.h[i]))
                          }, seq_len(s.n)))
                      })

    v.h <- l.inf[["hat"]]
    v.cook <- cooks.distance(mo)

    # Per-person percentage influence on each coefficient:
    # |dfbeta_ij / beta_j| * 100
    # For simulation (n=2000, p=6) the paper reports 1-7%.
    # For PHI with singleton housing strata, this exceeds 1000%.
    m.pct.inf <- abs(sweep(m.dfb.raw, 2, v.beta, `/`)) * 100

    # For each coefficient: max influence by any single person
    v.max.inf <- apply(m.pct.inf, 2, max)

    # Number of colluding people needed to shift a coefficient by 5%.
    # Optimistic (influences perfectly aligned): 5 / max_single
    # Typical (random alignment): 5 / median_single
    v.med.inf <- apply(m.pct.inf, 2, median)
    v.n.for.5pct.opt <- ceiling(5 / v.max.inf)
    v.n.for.5pct.typ <- ceiling(5 / v.med.inf)

    cat(sprintf("\n%s  (n=%d, p=%d)\n", s.label, s.n, s.p))
    cat(sprintf("  Leverage:   max=%.4f  median=%.4f  (1/n=%.4f)\n",
                max(v.h), median(v.h), 1 / s.n))
    cat(sprintf("  Cook's D:   max=%.4f  median=%.6f\n",
                max(v.cook), median(v.cook)))
    cat(sprintf("  Most influential observation: i=%d (Cook=%.4f, h=%.4f)\n",
                which.max(v.cook), max(v.cook), v.h[which.max(v.cook)]))

    cat(sprintf("\n  Per-coefficient: max single-person influence\n"))
    cat(sprintf("  %-22s %10s %10s %10s %8s %8s\n",
                "coeff", "beta_hat", "max|dfb|", "max inf%",
                "n(5%,opt)", "n(5%,typ)"))
    cat("  ", strrep("-", 72), "\n", sep = "")

    Map(function(j) {
        cat(sprintf("  %-22s %10.4f %10.6f %9.3f%% %8d %8d\n",
                    names(v.beta)[j], v.beta[j],
                    max(abs(m.dfb.raw[, j])),
                    v.max.inf[j],
                    v.n.for.5pct.opt[j],
                    v.n.for.5pct.typ[j]))
    }, seq_along(v.beta))

    # Summary: the core claim for the paper.
    # Simulation: max influence 1-7%, so 3-5 colluders needed for 5%.
    # PHI: max influence >1000% for housing dummies (singleton strata),
    # meaning a single person can dominate a coefficient.
    cat(sprintf("\n  SUMMARY:\n"))
    # index into names(v.beta): v.max.inf is unnamed on the glm path, and
    # sprintf over the NULL name collapses to character(0), so cat printed
    # Nothing for those blocks (the per-scenario maxima trace to the
    # per-coefficient rows above)
    cat(sprintf("    Max single-person influence: %.3f%% (on %s)\n",
                max(v.max.inf), names(v.beta)[which.max(v.max.inf)]))
    cat(sprintf("    People needed to shift any coeff by 5%%: %d (optimistic) / %d (typical)\n",
                min(v.n.for.5pct.opt), min(v.n.for.5pct.typ)))
    cat(sprintf("    Equivalent: the projection reveals what %d colluding\n", min(v.n.for.5pct.opt)))
    cat(sprintf("    people could already infer from the published regression.\n"))

    invisible(list(
        max.inf = max(v.max.inf),
        n.5pct.opt = min(v.n.for.5pct.opt),
        n.5pct.typ = min(v.n.for.5pct.typ),
        max.cook = max(v.cook),
        max.leverage = max(v.h)))
}

# +------------------------------------------------------------------+
# | Run: simulation scenarios                                         |
# +------------------------------------------------------------------+
# For the simulation scenarios (n=2000, p=6), per-person influence is
# small (1-7%), confirming that the projection's inherent privacy is
# adequate when data are well-distributed.
cat("=" |> rep(64) |> paste(collapse = ""), "\n")
cat("INHERENT PRIVACY OF THE PROJECTION GENERATOR\n")
cat("=" |> rep(64) |> paste(collapse = ""), "\n")
cat("\nThe projection discloses p+1 aggregate statistics (= regression\n")
cat("output). We measure per-person influence via DFBETA.\n")
cat("  max inf% = max over all persons of |dfbeta_ij / beta_j|\n")
cat("  n(5%) = people needed to shift a coefficient by 5%%\n\n")

quantify.influence("Scenario 1: gaussian(identity)",
                   l.scenarios[[1]][["data"]], "response",
                   "gaussian", "identity")

quantify.influence("Scenario 3: gamma(log)",
                   l.scenarios[[3]][["data"]], "response",
                   "quasi-gamma", "log")

quantify.influence("Scenario 4: ZI quasi-gamma(log)",
                   l.scenarios[[4]][["data"]], "response",
                   "quasi-gamma", "log")

# +------------------------------------------------------------------+
# | Run: PHI application data                                         |
# +------------------------------------------------------------------+
# The PHI dataset reveals the vulnerability: housing has 9 levels,
# some forming singleton strata when crossed with other predictors.
# Singleton strata give h_ii -> 1 and inflate DFBETA to >1000% on
# housing coefficients. This motivates the pre-masking strategy
# explored in 22-premasking-demo.r.
s.phi.path <- "phi_phase1_baseline.dta"
if(file.exists(s.phi.path)) {
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
    quantify.influence("PHI waist OLS", d.app, "waist",
                       "gaussian", "identity")
} else {
    cat("PHI data not available (restricted). Skipping.\n")
}
