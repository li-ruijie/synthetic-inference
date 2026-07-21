#!/usr/bin/env r
# +------------------------------------------------------------------------+
# |  64-dp-design-constants.r -- the public constants the DP boxes rest on  |
# +------------------------------------------------------------------------+
#
# Why this script exists.  Every DP mechanism here needs a public box: the box
# width sets s_j = (u_j - l_j)/m, which sets the noise scale, so a box derived
# from the sensitive sample makes the noise scale a function of the private data
# and no DP theorem can be instantiated (the sensitivity argument binds one box
# and applies it to both neighbours).  This project has now been
# bitten by that seven times.  This script is the single place the public
# constants come from, so a driver does not have to invent one.
#
# The trap this script exists to avoid.  It is tempting to read the "true
# coefficients" off Supplement S1.  Do not.  S1's column is headed "Original"
# and its note reads "Seed 6489, n = 2000": it is coef(lm(...)) on the analysis
# Sample, a statistic of the very data being protected.  Sourcing a box from it
# looks impeccable (the provenance reads "Supplement S1") and is exactly the
# error, carried through the paper.
#
# What is actually public.  00-simulate.r's DGP is a Gaussian copula: five
# correlated latent normals (fixed correlation matrix m.cor) pushed through fixed
# marginal transforms.  There is no beta parameter in it, but the population
# regression coefficient -- the L2 projection of Y on (1, X) under the true joint
# law -- is a deterministic functional of those published constants, and can be
# computed to arbitrary accuracy from an independent draw at a public seed.  That
# value depends on the DGP alone, is reproducible by any reader from the paper's
# own description of it, and touches the n = 2000 seed-6489 analysis sample
# nowhere.  That is the quantity a public box may rest on.
#
# The same applies to the population residual variance (for the sigma^2 box) and
# to the response bound B_y (for Route 1).
#
# Output: results/dp-design-constants.rds, results/dp-design-constants.log
#         Consumed by dp-public-bounds.r, which turns them into bounds objects.
#
# +------------------------------------------------------------------------+
.sobol.main <- TRUE
source("lib/00-shared.r")
rm(.sobol.main)

# +------------------------------------------------------------------------+
# |  The DGP, and only the DGP                                             |
# +------------------------------------------------------------------------+
# Lift sim.scenario() out of 00-simulate.r without executing its main block,
# which would regenerate and overwrite d.sim.scenarios.RData.  Note
# what is not loaded: d.sim.scenarios.RData itself.  This script must not see
# the analysis sample, and cannot, the object being out of scope.
v.exprs <- parse("simulate/00-simulate.r")
b.simdef <- Map(function(e) {
    is.call(e) && identical(as.character(e[[1L]]), "<-") &&
        identical(as.character(e[[2L]]), "sim.scenario")
}, v.exprs) |> unlist()
invisible(Map(function(e) eval(e, envir = globalenv()), v.exprs[b.simdef]))
stopifnot(is.function(sim.scenario))

N.POP    <- 1e6L        # population draw; MC error below is reported, not assumed
PUB.SEED <- 20260713L   # PUBLIC, and deliberately != 6489 (the analysis seed)
ANALYSIS.SEED <- 6489L  # named here only so the assertion below can forbid it
stopifnot(PUB.SEED != ANALYSIS.SEED)

fml <- response ~ age + female + race + income

# Coerce the two categoricals to the levels every driver uses.
prep <- function(d) {
    d <- `[[<-`(d, "female", value = factor(as.character(d[["female"]]),
                                            levels = c("0", "1")))
    `[[<-`(d, "race", value = factor(as.character(d[["race"]]),
                                     levels = c("1", "2", "3")))
}

# Which families are genuinely zero-inflated, hence two-part.  Key off the
# Family, not off "the response happens to contain zeros": scenario 6 is
# Poisson(8), whose zeros are the Poisson's own mass at 0 (P = e^-8, about 0.03%)
# and not structural inflation.  Testing for zeros instead would fit scenario 6's
# coefficient on the positive subset (a truncated Poisson, which no driver fits)
# and hang a near-separated logistic zero-model off it.  Only quasi-gamma
# (scenarios 4, 5, 8) is two-part, and only those reach gen.syn.dp.projected.zi.
is.zi <- function(family) family == "quasi-gamma"

# The analyst's model, per family.  Must match what the drivers fit: 62-'s
# fit.model puts Poisson and Gamma on the full data, and only the two-part
# quasi-gamma is fitted on the positive subset.
fit.pop <- function(d, family, link) {
    d.fit <- `if`(is.zi(family), d[d[["response"]] > 0, , drop = FALSE], d)
    switch(family,
        gaussian = lm(fml, d.fit),
        binomial = glm(fml, d.fit, family = binomial()),
        poisson  = glm(fml, d.fit, family = poisson()),
        gamma    = glm(fml, d.fit, family = Gamma(link = "log")),
        `quasi-gamma` = glm(fml, d.fit,
                            family = quasi(link = "log", variance = "mu^2"),
                            start = c(1, rep(0, 5)),
                            control = glm.control(maxit = 200)))
}

# The logistic zero-model, for the ZI generator's second box.
fit.pop.zero <- function(d) {
    d.z <- prep(d)
    d.z <- `[[<-`(d.z, ".Z", value = as.integer(d[["response"]] > 0))
    glm(.Z ~ age + female + race + income, d.z, family = binomial("logit"))
}

# +------------------------------------------------------------------------+
# |  Public response bounds (Route 1 only, subagg needs none)              |
# +------------------------------------------------------------------------+
# Route 1's sensitivity is defined through B_y (Delta_beta = 2 B_y max_i
# ||(X'X)^-1 x_i||), because OLS beta_hat has unbounded global sensitivity when y
# is unbounded.  So Route 1 needs a public bound, and 1.5 * max|y| was not one.
#
# Scenarios 1, 2 and 10 draw y = qnorm(u, 50, 15), which is genuinely unbounded,
# so no finite B_y exists and the mechanism needs an explicit public clipping
# Step.  We clip to the DGP's own marginal mean +/- 6 sd.  That is public (it
# reads the DGP, not the sample), it moves the estimand to the clipped-OLS one
# (which must be declared in the public-covariate exact release proposition),
# and at 6 sd it binds with probability ~1e-9 per row, so on any realised sample
# it is the identity map.
#
# Scenario 9 needs no clipping: its DGP is round(pmax(pmin(y, 100), 0)), so
# y in [0, 100] structurally.  B_y = 100 is exact, public, and spends no budget.
y.bounds.of <- function(sc) {
    if(sc %in% c(1L, 2L, 10L)) {
        v.clip <- c(50 - 6 * 15, 50 + 6 * 15)          # N(50, 15) marginal
        list(y.clip = v.clip, B.y = max(abs(v.clip)),
             source = "DGP marginal N(50, 15), clipped at mean +/- 6 sd")
    } else if(sc == 9L) {
        list(y.clip = NULL, B.y = 100,
             source = "DGP structural bound: y = round(pmax(pmin(., 100), 0))")
    } else list(y.clip = NULL, B.y = NA_real_,
                source = "not applicable (subagg needs no bound on |y|)")
}

# +------------------------------------------------------------------------+
# |  Compute                                                               |
# +------------------------------------------------------------------------+
cat("PUBLIC DESIGN CONSTANTS from the DGP (not from the analysis sample)\n")
cat(strrep("=", 78), "\n")
cat(sprintf("N = %g, public seed = %d (analysis seed %d is NOT used)\n\n",
            N.POP, PUB.SEED, ANALYSIS.SEED))

l.const <- Map(function(sc) {
    s <- sim.scenario(sc, n = N.POP, seed = PUB.SEED)
    d <- prep(s[["data"]])
    fam <- s[["family"]]; lnk <- s[["link"]]

    mo <- fit.pop(d, fam, lnk)
    v.beta <- coef(mo)

    # Population dispersion.  Gaussian: the residual variance.  Gamma and
    # quasi-gamma: the GLM dispersion phi (so shape = 1/phi), which is what
    # synth.response draws with.  Poisson and binomial have no free dispersion --
    # their variance is fixed by the mean -- so they need none and cannot leak one.
    s.sigma2 <- switch(fam,
        gaussian = sigma(mo)^2,
        gamma = ,
        `quasi-gamma` = summary(mo)[["dispersion"]],
        NA_real_)

    # The ZI generator privatises a logistic zero-model too, so it needs its own
    # public box.  Only the two-part (quasi-gamma) scenarios have one.
    b.zi <- is.zi(fam)
    mo.zero <- `if`(b.zi, fit.pop.zero(s[["data"]]), NULL)
    v.gamma <- `if`(b.zi, coef(mo.zero), NULL)

    # Per-fitted-observation coefficient variance, diag(vcov) x nobs, for the
    # two-stage box floor.  A block fit on b rows has sampling variance close to
    # beta.var1 / b, so the public block-spread scale at block count m is
    # tau_j(m) = sqrt(beta.var1_j * m / n).  Both read the vcov of the population
    # fits already computed above, so no RNG is touched and every pre-existing
    # field in this file stays byte-identical.  nobs
    # counts the observations in the fit, which for the two-part families is the
    # positive subset, and the positive block size is then recovered through
    # pos.frac below.
    v.var1 <- diag(vcov(mo)) * nobs(mo)
    v.var1.z <- `if`(b.zi, diag(vcov(mo.zero)) * nobs(mo.zero), NULL)

    # The public positive fraction, setting the ZI positive-channel block size
    # (n / m) * pos.frac.  A functional of the same public draw, again no RNG.
    s.pos <- `if`(b.zi, mean(d[["response"]] > 0), NULL)

    # The population second-moment matrix Sigma = E[x x'] over the public support
    # box, for the whitened sufficient-statistic release.  Whitening by
    # Sigma^{-1/2} is a fixed public affine map, so it costs no budget, and it
    # conditions X'X near the identity.  Computed from the same public draw by
    # crossprod alone, so the RNG stream above is untouched and every other
    # constant in this file is unchanged.  The clip to
    # the public box is the sensitivity domain's premise.  On these DGPs it does not
    # bind, so it is the identity map here and the provenance statement anyway.
    d.clip <- d
    for(nm in c("age", "income")) {
        v.r <- switch(nm, age = c(20, 90), income = c(0, 4000))
        d.clip[[nm]] <- pmin(pmax(d.clip[[nm]], v.r[[1L]]), v.r[[2L]])
    }
    m.Sigma <- crossprod(model.matrix(~ age + female + race + income, d.clip)) /
        nrow(d.clip)

    l.y <- y.bounds.of(sc)

    cat(sprintf("--- Scenario %2d  %s(%s)%s\n", sc, fam, lnk,
                `if`(b.zi, "  [zero-inflated]", "")))
    cat(sprintf("    beta   %s\n",
                paste(sprintf("%s=%.4f", names(v.beta), v.beta), collapse = "  ")))
    if(!is.na(s.sigma2))
        cat(sprintf("    %s %.4f\n",
                    `if`(fam == "gaussian", "sigma2", "disp  "), s.sigma2))
    if(b.zi) cat(sprintf("    gamma  %s\n",
                paste(sprintf("%s=%.4f", names(v.gamma), v.gamma), collapse = "  ")))
    cat(sprintf("    var1   %s\n",
                paste(sprintf("%s=%.4g", names(v.var1), v.var1), collapse = "  ")))
    if(b.zi) {
        cat(sprintf("    zvar1  %s\n",
                    paste(sprintf("%s=%.4g", names(v.var1.z), v.var1.z),
                          collapse = "  ")))
        cat(sprintf("    posfr  %.4f\n", s.pos))
    }
    if(!is.na(l.y[["B.y"]]))
        cat(sprintf("    B.y    %.1f   (%s)\n", l.y[["B.y"]], l.y[["source"]]))
    v.sig.eig <- eigen(m.Sigma, symmetric = TRUE, only.values = TRUE)[["values"]]
    cat(sprintf("    Sigma  E[xx'] eigenvalues in [%.4g, %.4g]\n",
                min(v.sig.eig), max(v.sig.eig)))

    list(scenario = sc, family = fam, link = lnk,
         beta = v.beta, sigma2 = s.sigma2, gamma = v.gamma,
         beta.var1 = v.var1, gamma.var1 = v.var1.z, pos.frac = s.pos,
         Sigma = m.Sigma,
         B.y = l.y[["B.y"]], y.clip = l.y[["y.clip"]],
         y.source = l.y[["source"]],
         source = sprintf(paste0("population value under the published DGP ",
                                 "(00-simulate.r, Gaussian copula), N = %g, ",
                                 "public seed %d, independent of the n = 2000 ",
                                 "seed-%d analysis sample"),
                          N.POP, PUB.SEED, ANALYSIS.SEED))
}, 1:11)
l.const <- `names<-`(l.const, paste0("sc", 1:11))

saveRDS(l.const, "results/dp-design-constants.rds")
cat(sprintf("\nWritten to results/dp-design-constants.rds (%d scenarios)\n",
            length(l.const)))
cat("\nEvery value above is a functional of the DGP's published parameters.\n")
cat("None of them is a statistic of the analysis sample.  That is what makes them public.\n")
