#!/usr/bin/env r
# +------------------------------------------------------------------------+
# | 16-report-privacy.r -- Privacy assessment report                       |
# +------------------------------------------------------------------------+
#
# Assesses the privacy properties of the projection generator against
# the copula baseline. This addresses Remark (projection privacy) in
# the manuscript, which decomposes the projected data into:
#   1. Fitted component X*beta_hat: embeds the sufficient statistic
#      X'y (p numbers), with the same disclosure risk as publishing
#      the regression output itself.
#   2. Residual component c*(I-H)*y_cop: carries n-p-1 degrees of
#      freedom of genuine novelty from the copula. Not copies of
#      original values.
#
# Metrics reported:
#   DCR (Distance to Closest Record): Gower distance from each synthetic
#     row to its nearest original row. The paper reports no DCR figure
#     and does not name the metric, so "no exact copies" and the
#     copula-projection comparison below are this script's own results.
#   Hitting rate: fraction of synthetic rows that are exact copies
#     (should be 0% for both generators).
#   k-map risk: quasi-identifier overlap -- measures how many original
#     records share the same categorical pattern as each synthetic record.
#   Membership inference: 5-fold CV logistic classifier AUC for
#     distinguishing original from synthetic. AUC ~ 0.5 means
#     indistinguishable. The paper puts the membership-inference area
#     under the receiver operating characteristic curve at "approximately
#     0.49, indistinguishable from chance" [Shokri et al. 2017,
#     shokri-membersh-2017].
#
# The projection is not differentially private (Remark: projection
# privacy). Adding one observation changes X'y and the entire fitted
# component. The natural DP route is Gaussian mechanism output
# perturbation [dwork-algorith-2014], but the analytical pre-check
# (20-dp-precheck.r) shows this is not viable for regression inference.
#
# +------------------------------------------------------------------------+
# | Setup                                                                  |
# +------------------------------------------------------------------------+
.sobol.main <- TRUE
source("lib/00-shared.r")
s.wd <- tryCatch(
    normalizePath(dirname(sys.frame(1)[["ofile"]]), mustWork = TRUE),
    error = function(e) getwd())
setwd(s.wd)
source("lib/gen-syn-family.r")
source("lib/07-calibrated-generator.r")
source("lib/09-symmetric-generator.r")
source("lib/11-projection-generator.r")
rm(.sobol.main)
# +------------------------------------------------------------------------+
# | DCR (Gower distance)                                                   |
# +------------------------------------------------------------------------+
# Gower distance normalises each variable to [0,1] and averages.
# For the projection, the DCR is bounded below by the copula's
# per-stratum conditional variance (Remark: projection privacy),
# which is positive whenever the stratum has >= 2 observations.
dcr <- function(d.syn, d.orig, n.sample = 500) {
    # sample for speed (O(n^2) otherwise)
    s.ns <- min(n.sample, nrow(d.syn))
    v.idx <- sample(nrow(d.syn), s.ns)

    v.types <- mapply(id.cat, d.orig)
    v.cnt <- names(v.types[v.types == "cnt"])
    v.cat <- names(v.types[v.types != "cnt"])

    # range for normalisation
    l.range <- Map(function(nm) {
        s.r <- diff(range(d.orig[[nm]], na.rm = TRUE))
        `if`(s.r > 0, s.r, 1)
    }, v.cnt)
    s.p <- length(v.cnt) + length(v.cat)

    mapply(function(i) {
        s.min <- Inf
        Map(function(j) {
            s.d <- 0
            Map(function(nm) {
                s.d <<- s.d + abs(d.syn[[nm]][i] - d.orig[[nm]][j]) /
                        l.range[[nm]]
            }, v.cnt)
            Map(function(nm) {
                s.d <<- s.d + as.integer(
                    as.character(d.syn[[nm]][i]) !=
                    as.character(d.orig[[nm]][j]))
            }, v.cat)
            s.min <<- min(s.min, s.d / s.p)
        }, seq_len(nrow(d.orig)))
        s.min
    }, v.idx)
}
# +------------------------------------------------------------------------+
# | k-map risk                                                             |
# +------------------------------------------------------------------------+
# Quasi-identifier overlap: for each synthetic record, the number of original
# records sharing the same categorical pattern. Low k = high re-identification
# risk. This is the stratum-level privacy metric. Sparse strata (n_s=1)
# from Theorem 5 are the most vulnerable -- the pre-masking strategy
# (22-premasking-demo.r) addresses this by collapsing rare levels.
kmap.risk <- function(d.syn, d.orig) {
    v.types <- mapply(id.cat, d.orig)
    v.cat <- names(v.types[v.types != "cnt"])
    if(length(v.cat) == 0) return(rep(nrow(d.orig), nrow(d.syn)))

    v.key.orig <- interaction(d.orig[v.cat], drop = TRUE)
    v.key.syn  <- interaction(d.syn[v.cat], drop = TRUE)
    v.tab <- table(v.key.orig)

    mapply(function(k) {
        s.k <- as.character(k)
        `if`(s.k %in% names(v.tab), v.tab[[s.k]], 0)
    }, v.key.syn)
}
# +------------------------------------------------------------------------+
# | Membership inference                                                   |
# +------------------------------------------------------------------------+
# Train a logistic classifier to distinguish original (1) from synthetic
# (0) records. AUC ~ 0.5 means the classifier cannot distinguish them,
# implying the synthetic data is statistically indistinguishable.
# This is the membership inference attack framework of Shokri et al.
# (2017) [shokri-membersh-2017]. The paper reports AUC ~0.49.
membership.auc <- function(d.syn, d.orig) {
    # train logistic classifier on the label: original (1) or synthetic (0).
    # use all variables. AUC ~ 0.5 = indistinguishable.
    s.n <- min(nrow(d.syn), nrow(d.orig))
    d.both <- rbind(
        cbind(d.orig[sample(nrow(d.orig), s.n), ], .real = 1),
        cbind(d.syn[sample(nrow(d.syn), s.n), ],   .real = 0))
    # simple 5-fold cross-validated AUC
    v.pred <- numeric(nrow(d.both))
    v.fold <- sample(rep(1:5, length.out = nrow(d.both)))
    Map(function(k) {
        b.test  <- v.fold == k
        b.train <- !b.test
        mo <- tryCatch(
            glm(.real ~ ., data = d.both[b.train, ],
                family = binomial(), control = glm.control(maxit = 50)),
            error = function(e) NULL)
        if(!is.null(mo)) {
            v.pred[b.test] <<- tryCatch(
                predict(mo, d.both[b.test, ], type = "response"),
                error = function(e) rep(0.5, sum(b.test)))
        } else {
            v.pred[b.test] <<- 0.5
        }
    }, 1:5)
    # AUC by concordance
    v.real <- d.both[[".real"]]
    v.pos <- v.pred[v.real == 1]
    v.neg <- v.pred[v.real == 0]
    s.n.pairs <- min(2000, length(v.pos)) * min(2000, length(v.neg))
    v.p.samp <- sample(v.pos, min(2000, length(v.pos)))
    v.n.samp <- sample(v.neg, min(2000, length(v.neg)))
    mean(outer(v.p.samp, v.n.samp, ">") +
         0.5 * outer(v.p.samp, v.n.samp, "=="))
}
# +------------------------------------------------------------------------+
# | Run assessments                                                        |
# +------------------------------------------------------------------------+
run.privacy <- function(d.orig, d.syn, label) {
    cat(sprintf("--- %s ---\n", label))

    # 1. DCR
    v.dcr <- dcr(d.syn, d.orig, n.sample = 300)
    cat(sprintf("  DCR:  median=%.4f  5th pct=%.4f  min=%.4f\n",
                median(v.dcr), quantile(v.dcr, 0.05), min(v.dcr)))

    # 2. Hitting rate (exact copies)
    s.hit <- mean(v.dcr < 1e-10)
    cat(sprintf("  Hitting rate: %.1f%% (exact copies)\n", 100 * s.hit))

    # 3. k-map risk
    v.kmap <- kmap.risk(d.syn, d.orig)
    cat(sprintf("  k-map: median k=%d  min k=%d  1-anon=%d%%\n",
                median(v.kmap), min(v.kmap),
                round(100 * mean(v.kmap == 1))))

    # 4. Membership inference
    s.auc <- membership.auc(d.syn, d.orig)
    cat(sprintf("  Membership AUC: %.3f (0.5=indistinguishable)\n", s.auc))

    cat("\n")
}

cat("PRIVACY ASSESSMENT REPORT\n")
cat(strrep("=", 80), "\n\n")

## simulation Scenario 1 +--------------------------------------+
load("results/d.sim.scenarios.RData")
d.sc1 <- l.scenarios[[1]][["data"]]
sup1  <- c(response = "real", age = "real", income = "real")

cat("=== Scenario 1: Gaussian (n=2000) ===\n\n")
d.c1 <- gen.syn.family(d.sc1, order = 4, supports = sup1)
d.p1 <- gen.syn.projected(d.sc1, y.nm = "response",
                          family = "gaussian", link = "identity",
                          supports = sup1)
run.privacy(d.sc1, d.c1, "copula")
run.privacy(d.sc1, d.p1, "projection")

## simulation Scenario 4 (ZI) +----------------------------------+
d.sc4 <- l.scenarios[[4]][["data"]]
sup4  <- c(response = "zeroinfl", age = "real", income = "positive")

cat("=== Scenario 4: ZI quasi-gamma (n=2000) ===\n\n")
d.c4 <- gen.syn.family(d.sc4, order = 4, supports = sup4)
d.p4 <- gen.syn.projected(d.sc4, y.nm = "response",
                          family = "quasi-gamma", link = "log",
                          supports = sup4)
run.privacy(d.sc4, d.c4, "copula")
run.privacy(d.sc4, d.p4, "projection")

# +------------------------------------------------------------------------+
# | Real data: PHI waist circumference                                     |
# +------------------------------------------------------------------------+
# The PHI dataset [yap-developm-2020] with its 26 singleton strata
# is the hardest privacy test. Remark (projection privacy) notes:
# "sparse housing categories produce singleton cross-strata with
# leverage approaching 1, yielding single-person influence exceeding
# 1000% on housing coefficients." The DCR and membership metrics
# here quantify the residual privacy despite this vulnerability.
# The PHI dataset is restricted.  Everything above runs without it.  The
# application section below needs it, so stop cleanly when it is absent.
if (!file.exists("phi_phase1_baseline.dta")) {
    cat("\nPHI data not available (restricted). Skipping the application section.\n")
    quit(save = "no", status = 0)
}

library(haven)
d.raw <- read_dta("phi_phase1_baseline.dta")
d.app <- data.frame(
    waist   = as.numeric(d.raw[["waist_circum"]]),
    age     = 2016 - as.numeric(d.raw[["birthyear"]]),
    female  = factor(as.numeric(d.raw[["female"]])),
    race    = factor(as.numeric(d.raw[["race"]])),
    edu     = factor(findInterval(as.numeric(d.raw[["education"]]),
                                  c(5, 8)) + 1),
    housing = factor(as.numeric(d.raw[["housingtype"]])))
d.app <- d.app[complete.cases(d.app), ]

sup.app <- c(waist = "real", age = "real")

cat("=== PHI Baseline: waist circumference (n=1927) ===\n\n")
d.c.app <- gen.syn.family(d.app, order = 4, supports = sup.app)
d.p.app <- gen.syn.projected(d.app, y.nm = "waist",
                             family = "gaussian", link = "identity",
                             supports = sup.app)
run.privacy(d.app, d.c.app, "copula")
run.privacy(d.app, d.p.app, "projection")
# +------------------------------------------------------------+
