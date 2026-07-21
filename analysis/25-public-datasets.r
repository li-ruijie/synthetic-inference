#!/usr/bin/env Rscript
# Reproducible demonstration on public datasets.
#
# Datasets:
#   1. Adult Census Income (UCI, n=32,561) -- logistic regression
#   2. diamonds (ggplot2, n=53,940) -- gamma/log regression
#
# These allow full reproducibility without access to PHI data.
#
# This is the only script that reads results/adult.data.csv, and that file is a
# cache of the UCI download below rather than a hand-made input.  It is committed,
# so the `if(!file.exists(s.cache))` guard means a clone does not re-download and
# the run does not depend on a third-party URL whose bytes are outside our
# control.  UCI has since restructured its repository, so committing the cache is
# what pins the exact 32,561 rows this paper's real-data validation used.
#
# The dataset is redistributed here under CC by 4.0, its UCI licence.  Credit:
# Becker, B. and Kohavi, R. (1996), Adult, UCI Machine Learning Repository,
# https://doi.org/10.24432/C5XW20.  Derived from the 1994 US Census Bureau data.
#
# Usage:
#   Rscript analysis/25-public-datasets.r
# setup +------------------------------------------------------+
.sobol.main <- TRUE
source("lib/00-shared.r")
s.wd <- tryCatch(
    normalizePath(dirname(sys.frame(1)[["ofile"]]), mustWork = TRUE),
    error = function(e) getwd())
setwd(s.wd)
source("lib/gen-syn-family.r")
source("lib/11-projection-generator.r")
rm(.sobol.main)
# helpers +----------------------------------------------------+
max.abs.pct <- function(v.syn, v.orig) {
    v.pct <- abs(100 * (v.syn - v.orig) / pmax(abs(v.orig), 1e-10))
    max(v.pct, na.rm = TRUE)
}

report <- function(s.label, mo.orig, mo.cop, mo.proj) {
    v.o <- coef(mo.orig)
    v.c <- coef(mo.cop)[seq_along(v.o)]
    v.p <- coef(mo.proj)[seq_along(v.o)]
    cat(sprintf("\n%s  (p=%d)\n", s.label, length(v.o)))
    cat(sprintf("  max|%%|: copula=%.1f%%  projection=%.1f%%\n",
                max.abs.pct(v.c, v.o), max.abs.pct(v.p, v.o)))
    cat(sprintf("\n  %-25s %10s %12s %12s\n",
                "coeff", "orig", "copula", "projection"))
    cat("  ", strrep("-", 62), "\n", sep = "")
    Map(function(j) {
        s.c <- sprintf("%+.1f%%", 100 * (v.c[j] - v.o[j]) / abs(v.o[j]))
        s.p <- sprintf("%+.1f%%", 100 * (v.p[j] - v.o[j]) / abs(v.o[j]))
        cat(sprintf("  %-25s %10.4f %12s %12s\n",
                    names(v.o)[j], v.o[j], s.c, s.p))
    }, seq_along(v.o))
}
# 1. Adult Census Income +-------------------------------------+
cat("=" |> rep(64) |> paste(collapse = ""), "\n")
cat("ADULT CENSUS INCOME (UCI)\n")
cat("=" |> rep(64) |> paste(collapse = ""), "\n")

s.url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/adult/adult.data"
s.cache <- "results/adult.data.csv"
if(!file.exists(s.cache)) {
    # 60s, R's default, is far too short for a multi-MB third-party download
    op <- options(timeout = max(900, getOption("timeout")))
    on.exit(options(op), add = TRUE)
    download.file(s.url, s.cache, mode = "wb")
}
v.names <- c("age", "workclass", "fnlwgt", "education", "education_num",
             "marital", "occupation", "relationship", "race", "sex",
             "capital_gain", "capital_loss", "hours", "country", "income")
d.adult <- read.csv(s.cache, header = FALSE, strip.white = TRUE,
                    col.names = v.names, na.strings = "?")
d.adult <- d.adult[complete.cases(d.adult), ]
cat(sprintf("n=%d, %d columns\n", nrow(d.adult), ncol(d.adult)))

# Subsample for tractability (copula on 32k rows is slow)
set.seed(6489)
d.ad <- d.adult[sample(nrow(d.adult), 5000), ]

# Simplify: income ~ age + sex + race + education_num + hours
d.ad <- data.frame(
    income = factor(ifelse(d.ad[["income"]] == ">50K", 1, 0)),
    age    = d.ad[["age"]],
    sex    = factor(d.ad[["sex"]]),
    race   = factor(ifelse(d.ad[["race"]] == "White", "White",
                    ifelse(d.ad[["race"]] == "Black", "Black", "Other"))),
    edu    = d.ad[["education_num"]],
    hours  = d.ad[["hours"]])

# logistic regression
fml.ad <- income ~ age + sex + race + edu + hours
mo.ad <- glm(fml.ad, d.ad, family = binomial())
cat("Model: income ~ age + sex + race + edu + hours (logistic)\n")
cat(sprintf("n=%d, p=%d\n", nrow(d.ad), length(coef(mo.ad))))

# copula
sup.ad <- c(income = "binary", age = "real", edu = "real", hours = "real")
d.ad.cop <- gen.syn.family(d.ad, order = 4, supports = sup.ad, seed = 6489)

# projection (binomial/logit)
v.x.ad <- setdiff(names(d.ad), "income")
# convert income to numeric for projection
d.ad.num <- d.ad
d.ad.num[["income"]] <- as.numeric(as.character(d.ad[["income"]]))
d.ad.cop.num <- d.ad.cop
d.ad.cop.num[["income"]] <- as.numeric(as.character(d.ad.cop[["income"]]))
# for binary data, projection = matching score X'(y-mu) = 0
# use project.glm with binomial family
v.y.proj <- tryCatch(
    project.glm(d.ad.cop.num, d.ad.num, "income", v.x.ad,
                "binomial", "logit"),
    error = function(e) { cat("  projection error:", e[["message"]], "\n"); NULL })

if(!is.null(v.y.proj)) {
    # The release is the fractional response in [0, 1] that project.glm returns,
    # which satisfies X'y* = X'y exactly (11-projection-generator.r).  Do not
    # threshold it at 0.5: a deterministic round is not measure-preserving, it
    # breaks that constraint, and thresholding here scores the projection at
    # 7243.9% on a canonical link.  glm(binomial) warns
    # "non-integer #successes" on a fractional response, which is expected and
    # harmless at weights = 1, so suppress it, as 60-mc-binary.r and the package's
    # test-projection.R both do.
    d.ad.proj <- d.ad.cop.num
    d.ad.proj[["income"]] <- v.y.proj
    mo.ad.cop <- tryCatch(glm(fml.ad, d.ad.cop, family = binomial()),
                          error = function(e) NULL)
    mo.ad.proj <- tryCatch(
        suppressWarnings(glm(fml.ad, d.ad.proj, family = binomial())),
        error = function(e) NULL)
    if(!is.null(mo.ad.cop) && !is.null(mo.ad.proj))
        report("Adult Income (logistic)", mo.ad, mo.ad.cop, mo.ad.proj)
} else {
    cat("  Projection failed for binary outcome. Using copula only.\n")
    mo.ad.cop <- tryCatch(glm(fml.ad, d.ad.cop, family = binomial()),
                          error = function(e) NULL)
    if(!is.null(mo.ad.cop)) {
        v.o <- coef(mo.ad)
        v.c <- coef(mo.ad.cop)[seq_along(v.o)]
        cat(sprintf("  copula max|%%|: %.1f%%\n", max.abs.pct(v.c, v.o)))
    }
}

# 2. diamonds +------------------------------------------------+
cat("\n\n")
cat("=" |> rep(64) |> paste(collapse = ""), "\n")
cat("DIAMONDS (ggplot2)\n")
cat("=" |> rep(64) |> paste(collapse = ""), "\n")

if(!requireNamespace("ggplot2", quietly = TRUE)) {
    cat("ggplot2 not installed, skipping diamonds\n")
} else {
    d.dia <- ggplot2::diamonds
    # Subsample for tractability
    set.seed(6489)
    d.dia <- d.dia[sample(nrow(d.dia), 5000), ]

    # Simplify: price ~ carat + cut + color + clarity + depth
    d.dia <- data.frame(
        price   = d.dia[["price"]],
        carat   = d.dia[["carat"]],
        cut     = factor(d.dia[["cut"]]),
        color   = factor(d.dia[["color"]]),
        clarity = factor(d.dia[["clarity"]]),
        depth   = d.dia[["depth"]])

    fml.dia <- price ~ carat + cut + color + clarity + depth
    # gamma/log is natural for positive skewed price
    mo.dia <- glm(fml.dia, d.dia, family = Gamma(link = "log"))
    cat("Model: price ~ carat + cut + color + clarity + depth (gamma/log)\n")
    cat(sprintf("n=%d, p=%d\n", nrow(d.dia), length(coef(mo.dia))))

    sup.dia <- c(price = "positive", carat = "positive", depth = "real")
    d.dia.cop <- gen.syn.family(d.dia, order = 4, supports = sup.dia,
                                seed = 6489)

    v.x.dia <- setdiff(names(d.dia), "price")
    v.y.proj.dia <- project.glm(d.dia.cop, d.dia, "price", v.x.dia,
                                "gamma", "log")
    d.dia.proj <- d.dia.cop
    d.dia.proj[["price"]] <- v.y.proj.dia

    mo.dia.cop <- tryCatch(
        glm(fml.dia, d.dia.cop, family = Gamma(link = "log")),
        error = function(e) NULL)
    mo.dia.proj <- tryCatch(
        glm(fml.dia, d.dia.proj, family = Gamma(link = "log")),
        error = function(e) NULL)

    if(!is.null(mo.dia.cop) && !is.null(mo.dia.proj))
        report("Diamonds (gamma/log)", mo.dia, mo.dia.cop, mo.dia.proj)
}
