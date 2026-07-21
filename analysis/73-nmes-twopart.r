# +------------------------------------------------------------------------+
# | 73-nmes-twopart.r -- public two-part application on AER::NMES1988       |
# +------------------------------------------------------------------------+
#
# Why this script exists.
#
# The paper's real-data two-part exhibit (supplement S4) runs on NMES1988, Deb
# and Trivedi's demand-for-medical-care data from the AER package: outpatient
# visits carry a 77% zero fraction and the demographic cross-classification
# carries singleton strata, the sparse-strata regime of Section 5.
#
# Two-part model of outpatient visits (ovisits) on health status, chronic
# conditions, insurance, and demographics.  Three arms: original fit, copula
# release (gen.syn.family), and the two-part projection release
# (gen.syn.projected.twopart), which keeps the zero pattern and projects the
# positive block onto the original positive design.  The grid records, per
# coefficient and block, the sign and the 5% significance decision under each
# arm; the substantive parameters read in the paper are the insurance and
# health-status coefficients.
#
# Output:  results/nmes-twopart.csv         (per-term grid, all three arms)
#          results/nmes-twopart-strata.csv  (cross-strata counts)
# Usage:   Rscript analysis/73-nmes-twopart.r
#
# Both source-guard flags are load-bearing (see montecarlo/25-mc-allmethod.r).
.shared.loaded <- TRUE
.sobol.main <- TRUE
source("lib/00-shared.r")

s.wd <- tryCatch(
    normalizePath(dirname(sys.frame(1)[["ofile"]]), mustWork = TRUE),
    error = function(e) getwd())
setwd(s.wd)

source("lib/gen-syn-family.r")
source("lib/11-projection-generator.r")
rm(.sobol.main)

SEED <- 6489L

# data +-------------------------------------------------------+
data("NMES1988", package = "AER")
# The AER factors carry explicit contrasts attributes (reference "average" for
# health, "other" for region), which gen.syn.family drops with a warning, and a
# dropped attribute moves the reference level and renames coefficients.  Every
# factor is normalised here to plain default coding, keeping the dataset's own
# reference choices by level order instead.
d <- data.frame(
    ovisits   = NMES1988[["ovisits"]],
    health    = factor(as.character(NMES1988[["health"]]),
                       levels = c("average", "poor", "excellent")),
    chronic   = NMES1988[["chronic"]],
    insurance = factor(as.character(NMES1988[["insurance"]])),
    age       = NMES1988[["age"]],
    gender    = factor(as.character(NMES1988[["gender"]])),
    married   = factor(as.character(NMES1988[["married"]])),
    employed  = factor(as.character(NMES1988[["employed"]])),
    medicaid  = factor(as.character(NMES1988[["medicaid"]])),
    region    = factor(as.character(NMES1988[["region"]]),
                       levels = c("other", "northeast", "midwest", "west")),
    school    = NMES1988[["school"]],
    income    = NMES1988[["income"]])
d <- d[complete.cases(d), ]
cat(sprintf("n = %d, zero fraction of ovisits = %.3f\n",
            nrow(d), mean(d[["ovisits"]] == 0)))

# A release may return a factor re-levelled alphabetically, which moves the
# reference level and renames coefficients.  Re-impose the original level sets.
relevel.like <- function(d.rel, d.ref) {
    v.fac <- Filter(function(nm) is.factor(d.ref[[nm]]), names(d.ref))
    Reduce(function(dd, nm)
        `[[<-`(dd, nm, value = factor(as.character(dd[[nm]]),
                                      levels = levels(d.ref[[nm]]))),
        v.fac, d.rel)
}

# cross-strata of the categorical predictors, for the sparse-strata sentence
v.strata.nm <- c("health", "insurance", "gender", "married", "employed",
                 "medicaid", "region")
t.cross <- table(do.call(interaction, c(d[v.strata.nm], drop = TRUE)))
d.strata <- data.frame(
    vars       = paste(v.strata.nm, collapse = ":"),
    cells      = length(t.cross),
    singletons = sum(t.cross == 1),
    doubletons = sum(t.cross == 2),
    min.cell   = min(t.cross))
write.csv(d.strata, "results/nmes-twopart-strata.csv", row.names = FALSE)
cat(sprintf("cross-strata %s: %d cells, %d singletons, %d doubletons, min %d\n",
            d.strata[["vars"]], d.strata[["cells"]], d.strata[["singletons"]],
            d.strata[["doubletons"]], d.strata[["min.cell"]]))

# fits +-------------------------------------------------------+
fml.pos  <- ovisits ~ health + chronic + insurance + age + gender +
    married + employed + medicaid + region + school + income
fml.zero <- z ~ health + chronic + insurance + age + gender +
    married + employed + medicaid + region + school + income

fit.pos <- function(d) {
    dp <- d[d[["ovisits"]] > 0, ]
    glm(fml.pos, dp, family = quasi(link = "log", variance = "mu^2"),
        start = c(1, rep(0, 14)), control = glm.control(maxit = 100L))
}
fit.zero <- function(d) {
    d[["z"]] <- as.integer(d[["ovisits"]] > 0)
    glm(fml.zero, d, family = binomial())
}

# releases +---------------------------------------------------+
d.cop <- gen.syn.family(d, order = 4, seed = SEED) |> relevel.like(d)
d.proj <- gen.syn.projected.twopart(d, y.nm = "ovisits",
                                    family = "quasi-gamma", link = "log",
                                    seed = SEED) |> relevel.like(d)
stopifnot(all((d.proj[["ovisits"]] > 0) == (d[["ovisits"]] > 0)))

# grid +-------------------------------------------------------+
grid.block <- function(s.block, fit) {
    mo.o <- fit(d); mo.c <- fit(d.cop); mo.p <- fit(d.proj)
    tab  <- function(mo) summary(mo)[["coefficients"]]
    t.o <- tab(mo.o); t.c <- tab(mo.c); t.p <- tab(mo.p)
    stopifnot(identical(rownames(t.o), rownames(t.c)),
              identical(rownames(t.o), rownames(t.p)))
    data.frame(
        block     = s.block,
        term      = rownames(t.o),
        est.orig  = t.o[, 1], se.orig = t.o[, 2], p.orig = t.o[, 4],
        est.cop   = t.c[, 1], se.cop  = t.c[, 2], p.cop  = t.c[, 4],
        est.proj  = t.p[, 1], se.proj = t.p[, 2], p.proj = t.p[, 4],
        dev.cop.pct  = 100 * abs(t.c[, 1] - t.o[, 1]) / pmax(abs(t.o[, 1]), 1e-10),
        dev.proj.pct = 100 * abs(t.p[, 1] - t.o[, 1]) / pmax(abs(t.o[, 1]), 1e-10),
        sign.agree.cop  = sign(t.c[, 1]) == sign(t.o[, 1]),
        sig.agree.cop   = (t.c[, 4] < 0.05) == (t.o[, 4] < 0.05),
        sign.agree.proj = sign(t.p[, 1]) == sign(t.o[, 1]),
        sig.agree.proj  = (t.p[, 4] < 0.05) == (t.o[, 4] < 0.05),
        row.names = NULL)
}

d.out <- rbind(grid.block("zero", fit.zero),
               suppressWarnings(grid.block("positive", fit.pos)))
write.csv(d.out, "results/nmes-twopart.csv", row.names = FALSE)

# summary +----------------------------------------------------+
smy <- function(s.block) {
    g <- d.out[d.out[["block"]] == s.block, ]
    cat(sprintf(
        paste0("%s block: copula max|%%| = %.1f, projection max|%%| = %.6f, ",
               "copula flips: sign %d, sig %d of %d; projection flips: ",
               "sign %d, sig %d\n"),
        s.block, max(g[["dev.cop.pct"]]), max(g[["dev.proj.pct"]]),
        sum(!g[["sign.agree.cop"]]), sum(!g[["sig.agree.cop"]]), nrow(g),
        sum(!g[["sign.agree.proj"]]), sum(!g[["sig.agree.proj"]])))
    g.sub <- g[grepl("health|insurance", g[["term"]]), ]
    Map(function(i) {
        r <- g.sub[i, ]
        cat(sprintf(
            "  %-22s orig %+.4f (p=%.4f)  cop %+.4f (p=%.4f)  proj %+.4f (p=%.4f)\n",
            r[["term"]], r[["est.orig"]], r[["p.orig"]],
            r[["est.cop"]], r[["p.cop"]], r[["est.proj"]], r[["p.proj"]]))
    }, seq_len(nrow(g.sub)))
    invisible(NULL)
}
smy("zero")
smy("positive")
cat("\nWritten to results/nmes-twopart.csv and results/nmes-twopart-strata.csv\n")
