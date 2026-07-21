# +------------------------------------------------------------------------+
# | 18c-twopart-se.r -- the two-part projection: exact standard errors on  |
# |                     zero-inflated data when the zero pattern is fixed   |
# +------------------------------------------------------------------------+
#
# Why this script exists.
#
# 18b-fixed-design-se.r showed the projection's standard errors go exact under
# the theorem's fixed design for every scenario except the two zero-inflated
# ones (Sc.4, Sc.8), which sit at 7 to 10%.  The reason is that the response is
# fitted as a two-part (hurdle) model on the positive subset, and the copula
# resamples which rows are zero, so the positive-subset design shifts even with
# the full design fixed.
#
# The two-part projection closes this.  Preserve the zero pattern (release which
# rows are zero unchanged) and project the positive component onto the original
# positive design.  Then the zero block is exact by identity (same indicator,
# same X), and the positive block is an ordinary sufficient-statistics projection on
# the shared positive design, exact in coefficients and standard errors.  This is the
# "Proposition (two-part projection)" of S13, proved there and implemented by
# gen.syn.projected.twopart (11-projection-generator.r).  This
# script drives that generator and checks the positive-block SE goes exact.
#
# Output:  twopart-se.csv
# Usage:   Rscript 18c-twopart-se.r
#
# Both source-guard flags are required (see 18b-, 25-mc-allmethod.r).
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

fit.pos <- function(d) {
    dp <- d[d[["response"]] > 0, ]
    glm(response ~ age + female + race + income, dp,
        family = quasi(link = "log", variance = "mu^2"),
        start = c(1, rep(0, 5)), control = glm.control(maxit = 100L))
}
fit.zero <- function(d) {
    d[["z"]] <- as.integer(d[["response"]] > 0)
    glm(z ~ age + female + race + income, d, family = binomial())
}
maxpct <- function(a, b) max(abs(100 * (a - b) / b), na.rm = TRUE)

sup <- c(response = "zeroinfl", age = "real", income = "positive")

load("results/d.sim.scenarios.RData")

l.rows <- Map(function(sc) {
    d <- l.scenarios[[sc]][["data"]]
    P <- d[["response"]] > 0                         # ORIGINAL zero/positive partition

    mo.pos.o  <- fit.pos(d);  se.pos.o  <- summary(mo.pos.o)[["coefficients"]][, 2]
    mo.zero.o <- fit.zero(d); se.zero.o <- summary(mo.zero.o)[["coefficients"]][, 2]

    # The two-part release now lives in gen.syn.projected.twopart, which keeps
    # the original zero pattern, injects copula novelty on the positive rows,
    # and projects the positive component onto the original X_+.
    d.fix <- gen.syn.projected.twopart(d, y.nm = "response", family = "quasi-gamma",
                                       link = "log", supports = sup, seed = SEED)
    stopifnot(all((d.fix[["response"]] > 0) == P))   # zero pattern preserved exactly

    mo.pos.x  <- fit.pos(d.fix)
    mo.zero.x <- fit.zero(d.fix)

    row <- data.frame(
        scenario = sc,
        pos_coef_dev = maxpct(coef(mo.pos.x), coef(mo.pos.o)),
        pos_se_dev   = maxpct(summary(mo.pos.x)[["coefficients"]][, 2], se.pos.o),
        zero_coef_dev = maxpct(coef(mo.zero.x), coef(mo.zero.o)),
        zero_se_dev   = maxpct(summary(mo.zero.x)[["coefficients"]][, 2], se.zero.o))
    cat(sprintf("sc%d  positive block: coef %.6f%% SE %.6f%%   zero block: coef %.6f%% SE %.6f%%\n",
                sc, row$pos_coef_dev, row$pos_se_dev, row$zero_coef_dev, row$zero_se_dev))
    row
}, c(4L, 8L))

d.out <- do.call(rbind, l.rows)
write.csv(d.out, "results/twopart-se.csv", row.names = FALSE)
cat(sprintf("\nWritten to twopart-se.csv (%d scenarios)\n", nrow(d.out)))
