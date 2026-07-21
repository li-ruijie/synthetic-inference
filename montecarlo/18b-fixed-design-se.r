# +------------------------------------------------------------------------+
# | 18b-fixed-design-se.r -- the projection's standard errors under the     |
# |                          theorem's fixed design versus full synthesis   |
# +------------------------------------------------------------------------+
#
# Why this script exists.
#
# The projection theorem's part (c) states that it preserves the analyst's standard
# errors exactly.  The simulation reports the projection's SE deviation as 1.0
# to 6.9%, which reads as the simulation refuting the theorem.  It does not.
# The theorem fixes the design X and synthesizes only the response y*, whereas
# gen.syn.projected draws X* from the copula too, so the analyst refits on
# (X*, y*) and the information matrix (X*'W X*)^-1 is formed on a different
# design.  The coefficients survive, since the projection targets the original
# sufficient statistic, but the standard errors move with the design.
#
# This script quantifies the two regimes side by side, single run at seed 6489,
# to match the SE column of section 1.  "full" is the released fully-synthetic
# table.  "fixed" releases the original covariates and projects only the
# response, the fixed-design public-covariate setting of the projection theorem.
#
# Output:  fixed-design-se.csv
# Usage:   Rscript 18b-fixed-design-se.r
#
# Both source-guard flags are required.  07- is not needed here, but 11-
# re-sources 00-shared.r, which clears the global environment unless
# .shared.loaded is set first (see 25-mc-allmethod.r, which sets both).
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

fit.model <- function(d, fam, lnk) {
    # zero-drop for the two-part quasi families only, and Poisson keeps its zeros
    d.fit <- if (lnk == "log" && any(d[["response"]] == 0)) d[d[["response"]] > 0, ] else d
    if (fam == "gaussian")
        lm(response ~ age + female + race + income, d.fit)
    else if (fam == "poisson")
        glm(response ~ age + female + race + income, d, family = poisson())
    else
        glm(response ~ age + female + race + income, d.fit,
            family = quasi(link = "log", variance = "mu^2"),
            start = c(1, rep(0, 5)), control = glm.control(maxit = 100L))
}

maxpct <- function(a, b) max(abs(100 * (a - b) / b), na.rm = TRUE)

cfg <- list(
    list(i = 1L, sup = c(response = "real",     age = "real", income = "real")),
    list(i = 3L, sup = c(response = "positive", age = "real", income = "positive")),
    list(i = 4L, sup = c(response = "zeroinfl", age = "real", income = "positive")),
    list(i = 6L, sup = c(response = "count",    age = "real", income = "positive")),
    list(i = 8L, sup = c(response = "zeroinfl", age = "real", income = "positive")),
    list(i = 9L, sup = c(response = "bounded",  age = "real", income = "positive")))

load("results/d.sim.scenarios.RData")
v.x <- c("age", "female", "race", "income")

l.rows <- Map(function(cc) {
    s   <- l.scenarios[[cc$i]]
    d   <- s[["data"]]
    fam <- s[["family"]]
    lnk <- s[["link"]]
    mo.o <- fit.model(d, fam, lnk)
    se.o <- summary(mo.o)[["coefficients"]][, 2]
    b.o  <- coef(mo.o)

    # Full synthesis: covariates X* drawn from the copula (the released table)
    set.seed(SEED)
    d.full <- gen.syn.projected(d, y.nm = "response", family = fam, link = lnk,
                                supports = cc$sup, seed = SEED)
    mo.f <- fit.model(d.full, fam, lnk)

    # Fixed design: original covariates, copula-drawn response projected onto
    # the original X's sufficient-statistic constraint (public-covariate release)
    set.seed(SEED)
    d.cop <- gen.syn.family(d, order = 4, supports = cc$sup, seed = SEED)
    d.mix <- d
    d.mix[["response"]] <- d.cop[["response"]]
    y.fix <- project.glm(d.mix, d, "response", v.x, fam, lnk)
    d.fix <- d
    d.fix[["response"]] <- y.fix
    mo.x <- fit.model(d.fix, fam, lnk)

    se.f <- summary(mo.f)[["coefficients"]][, 2]
    se.x <- summary(mo.x)[["coefficients"]][, 2]

    row <- data.frame(
        scenario  = cc$i, family = fam, link = lnk,
        zero_infl = grepl("zeroinfl", cc$sup[["response"]]),
        coef_dev_full = maxpct(coef(mo.f), b.o),
        se_dev_full   = maxpct(se.f, se.o),
        coef_dev_fixed = maxpct(coef(mo.x), b.o),
        se_dev_fixed   = maxpct(se.x, se.o))
    cat(sprintf("sc%d %-12s/%-8s  SE full %6.3f%%  SE fixed %8.6f%%%s\n",
                cc$i, fam, lnk, row$se_dev_full, row$se_dev_fixed,
                if (row$zero_infl) "  (zero-inflated, two-part fit)" else ""))
    row
}, cfg)

d.out <- do.call(rbind, l.rows)
write.csv(d.out, "results/fixed-design-se.csv", row.names = FALSE)
cat(sprintf("\nWritten to fixed-design-se.csv (%d scenarios)\n", nrow(d.out)))
