# +------------------------------------------------------------------------+
# | 75-selection-commutation.r -- all-subset model selection on the         |
# |                               original data against the release         |
# +------------------------------------------------------------------------+
#
# Why this script exists.
#
# The paper states that the release commutes with model selection, resting
# on the sufficient statistics projection theorem (all submodels exact) and
# the prediction remark (information criteria match exactly for the Gaussian
# family and up to a data constant otherwise, so the ranking is preserved).
# This script backs the sentence empirically: all-subset selection over the
# 2^4 predictor subsets by AIC and BIC, original data against the release,
# on the canonical-link scenarios (1, 6, 9, 11) and the non-canonical gamma
# scenario (3), under both regimes.  "fixed" releases the original
# covariates with the projected response (the public-covariate regime of the
# projection theorem), "full" is the released fully synthetic table.  The
# strong check is ordering identity over all 16 subsets, which covers every
# selection rule based on the criterion, beyond the argmin alone.  The
# centred deviation removes the data constant by differencing each subset's
# criterion against the full model's within arm.
#
# The zero-inflated scenarios (4, 5, 8) are fitted as quasi models, which
# carry no likelihood and hence no AIC, so they are outside the sentence and
# this run.
#
# Output:  results/selection-commutation.csv
#          results/tab-supp-selection-commutation.tex
# Usage:   Rscript analysis/75-selection-commutation.r
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
V.X  <- c("age", "female", "race", "income")

# all 16 predictor subsets, intercept-only included
l.subsets <- Reduce(function(acc, nm) c(acc, Map(function(v) c(v, nm), acc)),
                    V.X, list(character(0)))
lab.sub <- function(v) if (length(v) == 0) "1" else paste(v, collapse = "+")

fam.lik <- function(fam, lnk) {
    switch(fam,
           gaussian = gaussian(),
           gamma    = Gamma(link = lnk),
           poisson  = poisson(),
           binomial = binomial())
}

fit.sub <- function(d, v.sub, fam.obj) {
    fml <- reformulate(if (length(v.sub) == 0) "1" else v.sub,
                       response = "response")
    suppressWarnings(glm(fml, d, family = fam.obj,
                         control = glm.control(maxit = 100L)))
}

crit.grid <- function(d, fam.obj) {
    l.mo <- Map(function(v) fit.sub(d, v, fam.obj), l.subsets)
    data.frame(subset = vapply(l.subsets, lab.sub, character(1)),
               aic = vapply(l.mo, AIC, numeric(1)),
               bic = vapply(l.mo, BIC, numeric(1)),
               dev = vapply(l.mo, deviance, numeric(1)))
}

cfg <- list(
    list(i = 1L,  sup = c(response = "real",     age = "real", income = "real")),
    list(i = 3L,  sup = c(response = "positive", age = "real", income = "positive")),
    list(i = 6L,  sup = c(response = "count",    age = "real", income = "positive")),
    list(i = 9L,  sup = c(response = "bounded",  age = "real", income = "positive")),
    list(i = 11L, sup = c(age = "real", income = "positive")))

load("results/d.sim.scenarios.RData")

one.arm <- function(s.arm, d.rel, g.o, fam.obj, s.i, fam, lnk) {
    g.r <- crit.grid(d.rel, fam.obj)
    stopifnot(identical(g.o[["subset"]], g.r[["subset"]]))
    s.full <- lab.sub(V.X)
    ctr <- function(g, s.col) g[[s.col]] - g[[s.col]][g[["subset"]] == s.full]
    data.frame(
        scenario  = s.i, family = fam, link = lnk, arm = s.arm,
        sel.aic.orig = g.o[["subset"]][which.min(g.o[["aic"]])],
        sel.aic.rel  = g.r[["subset"]][which.min(g.r[["aic"]])],
        sel.bic.orig = g.o[["subset"]][which.min(g.o[["bic"]])],
        sel.bic.rel  = g.r[["subset"]][which.min(g.r[["bic"]])],
        aic.same   = identical(which.min(g.o[["aic"]]), which.min(g.r[["aic"]])),
        bic.same   = identical(which.min(g.o[["bic"]]), which.min(g.r[["bic"]])),
        order.same = identical(order(g.o[["aic"]]), order(g.r[["aic"]])),
        max.centred.aic.dev = max(abs(ctr(g.r, "aic") - ctr(g.o, "aic"))),
        # the centred deviance is implementation-independent: for a fixed
        # dispersion the criterion ordering is the ordering of deviance plus
        # penalty, and centring removes the saturated terms, which are a
        # subset-independent function of the response alone
        max.centred.dev.dev = max(abs(ctr(g.r, "dev") - ctr(g.o, "dev"))))
}

l.rows <- Map(function(cc) {
    s   <- l.scenarios[[cc[["i"]]]]
    d   <- s[["data"]]
    fam <- s[["family"]]
    lnk <- s[["link"]]
    fam.obj <- fam.lik(fam, lnk)
    g.o <- crit.grid(d, fam.obj)

    # FIXED design: original covariates, copula response projected onto the
    # original X's constraint set (the public-covariate release)
    set.seed(SEED)
    d.cop <- gen.syn.family(d, order = 4, supports = cc[["sup"]], seed = SEED)
    d.mix <- `[[<-`(d, "response", value = d.cop[["response"]])
    d.fix <- `[[<-`(d, "response",
                    value = project.glm(d.mix, d, "response", V.X, fam, lnk))

    # FULL synthesis: the released fully synthetic table
    set.seed(SEED)
    d.full <- gen.syn.projected(d, y.nm = "response", family = fam, link = lnk,
                                supports = cc[["sup"]], seed = SEED)

    r.fix  <- one.arm("fixed", d.fix,  g.o, fam.obj, cc[["i"]], fam, lnk)
    r.full <- one.arm("full",  d.full, g.o, fam.obj, cc[["i"]], fam, lnk)
    cat(sprintf(paste0("sc%-2d %-8s/%-8s fixed: AIC %s BIC %s order %s ",
                       "aicdev %.3g devdev %.3g | full: AIC %s BIC %s ",
                       "order %s aicdev %.3g devdev %.3g\n"),
                cc[["i"]], fam, lnk,
                r.fix[["aic.same"]],  r.fix[["bic.same"]],
                r.fix[["order.same"]],  r.fix[["max.centred.aic.dev"]],
                r.fix[["max.centred.dev.dev"]],
                r.full[["aic.same"]], r.full[["bic.same"]],
                r.full[["order.same"]], r.full[["max.centred.aic.dev"]],
                r.full[["max.centred.dev.dev"]]))
    rbind(r.fix, r.full)
}, cfg)

d.out <- do.call(rbind, l.rows)
write.csv(d.out, "results/selection-commutation.csv", row.names = FALSE)

# emit the supplement fragment from the CSV just written
d.emit <- read.csv("results/selection-commutation.csv")
yn <- function(b) ifelse(b, "same", "differs")
dev.fmt <- function(x) ifelse(x < 1e-9, "$<10^{-9}$",
                              sprintf("$%.2g$", x)) |>
    gsub("e-0?(\\d+)", " \\\\times 10^{-\\1}", x = _) |>
    gsub("e\\+0?(\\d+)", " \\\\times 10^{\\1}", x = _)
row.tex <- function(k) {
    r <- d.emit[k, ]
    sprintf("        %d & %s (%s) & %s & %s & %s & %s & %s & %s \\\\",
            r[["scenario"]],
            switch(r[["family"]], gaussian = "Gaussian", gamma = "Gamma",
                   poisson = "Poisson", binomial = "Binomial"),
            r[["link"]], r[["arm"]], yn(r[["aic.same"]]), yn(r[["bic.same"]]),
            yn(r[["order.same"]]), dev.fmt(r[["max.centred.aic.dev"]]),
            dev.fmt(r[["max.centred.dev.dev"]]))
}
v.lines <- c(
    "% Emitted by analysis/75-selection-commutation.r from",
    "% results/selection-commutation.csv.  Do not edit by hand; rerun the script.",
    "\\begin{table}[htbp]",
    "    \\centering",
    "    \\caption{All-subset model selection on the original data against the",
    "        release: agreement of the AIC- and BIC-selected subsets, identity",
    "        of the full AIC ordering over the 16 predictor subsets, and the",
    "        maximum centred deviations of the AIC and of the deviance (each",
    "        subset's value differenced against the full model's within arm).",
    "        Fixed, the original covariates with the projected response; full,",
    "        the fully synthetic table.}",
    "    \\label{tab:supp-selection-commutation}",
    "    \\small",
    "    \\begin{tabular}{@{}clllllrr@{}}",
    "        \\toprule",
    "        Sc. & Family (link) & Release & AIC & BIC & Ordering & AIC dev & Dev.\\ dev \\\\",
    "        \\midrule",
    unlist(Map(row.tex, seq_len(nrow(d.emit)))),
    "        \\bottomrule",
    "    \\end{tabular}",
    "\\end{table}")
writeLines(v.lines, "results/tab-supp-selection-commutation.tex")
cat("Written results/selection-commutation.csv and tab-supp-selection-commutation.tex\n")
