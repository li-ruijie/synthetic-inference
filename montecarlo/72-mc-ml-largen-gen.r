# +------------------------------------------------------------------------+
# | 72-mc-ml-largen-gen.r -- generate scenario-1 draws at n in {2000, 20000}|
# |                          over 3 reps, for the larger-n ML robustness    |
# +------------------------------------------------------------------------+
# Phase 1 of the larger-n ML check (see 72-mc-ml-largen-report.r).  Redraws
# scenario 1 at each n for reps 1-3 (seeds 6489+k), writes the orig CSVs the
# Python trainer consumes, and records the reference OLS coefficients.  Output
# goes to $LARGEN_DIR (default largen/), transient and regenerable.
# Usage: Rscript 72-mc-ml-largen-gen.r
# +------------------------------------------------------------------------+
# extract sim.scenario without 00-simulate.r's bottom-of-file Map/save (00-shared.r,
# sourced inside, clears the global env, so nothing above the eval survives).
local({
    src <- readLines("00-simulate.r")
    cut <- grep("^l.scenarios <- Map\\(sim.scenario", src)[1]
    eval(parse(text = paste(src[seq_len(cut - 1L)], collapse = "\n")), envir = globalenv())
})
stopifnot(is.function(sim.scenario))
WD <- Sys.getenv("LARGEN_DIR", unset = "results/largen")
dir.create(WD, showWarnings = FALSE, recursive = TRUE)
fit1 <- function(d) { d[["female"]] <- factor(d[["female"]]); d[["race"]] <- factor(d[["race"]])
    lm(response ~ age + female + race + income, d) }
ref <- list()
for (n in c(2000L, 20000L)) for (k in 1:3) {
    d <- sim.scenario(1, n = n, seed = 6489L + k)[["data"]]
    write.csv(d, file.path(WD, sprintf("results/largen_sc01_n%d_rep%d_orig.csv", n, k)), row.names = FALSE)
    ref[[sprintf("%d_%d", n, k)]] <- coef(fit1(d))
}
saveRDS(ref, file.path(WD, "largen_reps_refcoefs.rds"))
cat("gen done:", length(ref), "orig files in", WD, "\n")
