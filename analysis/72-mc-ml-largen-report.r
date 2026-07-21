# +------------------------------------------------------------------------+
# | 72-mc-ml-largen-report.r -- larger-n ML robustness, per-run summary     |
# +------------------------------------------------------------------------+
#
# Referees object that n = 2000 is too small for CTGAN/TVAE, so the failure
# might be a small-sample artefact.  This re-runs the four SDV generators on
# scenario 1 at n = 20000 (10x) over three reps and reports the per-run maximum
# coefficient deviation.  The generators are high variance, so the reported
# statistic is the minimum over runs (a robust lower bound), not a mean, and
# this experiment redraws the DGP each rep, so it is an independent check, not
# comparable to the fixed-data 100-replication MC of 25-mc-allmethod.r.
#
# Pipeline: 72-mc-ml-largen-gen.r -> 72-mc-ml-largen.py -> this report.
#
# Input:  $LARGEN_DIR/largen_sc01_n{2000,20000}_rep{1,2,3}_{method}.csv and
#         $LARGEN_DIR/largen_reps_refcoefs.rds (LARGEN_DIR default largen)
# Output: mc-ml-largen.csv (24 rows: n x method x rep -> maxdev)
# Usage:  Rscript 72-mc-ml-largen-report.r
# +------------------------------------------------------------------------+
WD <- Sys.getenv("LARGEN_DIR", unset = "results/largen")
ref <- readRDS(file.path(WD, "largen_reps_refcoefs.rds"))
ns <- c(2000L, 20000L); reps <- 1:3; meths <- c("ctgan", "tvae", "gauscop", "copgan")
fit1 <- function(d) {
    d[["female"]] <- factor(d[["female"]]); d[["race"]] <- factor(d[["race"]])
    lm(response ~ age + female + race + income, d)
}
rows <- list()
for (n in ns) for (m in meths) for (k in reps) {
    p <- file.path(WD, sprintf("results/largen_sc01_n%d_rep%d_%s.csv", n, k, m))
    if (!file.exists(p)) next
    mo <- tryCatch(fit1(read.csv(p)), error = function(e) NULL); if (is.null(mo)) next
    vs <- coef(mo); vo <- ref[[sprintf("%d_%d", n, k)]]; cc <- intersect(names(vo), names(vs))
    pct <- 100 * (vs[cc] - vo[cc]) / abs(vo[cc]); pct <- pct[is.finite(pct)]
    rows[[paste(n, m, k)]] <- data.frame(n = n, method = m, rep = k, maxdev = max(abs(pct)))
}
allv <- do.call(rbind, rows)
write.csv(allv, "results/mc-ml-largen.csv", row.names = FALSE)

neural <- allv[["method"]] %in% c("ctgan", "tvae", "copgan")
cat(sprintf("larger-n ML robustness: %d runs\n", nrow(allv)))
cat(sprintf("  MIN max-deviation at n = 2000  (any method): %.1f%%\n",
            min(allv[["maxdev"]][allv[["n"]] == 2000])))
cat(sprintf("  MIN max-deviation at n = 20000, GaussianCopula: %.1f%%\n",
            min(allv[["maxdev"]][allv[["n"]] == 20000 & allv[["method"]] == "gauscop"])))
cat(sprintf("  MIN max-deviation at n = 20000, neural (CTGAN/TVAE/CopulaGAN): %.1f%%\n",
            min(allv[["maxdev"]][allv[["n"]] == 20000 & neural])))
cat("Written to mc-ml-largen.csv\n")
