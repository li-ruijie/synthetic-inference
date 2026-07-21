#!/usr/bin/env r
# +------------------------------------------------------------------------+
# | 60b-binary-coeftable.r -- Scenario 11 per-coefficient supplement table  |
# +------------------------------------------------------------------------+
#
# The single representative run (seed 6489, n = 2000) behind the Monte Carlo
# summary in 60-mc-binary.r, in the format of the S1 supplement tables
# (the ten per-scenario tables): the original coefficient, a copula baseline,
# and the projection, each with percentage deviation.
#
# The empirical Bayes calibrated column shown in tables S1-S10 is not defined
# for the binomial family: gen.syn.calibrated calibrates conditional moments on
# the link scale of a continuous or count response, whereas the binary response
# is categorical and bootstrapped by the copula. The baseline shown here is
# therefore the copula (gen.syn.family with shrink = TRUE), the same baseline
# behind the 20.9% main-text figure and 60-mc-binary.r.
#
# Output: mc-binary-coeftable.csv and a console table.
#
# +------------------------------------------------------------------------+
# | Setup (mirror 60-mc-binary.r: dotted .sobol.main survives the           |
# | globalenv clear in 00-shared.r, suppressing the generator harness)      |
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
# | Scenario 11 (binomial/logit): original fit, copula, projection          |
# +------------------------------------------------------------------------+
load("results/d.sim.scenarios.RData")
s   <- l.scenarios[[11]]
d   <- s[["data"]]
stopifnot(s[["family"]] == "binomial", s[["link"]] == "logit")

sup    <- c(age = "real", income = "positive")
fml    <- response ~ age + female + race + income
s.seed <- 6489

mo.o <- glm(fml, d, family = binomial())

# copula baseline (identical call to 60-mc-binary.r's copula arm)
d.cop  <- gen.syn.family(d, order = 4, supports = sup, shrink = TRUE, seed = s.seed)
mo.cop <- suppressWarnings(glm(fml, d.cop, family = binomial()))

# projection (path B: exact logistic refit expected)
d.proj  <- gen.syn.projected(d, y.nm = "response", family = "binomial",
                             link = "logit", supports = sup, seed = s.seed)
mo.proj <- suppressWarnings(glm(fml, d.proj, family = binomial()))

# +------------------------------------------------------------------------+
# | Per-coefficient deviations (aligned by coefficient name)                |
# +------------------------------------------------------------------------+
v.o    <- coef(mo.o)
v.cop  <- coef(mo.cop)[names(v.o)]
v.proj <- coef(mo.proj)[names(v.o)]
pct    <- function(v) 100 * (v - v.o) / abs(v.o)

d.out <- data.frame(
    coeff      = names(v.o),
    orig       = as.numeric(v.o),
    copula     = as.numeric(v.cop),
    copula_pct = as.numeric(pct(v.cop)),
    proj       = as.numeric(v.proj),
    proj_pct   = as.numeric(pct(v.proj)),
    row.names  = NULL, stringsAsFactors = FALSE)
write.csv(d.out, "results/mc-binary-coeftable.csv", row.names = FALSE)

# +------------------------------------------------------------------------+
# | Console table                                                           |
# +------------------------------------------------------------------------+
cat("SCENARIO 11 binomial(logit): single run, seed 6489, n = 2000\n")
cat(strrep("=", 78), "\n")
Map(function(k) {
    cat(sprintf("%-12s %9.4f  %9.4f (%+6.1f%%)  %9.4f (%+10.6f%%)\n",
                d.out[["coeff"]][k], d.out[["orig"]][k],
                d.out[["copula"]][k], d.out[["copula_pct"]][k],
                d.out[["proj"]][k], d.out[["proj_pct"]][k]))
}, seq_len(nrow(d.out)))
cat(strrep("-", 78), "\n")
cat(sprintf("max|%%|  copula = %.2f%%   projection = %.3g%%\n",
            max(abs(d.out[["copula_pct"]])), max(abs(d.out[["proj_pct"]]))))
cat("Wrote mc-binary-coeftable.csv\n")
# +------------------------------------------------------------------------+
# | End of 60b-binary-coeftable.r                                           |
# +------------------------------------------------------------------------+
