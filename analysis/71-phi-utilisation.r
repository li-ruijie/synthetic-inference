# +------------------------------------------------------------------------+
# | 71-phi-utilisation.r -- PHI healthcare-utilisation application          |
# +------------------------------------------------------------------------+
#
# A two-part (logistic zero model + quasi-gamma log-link positive block) model
# of specialist-outpatient utilisation (util_soc, 74% zeros) on the PHI health
# index and demographics, the structure of the Yap et al. (2020, PLOS one)
# two-part spending model reproduced on the baseline survey.  The projection
# reproduces the analyst's two-part fit while the copula distorts the
# health-index effect.  Exercises the two-part projection, the
# "Proposition (two-part projection)" of S13, on its released-design path: the
# call below leaves design.exact at its FALSE default, so the full-synthesis
# corollary is not what runs here.  Also the non-canonical partial-moment
# bound (the partial moment matching theorem), on real data.
#
# The PHI dataset is restricted.  This script skips gracefully when it is absent.
#
# Input:  phi_phase1_baseline.dta (restricted)
# Output: phi-utilisation.csv (per-coefficient, two blocks x methods)
# Usage:  Rscript 71-phi-utilisation.r
# +------------------------------------------------------------------------+
.shared.loaded <- TRUE
.sobol.main <- TRUE
source("lib/00-shared.r")

if (!file.exists("phi_phase1_baseline.dta")) {
    cat("PHI data not available (restricted). Skipping.\n")
    q(save = "no")
}

source("lib/gen-syn-family.r")
source("lib/11-projection-generator.r")
rm(.sobol.main)
suppressWarnings(suppressMessages(library(haven)))

SEED <- 6489L

d.raw <- read_dta("phi_phase1_baseline.dta")
zap <- function(x) as.numeric(zap_labels(x))
d.app <- data.frame(
    response = zap(d.raw[["util_soc"]]),          # specialist-outpatient visits (count, ZI)
    phi      = zap(d.raw[["phi"]]),               # health index 0-100 (the key predictor)
    age      = 2016 - zap(d.raw[["birthyear"]]),
    female   = factor(zap(d.raw[["female"]])),
    race     = factor(zap(d.raw[["race"]])))
d.app <- d.app[complete.cases(d.app), ]

sup <- c(response = "zeroinfl", phi = "real", age = "real")

fit.zero <- function(d) {
    d[["z"]] <- as.integer(d[["response"]] > 0)
    glm(z ~ phi + age + female + race, d, family = binomial())
}
fit.pos <- function(d) {
    dp <- droplevels(d[d[["response"]] > 0, ])
    mm <- model.matrix(response ~ phi + age + female + race, dp)
    glm(response ~ phi + age + female + race, dp,
        family = quasi(link = "log", variance = "mu^2"),
        start = c(log(mean(dp[["response"]])), rep(0, ncol(mm) - 1L)),
        control = glm.control(maxit = 200L))
}

mz.o <- fit.zero(d.app)
mp.o <- fit.pos(d.app)

d.prj <- gen.syn.projected.twopart(d.app, y.nm = "response", family = "quasi-gamma",
                                   link = "log", supports = sup, seed = SEED)
d.cop <- gen.syn.family(d.app, order = 4, supports = sup, seed = SEED)

coef.row <- function(block, method, mo.s, mo.o) {
    vo <- coef(mo.o); vs <- coef(mo.s)
    cc <- intersect(names(vo), names(vs))
    data.frame(block = block, method = method, coefficient = cc,
               original = as.numeric(vo[cc]), estimate = as.numeric(vs[cc]),
               stringsAsFactors = FALSE)
}
d.out <- do.call(rbind, list(
    coef.row("zero",     "projection", fit.zero(d.prj), mz.o),
    coef.row("zero",     "copula",     fit.zero(d.cop), mz.o),
    coef.row("positive", "projection", fit.pos(d.prj),  mp.o),
    coef.row("positive", "copula",     fit.pos(d.cop),  mp.o)))
d.out[["pct.dev"]] <- ifelse(abs(d.out[["original"]]) < 1e-10, NA_real_,
    100 * (d.out[["estimate"]] - d.out[["original"]]) / abs(d.out[["original"]]))
write.csv(d.out, "results/phi-utilisation.csv", row.names = FALSE)
cat(sprintf("PHI utilisation (n=%d): written phi-utilisation.csv\n", nrow(d.app)))
