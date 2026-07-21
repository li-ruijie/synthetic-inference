# +------------------------------------------------------------------------+
# | 76-misspec-identity.r -- the analyst's misspecified fit on the release  |
# +------------------------------------------------------------------------+
#
# Why this script exists.
#
# The paper claims the projection reproduces the analyst's fitted
# coefficient whether or not the model is correctly specified, the
# assumption-lean reading.  This script fits deliberately misspecified
# analyst models on the original data and on the release and reports the
# maximum absolute percentage coefficient deviation between the two fits.
# The misspecifications: an omitted covariate (income dropped, a submodel,
# so part (a) of the projection theorem predicts identity), a wrong
# canonical family (Gaussian least squares on the count release, the
# cross-family proposition predicting identity when the covariates are
# public), and a wrong non-canonical link (probit on the logit release, and
# a log-link Gaussian analyst on the gamma release, where no exactness is
# claimed and the deviation is documented).  "fixed" releases the original
# covariates with the projected response, "full" the fully synthetic table.
#
# Output:  results/misspec-identity.csv
#          results/tab-supp-misspec-identity.tex
# Usage:   Rscript analysis/76-misspec-identity.r
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

load("results/d.sim.scenarios.RData")

releases <- function(s.i, sup) {
    s   <- l.scenarios[[s.i]]
    d   <- s[["data"]]
    fam <- s[["family"]]
    lnk <- s[["link"]]
    set.seed(SEED)
    d.cop <- gen.syn.family(d, order = 4, supports = sup, seed = SEED)
    d.mix <- `[[<-`(d, "response", value = d.cop[["response"]])
    d.fix <- `[[<-`(d, "response",
                    value = project.glm(d.mix, d, "response", V.X, fam, lnk))
    set.seed(SEED)
    d.full <- gen.syn.projected(d, y.nm = "response", family = fam, link = lnk,
                                supports = sup, seed = SEED)
    list(orig = d, fix = d.fix, full = d.full)
}

maxpct <- function(a, b) max(abs(100 * (a - b) / pmax(abs(b), 1e-10)))

one.case <- function(s.case, l.d, fit) {
    b.o <- coef(fit(l.d[["orig"]]))
    b.x <- coef(fit(l.d[["fix"]]))[names(b.o)]
    b.f <- coef(fit(l.d[["full"]]))[names(b.o)]
    row <- data.frame(case = s.case,
                      dev.fixed = maxpct(b.x, b.o),
                      dev.full  = maxpct(b.f, b.o))
    cat(sprintf("%-52s fixed %.3g%%  full %.3g%%\n",
                s.case, row[["dev.fixed"]], row[["dev.full"]]))
    row
}

sup1  <- c(response = "real",     age = "real", income = "real")
sup3  <- c(response = "positive", age = "real", income = "positive")
sup6  <- c(response = "count",    age = "real", income = "positive")
sup11 <- c(age = "real", income = "positive")

l1  <- releases(1L,  sup1)
l3  <- releases(3L,  sup3)
l6  <- releases(6L,  sup6)
l11 <- releases(11L, sup11)

fml.omit <- response ~ age + female + race
d.out <- rbind(
    one.case("Sc.1 omit income (Gaussian)", l1,
             function(d) glm(fml.omit, d, family = gaussian())),
    one.case("Sc.6 omit income (Poisson)", l6,
             function(d) suppressWarnings(glm(fml.omit, d, family = poisson()))),
    one.case("Sc.11 omit income (logistic)", l11,
             function(d) suppressWarnings(glm(fml.omit, d, family = binomial()))),
    one.case("Sc.6 Gaussian analyst on the count release", l6,
             function(d) glm(response ~ age + female + race + income, d,
                             family = gaussian())),
    one.case("Sc.11 probit analyst on the logit release", l11,
             function(d) suppressWarnings(
                 glm(response ~ age + female + race + income, d,
                     family = binomial(link = "probit")))),
    one.case("Sc.3 log-link Gaussian analyst on the gamma release", l3,
             function(d) glm(response ~ age + female + race + income, d,
                             family = gaussian(link = "log"),
                             start = c(log(mean(d[["response"]])), rep(0, 5)),
                             control = glm.control(maxit = 100L))))
write.csv(d.out, "results/misspec-identity.csv", row.names = FALSE)

# emit the supplement fragment from the CSV just written
d.emit <- read.csv("results/misspec-identity.csv")
dev.fmt <- function(x) ifelse(x < 1e-8, "$<10^{-8}$", sprintf("$%.2g$", x)) |>
    gsub("e-0?(\\d+)", " \\\\times 10^{-\\1}", x = _) |>
    gsub("e\\+0?(\\d+)", " \\\\times 10^{\\1}", x = _)
row.tex <- function(k) {
    r <- d.emit[k, ]
    sprintf("        %s & %s & %s \\\\", r[["case"]],
            dev.fmt(r[["dev.fixed"]]), dev.fmt(r[["dev.full"]]))
}
v.lines <- c(
    "% Emitted by analysis/76-misspec-identity.r from results/misspec-identity.csv.",
    "% Do not edit by hand; rerun the script.",
    "\\begin{table}[htbp]",
    "    \\centering",
    "    \\caption{Misspecified analyst fits on the release: maximum absolute",
    "        percentage coefficient deviation of the release refit from the",
    "        same misspecified fit on the original data.  Fixed, the original",
    "        covariates with the projected response; full, the fully synthetic",
    "        table.}",
    "    \\label{tab:supp-misspec-identity}",
    "    \\small",
    "    \\begin{tabular}{@{}lrr@{}}",
    "        \\toprule",
    "        Misspecified analyst model & Fixed (\\%) & Full (\\%) \\\\",
    "        \\midrule",
    unlist(Map(row.tex, seq_len(nrow(d.emit)))),
    "        \\bottomrule",
    "    \\end{tabular}",
    "\\end{table}")
writeLines(v.lines, "results/tab-supp-misspec-identity.tex")
cat("Written results/misspec-identity.csv and tab-supp-misspec-identity.tex\n")
