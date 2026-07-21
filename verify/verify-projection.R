#!/usr/bin/env r
# Assertion-based verification of the design-Gram projection.  The
# strata.exact permutation draw of gen-syn-family.r, project.designgram and the
# design.exact generator paths of 11-projection-generator.r. Six tasks. Every
# task stops the run on failure. Companion MC evidence is 69-design-gram.r.
.sobol.main <- TRUE
source("lib/00-shared.r")
source("lib/gen-syn-family.r")
source("lib/07-calibrated-generator.r")
source("lib/09-symmetric-generator.r")
source("lib/11-projection-generator.r")
rm(.sobol.main)

set.seed(6489)
d.t <- data.frame(
    response = rnorm(400, 10, 2),
    age      = rnorm(400, 45, 9),
    income   = rlnorm(400, 3, 0.4),
    group    = factor(sample(c("a", "b", "c"), 400, TRUE)),
    area     = factor(sample(c("u", "r"), 400, TRUE)))

# --- Task 1: strata.exact pins every cell count -----------------------------
d.s1 <- gen.syn.family(d.t, order = 4, shrink = TRUE, seed = 101,
                       strata.exact = TRUE)
stopifnot(identical(table(d.s1[["group"]], d.s1[["area"]]),
                    table(d.t[["group"]], d.t[["area"]])))
# the permutation must still randomise the order (a fixed identity would leak
# row alignment): at least one label differs positionally
stopifnot(any(as.character(d.s1[["group"]]) != as.character(d.t[["group"]])))
cat("verify-projection: Task 1 strata.exact cell counts OK\n")

# --- Task 2: OLS design-Gram identity and the exact triple ------------------
d.g <- gen.syn.projected(d.t, y.nm = "response", family = "gaussian",
                         link = "identity", seed = 303, design.exact = TRUE)
fml.t <- response ~ age + income + group + area
m.X.o <- model.matrix(fml.t, d.t)
m.X.s <- model.matrix(fml.t, d.g)
s.scl <- max(abs(crossprod(m.X.o)))
stopifnot(max(abs(crossprod(m.X.s) - crossprod(m.X.o))) < 1e-8 * s.scl)
mo.o <- lm(fml.t, d.t)
mo.s <- lm(fml.t, d.g)
stopifnot(max(abs(coef(mo.s) - coef(mo.o))) < 1e-8,
          max(abs(coef(summary(mo.s))[, 2] - coef(summary(mo.o))[, 2])) < 1e-8,
          abs(summary(mo.s)[["r.squared"]] - summary(mo.o)[["r.squared"]]) < 1e-8,
          abs(AIC(mo.s) - AIC(mo.o)) < 1e-6)
# submodel spot-check: the triple is exact on a strict submodel too
mo.o2 <- lm(response ~ age + group, d.t)
mo.s2 <- lm(response ~ age + group, d.g)
stopifnot(max(abs(coef(mo.s2) - coef(mo.o2))) < 1e-8,
          max(abs(coef(summary(mo.s2))[, 2] - coef(summary(mo.o2))[, 2])) < 1e-8)
cat("verify-projection: Task 2 OLS design-Gram triple OK\n")

# --- Task 3: GLM Fisher-weighted Gram within tolerance ----------------------
set.seed(6489)
v.age.g <- rnorm(400, 45, 9)
d.tp <- data.frame(
    response = pmax(rgamma(400, shape = 2,
                           rate = 2 / exp(1 + 0.02 * v.age.g)), 1e-3),
    age      = v.age.g,
    group    = factor(sample(c("a", "b", "c"), 400, TRUE)))
d.gg <- gen.syn.projected(d.tp, y.nm = "response", family = "gamma",
                          link = "log", seed = 404, design.exact = TRUE,
                          supports = c(response = "positive", age = "real"))
fml.g  <- response ~ age + group
# The release's exact score root is the original fit as computed by
# project.glm's internal call (default epsilon, maxit 200), so the reference
# fit uses that same control, and the refit is tightened to 1e-12 so the
# assertion measures the construction rather than the refit's own
# deviance-convergence slack (~2e-5 on coefficients at the default 1e-8).
mo.go  <- glm(fml.g, d.tp, family = Gamma(link = "log"),
              control = glm.control(maxit = 200))
m.Xo   <- model.matrix(fml.g, d.tp)
m.Xs   <- model.matrix(fml.g, d.gg)
fam.g  <- Gamma(link = "log")
w.info.chk <- function(m.X) {
    v.eta <- as.numeric(m.X %*% coef(mo.go))
    fam.g[["mu.eta"]](v.eta)^2 /
        fam.g[["variance"]](fam.g[["linkinv"]](v.eta))
}
m.I.o <- crossprod(m.Xo * sqrt(w.info.chk(m.Xo)))
m.I.s <- crossprod(m.Xs * sqrt(w.info.chk(m.Xs)))
stopifnot(max(abs(m.I.s - m.I.o)) < 1e-6 * max(abs(m.I.o)))
mo.gs <- glm(fml.g, d.gg, family = Gamma(link = "log"),
             control = glm.control(epsilon = 1e-12, maxit = 200))
stopifnot(max(abs(coef(mo.gs) - coef(mo.go))) < 1e-6,
          max(abs(coef(summary(mo.gs))[, 2] - coef(summary(mo.go))[, 2])
              / coef(summary(mo.go))[, 2]) < 1e-4)
cat("verify-projection: Task 3 GLM weighted Gram OK\n")

# --- Task 4: two-part full-synthesis corollary ------------------------------
set.seed(6489)
v.age.z <- rnorm(400, 45, 9)
v.eta.zi <- 0.8 + 0.03 * v.age.z
d.zi <- data.frame(
    response = ifelse(runif(400) < 0.35, 0,
                      rgamma(400, 2, 2 / exp(v.eta.zi))),
    age      = v.age.z,
    group    = factor(sample(c("a", "b"), 400, TRUE)))
d.z <- gen.syn.projected.twopart(d.zi, y.nm = "response",
                                 family = "quasi-gamma", link = "log",
                                 seed = 505, design.exact = TRUE,
                                 supports = c(response = "zeroinfl",
                                              age = "real"))
# all-rows and positive-rows Gram identities (the disjoint-block pins)
fml.z  <- ~ age + group
m.Xo.z <- model.matrix(fml.z, d.zi)
m.Xs.z <- model.matrix(fml.z, d.z)
s.scl.z <- max(abs(crossprod(m.Xo.z)))
stopifnot(max(abs(crossprod(m.Xs.z) - crossprod(m.Xo.z))) < 1e-8 * s.scl.z)
m.Xo.p <- m.Xo.z[d.zi[["response"]] > 0, , drop = FALSE]
m.Xs.p <- m.Xs.z[d.z[["response"]] > 0, , drop = FALSE]
stopifnot(max(abs(crossprod(m.Xs.p) - crossprod(m.Xo.p))) < 1e-8 * s.scl.z)
# zero block: coefficients near the original zero-model fit
mo.z0 <- glm(I(response > 0) ~ age + group, d.zi, family = binomial())
mo.zs <- glm(I(response > 0) ~ age + group, d.z,  family = binomial())
stopifnot(max(abs(coef(mo.zs) - coef(mo.z0))) < 1e-3)
# positive block: the exact triple on the pinned positive design (reference
# fit at project.glm's own control, refit tightened, per the Task 3 lesson)
d.zi.p <- d.zi[d.zi[["response"]] > 0, ]
d.z.p  <- d.z[d.z[["response"]] > 0, ]
mo.p0 <- glm(response ~ age + group, d.zi.p,
             family = quasi(link = "log", variance = "mu^2"),
             start = c(1, 0, 0), control = glm.control(maxit = 200))
mo.ps <- glm(response ~ age + group, d.z.p,
             family = quasi(link = "log", variance = "mu^2"),
             start = c(1, 0, 0),
             control = glm.control(epsilon = 1e-12, maxit = 200))
stopifnot(max(abs(coef(mo.ps) - coef(mo.p0))) < 1e-6,
          max(abs(coef(summary(mo.ps))[, 2] - coef(summary(mo.p0))[, 2])
              / coef(summary(mo.p0))[, 2]) < 1e-4)
cat("verify-projection: Task 4 two-part corollary OK\n")

# --- Task 5: engineered failures stop loudly --------------------------------
# (a) precondition: a bootstrap draw (multinomial counts) must be refused
d.boot <- gen.syn.family(d.t, order = 4, shrink = TRUE, seed = 606)
b.err.a <- tryCatch({
    project.designgram(d.boot, d.t, "response",
                       setdiff(names(d.t), "response"))
    FALSE
}, error = function(e) grepl("D'D not pinned", conditionMessage(e)))
stopifnot(isTRUE(b.err.a))
# (b) degenerate synthetic residual: a continuous column equal to the
# fitted values of its additive dummy regression lies in span(D) exactly,
# so the residual Gram is singular
d.dgn <- d.t
d.dgn[["age"]] <- fitted(lm(age ~ group + area, d.t))
b.err.b <- tryCatch({
    project.designgram(d.dgn, d.t, "response",
                       setdiff(names(d.t), "response"))
    FALSE
}, error = function(e) grepl("design-gram", conditionMessage(e)))
stopifnot(isTRUE(b.err.b))
cat("verify-projection: Task 5 loud failure OK\n")

# --- Task 6: byte-identity of the default path ------------------------------
# The default strata.exact = FALSE output must be unchanged.
d.def.a <- gen.syn.family(d.t, order = 4, shrink = TRUE, seed = 202)
d.def.b <- gen.syn.family(d.t, order = 4, shrink = TRUE, seed = 202,
                          strata.exact = FALSE)
stopifnot(identical(d.def.a, d.def.b))
cat("verify-projection: Task 6 default byte-identity OK\n")

cat("verify-projection: ALL TASKS OK\n")
