#!/usr/bin/env r
# +------------------------------------------------------------+
# |  Baseline comparison: DP-Projection vs Gaussian-on-X'y     |
# +------------------------------------------------------------+
#
# The paper's Supplementary Material S6 analyses the Gaussian
# mechanism applied to the sufficient statistic X'y (input
# perturbation):
#
#   (X'y)_DP = X'y + N(0, sigma^2 I_p)
#   beta_DP = (X'X)^{-1} (X'y)_DP
#
# S6 is analytical (ridge-bias formula).  This script runs an empirical
# Monte Carlo for the same mechanism and compares decision agreement
# against DP-Projection (output perturbation on beta).  Direct
# apples-to-apples comparison at matched (eps, delta).
#
# +------------------------------------------------------------+
.sobol.main <- TRUE
source("lib/00-shared.r")
source("lib/gen-syn-family.r")
source("lib/07-calibrated-generator.r")
source("lib/09-symmetric-generator.r")
source("lib/11-projection-generator.r")
source("lib/28-dp-projection.r")
source("lib/dp-public-bounds.r")            # the ONLY permitted source of a DP bound
rm(.sobol.main)

rlaplace <- function(n, location = 0, scale = 1) {
    v.u <- runif(n) - 0.5
    location - scale * sign(v.u) * log(1 - 2 * abs(v.u))
}

# +------------------------------------------------------------+
# |  Gaussian-on-X'y baseline                                  |
# +------------------------------------------------------------+
# Sensitivity of X'y under one-observation replacement:
#   || (x_i' y_i) - (x_i'' y_i'') ||_2 <= 2 * B_x * B_y (worst case)
# where B_x = max ||x_i||_2 (public bound).
#
# For our scenarios, X includes intercept + age + female + race + income.
# Typical ||x_i|| is dominated by income (max ~3778 for Sc.1), so
# sensitivity is large.
beta.via.xy.perturbation <- function(d, y.nm = "response", eps = 1,
                                       delta = 1e-6, B.x, B.y, seed = 42) {
    fml <- reformulate(setdiff(names(d), y.nm), response = y.nm)
    m.X <- model.matrix(fml, d)
    v.y <- d[[y.nm]]
    p <- ncol(m.X)
    v.Xty <- as.numeric(crossprod(m.X, v.y))
    Delta <- 2 * B.x * B.y
    sigma <- Delta * dp.gm.multiplier(eps, delta)
    set.seed(seed)
    v.noise <- rnorm(p, 0, sigma)
    v.Xty.dp <- v.Xty + v.noise
    m.XtX.inv <- solve(crossprod(m.X))
    v.beta.dp <- as.numeric(m.XtX.inv %*% v.Xty.dp)
    names(v.beta.dp) <- colnames(m.X)
    # Simulate SE from OLS-on-synthetic approach: compute residual
    # variance from predicted - actual using beta_DP (approximation).
    v.resid <- v.y - m.X %*% v.beta.dp
    s.sigma2 <- sum(v.resid^2) / (length(v.y) - p)
    v.se <- sqrt(s.sigma2 * diag(m.XtX.inv))
    list(beta.dp = v.beta.dp, se = v.se, sigma.xy = sigma)
}

decision.agreement.from.coefs <- function(b.real, se.real, p.real,
                                             b.syn, se.syn) {
    v.b.real <- b.real; v.p.real <- p.real
    v.b.syn  <- b.syn
    v.se.syn <- se.syn
    # approximate p-value for synthetic: 2*(1-Phi(|b|/se))
    v.t <- v.b.syn / v.se.syn
    v.p.syn <- 2 * (1 - pnorm(abs(v.t)))
    b.sign <- sign(v.b.real) == sign(v.b.syn)
    b.sig <- (v.p.real < 0.05) == (v.p.syn < 0.05)
    v.lo <- v.b.syn - 1.96 * v.se.syn
    v.hi <- v.b.syn + 1.96 * v.se.syn
    b.ci <- (v.b.real >= v.lo) & (v.b.real <= v.hi)
    list(sign.agree = mean(b.sign), sig.agree = mean(b.sig),
         both.agree = mean(b.sign & b.sig), ci.overlap = mean(b.ci))
}

# +------------------------------------------------------------+
# |  Run MC: Gaussian-on-X'y vs DP-Projection V1               |
# +------------------------------------------------------------+
load("results/d.sim.scenarios.RData")

s.B <- 15L
v.eps <- c(1, 2, 5, 10)
s.delta <- 1e-6

cat("BASELINE COMPARISON: Gaussian-on-X'y vs DP-Projection V1\n")
cat(strrep("=", 78), "\n\n")

l.results <- list()
for(s.sc in c(1, 2)) {
    d <- l.scenarios[[s.sc]][["data"]]
    d[["female"]] <- factor(as.character(d[["female"]]), levels = c("0","1"))
    d[["race"]] <- factor(as.character(d[["race"]]), levels = c("1","2","3"))
    fml <- response ~ age + female + race + income

    # Both bounds come from dp-public-bounds.r.  Taking them from this sample, as
    # in B.x <- sqrt(max(rowSums(m.X^2))) or B.y <- 1.5 * max(abs(response)),
    # would make Delta = 2 B.x B.y data-derived and void the epsilon reported
    # below.  The paper's central negative result is read against this baseline.
    # The baseline is generated
    # from, and scored against, the publicly clipped data its bounds are bounds
    # for -- the same discipline Route 1 already follows.
    d <- dp.public.clip.x(dp.public.clip.y(d, "response", s.sc))
    l.Bx <- dp.public.B.x()
    l.By <- dp.public.B.y(s.sc)
    B.x  <- l.Bx[["B.x"]]
    B.y  <- l.By[["B.y"]]

    mo.real <- lm(fml, d)
    v.beta.real <- coef(mo.real)
    v.se.real <- summary(mo.real)[["coefficients"]][, 2]
    v.p.real  <- summary(mo.real)[["coefficients"]][, 4]
    sup <- c(response = "real", age = "real", income = "real")

    cat(sprintf("=== Scenario %d ===\n", s.sc))
    dp.bounds.announce(B.x = l.Bx, B.y = l.By)
    cat(sprintf("B_x = %.1f, B_y = %.1f, sigma_Xy scale = %.1e\n\n",
                B.x, B.y, 2 * B.x * B.y / 10))

    cat(sprintf("%-10s  %-6s  %-8s  %-8s  %-8s  %-8s\n",
                "method", "eps", "sign", "sig", "both", "ci"))
    cat(strrep("-", 60), "\n")

    for(eps in v.eps) {
        # Gaussian-on-X'y
        v.seeds <- 6489L + seq_len(s.B)
        l.xy <- Map(function(sd) {
            res <- beta.via.xy.perturbation(d, "response", eps = eps,
                delta = s.delta, B.x = B.x, B.y = B.y, seed = sd)
            decision.agreement.from.coefs(v.beta.real, v.se.real,
                v.p.real, res[["beta.dp"]], res[["se"]])
        }, v.seeds)
        v.xy.sign <- sapply(l.xy, function(r) r[["sign.agree"]])
        v.xy.sig  <- sapply(l.xy, function(r) r[["sig.agree"]])
        v.xy.both <- sapply(l.xy, function(r) r[["both.agree"]])
        v.xy.ci   <- sapply(l.xy, function(r) r[["ci.overlap"]])
        cat(sprintf("%-10s  %-6g  %6.1f%%  %6.1f%%  %6.1f%%  %6.1f%%\n",
                    "Xy (S6)", eps,
                    100 * mean(v.xy.sign), 100 * mean(v.xy.sig),
                    100 * mean(v.xy.both), 100 * mean(v.xy.ci)))

        # DP-Projection V1
        l.v1 <- Map(function(sd) {
            res <- tryCatch(gen.syn.dp.projected(d, y.nm = "response",
                family = "gaussian", link = "identity", supports = sup,
                eps = eps, delta = s.delta, seed = sd),
                error = function(e) NULL)
            if(is.null(res)) return(NULL)
            mo.s <- tryCatch(lm(fml, res[["d.syn"]]),
                             error = function(e) NULL)
            if(is.null(mo.s)) return(NULL)
            s.s <- summary(mo.s)[["coefficients"]]
            v.nms <- intersect(names(v.beta.real), rownames(s.s))
            decision.agreement.from.coefs(
                v.beta.real[v.nms], v.se.real[v.nms], v.p.real[v.nms],
                s.s[v.nms, 1], s.s[v.nms, 2])
        }, v.seeds)
        l.v1 <- Filter(Negate(is.null), l.v1)
        if(length(l.v1) > 0) {
            v.v1.sign <- sapply(l.v1, function(r) r[["sign.agree"]])
            v.v1.sig  <- sapply(l.v1, function(r) r[["sig.agree"]])
            v.v1.both <- sapply(l.v1, function(r) r[["both.agree"]])
            v.v1.ci   <- sapply(l.v1, function(r) r[["ci.overlap"]])
            cat(sprintf("%-10s  %-6g  %6.1f%%  %6.1f%%  %6.1f%%  %6.1f%%\n",
                        "DP-Proj V1", eps,
                        100 * mean(v.v1.sign), 100 * mean(v.v1.sig),
                        100 * mean(v.v1.both), 100 * mean(v.v1.ci)))
        }

        l.results[[length(l.results) + 1L]] <- data.frame(
            scenario = s.sc, eps = eps,
            xy.sign = mean(v.xy.sign), xy.sig = mean(v.xy.sig),
            xy.both = mean(v.xy.both), xy.ci = mean(v.xy.ci),
            v1.sign = `if`(length(l.v1) > 0, mean(v.v1.sign), NA),
            v1.sig  = `if`(length(l.v1) > 0, mean(v.v1.sig), NA),
            v1.both = `if`(length(l.v1) > 0, mean(v.v1.both), NA),
            v1.ci   = `if`(length(l.v1) > 0, mean(v.v1.ci), NA))
    }
    cat("\n")
}

d.out <- do.call(rbind, l.results)
# output is written to the repository root
write.csv(d.out, "results/mc-baseline-comparison.csv", row.names = FALSE)

cat(strrep("=", 78), "\n")
cat("Results in mc-baseline-comparison.csv\n")
