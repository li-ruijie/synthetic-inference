#!/usr/bin/env r
# +------------------------------------------------------------------------+
# |  dp-gaussian-mechanism.r -- the one Gaussian-mechanism calibration      |
# +------------------------------------------------------------------------+
#
# Why this file exists.  This calibration must have exactly one definition.  Two
# copies of a constant drift apart as soon as one of them is corrected.
#
# The exact analytic Gaussian mechanism replaces the classical Dwork bound.  When
# a calibration exists in two places, correcting one and not re-running what the
# other produced leaves stale numbers calibrated to a bound already
# rejected in writing as invalid.
#
# **A script fixed is not a result refreshed.**  The way to make that failure
# unwriteable is to leave nowhere for a second copy to live.  So: one definition,
# one file, sourced by everything that needs it.  Do not re-derive it locally.
#
# dp.gm.multiplier(eps, delta) returns sigma / Delta, the noise multiplier per
# unit of L2 sensitivity, for an (eps, delta)-DP Gaussian mechanism.
#
#     sigma = Delta * dp.gm.multiplier(eps, delta)
#
# For eps < 1 the classical Dwork-Roth bound sqrt(2 log(1.25/delta))/eps applies
# (2014, Thm A.1).  That bound is valid only for eps < 1: at eps >= 1 it can
# under-noise (at eps = 10, delta = 1e-6 it certifies only delta = 1.9e-6, not
# 1e-6), so for eps >= 1 we use the exact analytic Gaussian mechanism of Balle &
# Wang (2018), equivalently the Gaussian-DP calibration of Dong, Roth & Su (2022),
# which is valid for all eps > 0.
#
# +------------------------------------------------------------------------+

# gdp.delta(eps, mu): the (eps, delta) profile of a Gaussian mechanism with
# multiplier mu = Delta / sigma.  Monotone increasing in mu.  The e^eps term is
# evaluated in log-space so it stays stable at large eps.
gdp.delta <- function(eps, mu) {
    pnorm(mu / 2 - eps / mu) -
        exp(eps + pnorm(-mu / 2 - eps / mu, log.p = TRUE))
}

dp.gm.multiplier <- function(eps, delta) {
    if(eps < 1)                                       # classical, valid eps < 1
        return(sqrt(2 * log(1.25 / delta)) / eps)
    mu <- uniroot(function(m) gdp.delta(eps, m) - delta,     # exact, all eps
                  interval = c(1e-8, 1e3), tol = 1e-12)[["root"]]
    1 / mu
}
