# Shared fixtures for the generalized ll() carry tests (#1004,
# test-focei-lincmt-carry-ll.R / -ll-fit.R): likelihood models embedding
# the linCmt() concentration, and their event table.

# The model functions are rxode2 ini()/model() DSL blocks, which lintr's
# object_usage cannot follow.
# nolint start: object_usage_linter.
.carryModLlNorm <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    eta.cl ~ 0.1
    add.sd <- 0.5
  })
  model({
    cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
    v <- exp(tv)
    cp <- linCmt()
    ll(err) ~ -0.5 * log(2 * pi) - log(add.sd) - 0.5 * ((DV - cp) / add.sd)^2
  })
}

.carryModLlProp <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    eta.cl ~ 0.1
    prop.sd <- 0.2
  })
  model({
    cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
    v <- exp(tv)
    cp <- linCmt()
    ll(err) ~ -0.5 * log(2 * pi) - log(prop.sd * cp) - 0.5 * ((DV - cp) / (prop.sd * cp))^2
  })
}

# count endpoint whose rate is driven by the concentration, with a second
# eta that only enters the likelihood outside the concentration
.carryModLlPois <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    tb <- 0.5
    slope <- 0.3
    eta.cl ~ 0.1
    eta.b ~ 0.1
  })
  model({
    cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
    v <- exp(tv)
    cp <- linCmt()
    lam <- exp(tb + eta.b + slope * cp)
    ll(err) ~ DV * log(lam) - lam - lgamma(DV + 1)
  })
}

# the same eta on the kernel and in the likelihood intercept
.carryModLlPoisShared <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    tb <- 0.5
    slope <- 0.3
    eta.cl ~ 0.1
  })
  model({
    cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
    v <- exp(tv)
    cp <- linCmt()
    lam <- exp(tb + eta.cl + slope * cp)
    ll(err) ~ DV * log(lam) - lam - lgamma(DV + 1)
  })
}

.carryModLlNoCov <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    eta.cl ~ 0.1
    add.sd <- 0.5
  })
  model({
    cl <- exp(tcl) * exp(eta.cl)
    v <- exp(tv)
    cp <- linCmt()
    ll(err) ~ -0.5 * log(2 * pi) - log(add.sd) - 0.5 * ((DV - cp) / add.sd)^2
  })
}
# nolint end

.carryLlEv <- function() {
  ev <- .carryEv() # nolint: object_usage_linter.
  ev$dv <- 0
  ev$dv[ev$evid == 0] <- c(3, 2.5, 2, 1.5)
  ev
}

.carryLlPars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.5,
                  `ETA[1]` = 0.3)
