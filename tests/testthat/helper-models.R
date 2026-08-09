# Common model definitions for testing
# This file is loaded before tests run (helper- prefix ensures early loading)
# These models are used across multiple test files to reduce duplication and
# improve test execution time by parsing models once

#' Basic One Compartment PK Model
#'
#' A standard one-compartment pharmacokinetic model with first-order absorption.
#' Parameters: tka (log Ka), tcl (log Cl), tv (log V)
#' Random effects: eta.ka, eta.cl, eta.v
#' Error model: additive (add.sd)
#'
#' This is the most commonly used model in the test suite, appearing in 15+ files.
#' Use this when you need a simple 1-compartment model without parameter labels.
#'
#' @return A model function suitable for nlmixr2
one.compartment <- function() {
  ini({
    tka <- 0.45; label("Ka")
    tcl <- 1; label("Cl")
    tv <- 3.45; label("V")
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d/dt(depot) = -ka * depot
    d/dt(center) = ka * depot - cl / v * center
    cp = center / v
    cp ~ add(add.sd)
  })
}

one.compartment <- one.compartment() # Pre-parse for faster tests

#' One Compartment Model with Lag Time (KA1Lode)
#'
#' A one-compartment model with absorption lag time.
#' Used with warfarin data in lag time testing.
#' Parameters: ltlag (log Tlag), lka (log ka), lcl (log Cl), lv (log V)
#' Random effects: eta.tlag, eta.ka, eta.cl, eta.v
#' Error model: proportional + additive
#'
#' @return A model function suitable for nlmixr2
one.compartment.with.lag <- function() {
  ini({
    ltlag <- log(0.2)  # log Tlag (h)
    lka  <- log(1.15)  # log ka (/h)
    lcl  <- log(0.135) # log Cl (L/h)
    lv   <- log(8)     # log V (L)
    prop.err <- 0.15   # proportional error (SD/mean)
    add.err  <- 0.6    # additive error (mg/L)
    eta.tlag ~ 0.5 # IIV tlag
    eta.ka ~ 0.5   # IIV ka
    eta.cl ~ 0.1   # IIV cl
    eta.v  ~ 0.1   # IIV v
  })
  model({
    tlag <- exp(ltlag + eta.tlag)
    ka <- exp(lka + eta.ka)
    cl <- exp(lcl + eta.cl)
    v  <- exp(lv + eta.v)
    d/dt(gut) <- -ka*gut
    d/dt(central) <- ka*gut - (cl/v)*central
    lag(gut) <- tlag
    cp <- central/v
    cp ~ prop(prop.err) + add(add.err)
  })
}

one.compartment.with.lag <- one.compartment.with.lag()

#' Two Compartment PK Model (cmt2)
#'
#' A two-compartment model with first-order absorption and inter-compartmental clearance.
#' Used with xgxr case1_pkpd dataset.
#' Parameters: lka (log Ka), lv (log Vc), lcl (log Cl), lq (log Q), lvp (log Vp)
#' Random effects: eta.ka, eta.v, eta.cl
#' Error model: log-normal
#'
#' @return A model function suitable for nlmixr2
two.compartment <- function() {
  ini({
    lka <- log(0.1) # log Ka
    lv <- log(10) # Log Vc
    lcl <- log(4) # Log Cl
    lq <- log(10) # log Q
    lvp <- log(20) # Log Vp
    eta.ka ~ 0.01
    eta.v ~ 0.1
    eta.cl ~ 0.1
    logn.sd = 10
  })
  model({
    ka <- exp(lka + eta.ka)
    cl <- exp(lcl + eta.cl)
    v <- exp(lv + eta.v)
    q <- exp(lq)
    vp <- exp(lvp)
    linCmt() ~ lnorm(logn.sd)
  })
}

two.compartment <- two.compartment()

#' Two-endpoint PK/PD model (direct inhibitory Emax on a fixed baseline)
#'
#' Deliberately well conditioned and well scaled, and the two endpoints have
#' residual SDs an order of magnitude apart (0.4 vs 4) so that anything which
#' pairs an endpoint with the wrong residual shows up immediately.  ka and e0 are
#' fixed so the fit stays small.
#'
#' @return A model function suitable for nlmixr2
pkpd.two.endpoint <- function() {
  ini({
    tcl <- -3.2
    tv <- -1
    tec50 <- 0.5
    eta.cl ~ 0.09
    eta.v ~ 0.09
    pk.sd <- 0.4
    pd.sd <- 4
  })
  model({
    ka <- exp(0.5)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    e0 <- 100
    ec50 <- exp(tec50)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    eff <- e0 * (1 - cp / (ec50 + cp))
    cp ~ add(pk.sd) | cp
    eff ~ add(pd.sd) | eff
  })
}

#' Simulate data for [pkpd.two.endpoint()] at its own ini() values.
#'
#' @param n number of subjects
#' @param seed RNG seed (both the etas and the residual noise)
#' @return a data frame with `cp` and `eff` observations
mkPkpdTwoEndpointData <- function(n = 20L, seed = 42L) {
  set.seed(seed)
  .m <- rxode2::rxode2({
    ka <- exp(0.5)
    cl <- exp(-3.2 + eta.cl)
    v <- exp(-1 + eta.v)
    e0 <- 100
    ec50 <- exp(0.5)
    d/dt(depot) <- -ka * depot
    d/dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    eff <- e0 * (1 - cp / (ec50 + cp))
  })
  .ev <- rxode2::et(amt = 100, cmt = "depot") |>
    rxode2::et(c(0.5, 1, 2, 4, 8, 12, 24)) |>
    rxode2::et(id = seq_len(n))
  .s <- rxode2::rxSolve(.m, .ev,
                        omega = lotri::lotri(eta.cl ~ 0.09, eta.v ~ 0.09),
                        returnType = "data.frame")
  .pk <- data.frame(id = .s$id, time = .s$time,
                    dv = .s$cp + stats::rnorm(nrow(.s), 0, 0.4),
                    dvid = "cp", amt = 0, evid = 0)
  .pd <- data.frame(id = .s$id, time = .s$time,
                    dv = .s$eff + stats::rnorm(nrow(.s), 0, 4),
                    dvid = "eff", amt = 0, evid = 0)
  .dose <- data.frame(id = seq_len(n), time = 0, dv = NA_real_,
                      dvid = "cp", amt = 100, evid = 1)
  .d <- rbind(.dose, .pk, .pd)
  .d[order(.d$id, .d$time, -.d$evid), ]
}
