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
    tka <- 0.45 # Log Ka
    tcl <- 1 # Log Cl
    tv <- 3.45    # Log V
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

#' One Compartment Model (add.err variant)
#'
#' Same as one.compartment() but uses add.err instead of add.sd for the error parameter.
#' Used in broom tests where the parameter name is checked explicitly.
#'
#' @return A model function suitable for nlmixr2
one.compartment.add.err <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.err <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.err)
  })
}

one.compartment.add.err <- one.compartment.add.err()

#' One Compartment Model with Parameter Labels
#'
#' Same as one.compartment() but with human-readable labels on parameters.
#' Used in tests that check label handling or require labeled output.
#'
#' @return A model function suitable for nlmixr2
one.compartment.with.labels <- function() {
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

one.compartment.with.labels <- one.compartment.with.labels()

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
