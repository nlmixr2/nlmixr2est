# An ODE model with a subject-level eta (eta.v) and an occasion eta
# (iov.cl ~ 0.1 | occ).  An ODE model, not linCmt(), so a solve on a finer
# time grid (augPred()) gives the fit's predictions at the observations:
# linCmt() takes an interval's parameters from its end time, so added grid
# points change its solution where a covariate changes.
.iovOdeMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.v ~ 0.1
    iov.cl ~ 0.1 | occ
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + iov.cl)
    v <- exp(tv + eta.v)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

# theo_sd with two occasions; "subject" puts the boundary at 3 h for even
# subjects and 8 h for odd ones, "same" at 6 h for all
.iovOdeData <- function(occ = c("same", "subject")) {
  occ <- match.arg(occ)
  .d <- nlmixr2data::theo_sd
  .cut <- if (occ == "same") 6 else ifelse(.d$ID %% 2 == 0, 3, 8)
  .d$occ <- ifelse(.d$TIME < .cut, 1L, 2L)
  .d
}

# A quiet fit with no covariance step
.iovOdeFit <- function(est, occ = "same") {
  suppressWarnings(suppressMessages(
    nlmixr2(.iovOdeMod, .iovOdeData(occ), est, control = list(print = 0, covMethod = ""))
  ))
}

# Whether rxode2's .rxParamsZero() zeroes the etas of every block of a nested
# (list) omega, as population predictions need
.rxode2ZeroesNestedOmega <- function() {
  .om <- lotri::lotri(lotri::lotri(eta.v ~ 0.1) | id(nu = 100), lotri::lotri(iov.cl ~ 0.1) | occ(nu = 100))
  all(c("eta.v", "iov.cl") %in% names(rxode2:::.rxParamsZero(c(tcl = 1), .om)))
}
