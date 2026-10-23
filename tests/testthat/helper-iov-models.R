# A model with a subject-level eta (eta.v) and an occasion eta
# (iov.cl ~ 0.1 | occ), for the tests of what a finished IOV fit can do
.iovFitMod <- function() {
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
    linCmt() ~ add(add.sd)
  })
}

# .iovFitMod() with a subject eta on cl as well, so saem's own two-level and
# collapsed IOV handling take it (without one, saem uses iovMethod = "theta")
.iovFitMuMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.cl ~ 0.3
    eta.v ~ 0.1
    iov.cl ~ 0.1 | occ
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

# .iovFitMod() as ODEs
.iovFitOdeMod <- function() {
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

# theo_sd with two occasions: before and after 6 h
.iovFitData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$occ <- ifelse(.d$TIME < 6, 1L, 2L)
  .d
}

# .iovFitData() with subject 3's observations all NA, so the fit drops it
# and its table reinserts it
.iovFitNoObsData <- function() {
  .d <- .iovFitData()
  .d$DV[.d$ID == 3 & .d$EVID == 0] <- NA
  .d
}

# A quiet fit with no covariance step; `control` adds to it
.iovFit <- function(
  est,
  control = list(),
  table = tableControl(),
  model = .iovFitMod,
  data = .iovFitData()
) {
  suppressWarnings(suppressMessages(nlmixr2(
    model,
    data,
    est,
    control = c(list(print = 0, covMethod = ""), control),
    table = table
  )))
}

# Whether rxode2's rxSolve() turns the caller's nesting (occasion) column into
# a factor in place, as rxode2 5.1.8 does
.rxode2WritesNestingData <- function() {
  .d <- .iovFitData()
  .d <- .d[.d$ID %in% 1:2, ]
  .ref <- .iovFitData()
  .ref <- .ref[.ref$ID %in% 1:2, ]
  invisible(suppressWarnings(suppressMessages(rxode2::rxSolve(.iovFitMod, .d))))
  !identical(.d, .ref)
}
