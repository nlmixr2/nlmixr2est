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

# theo_sd with two occasions: before and after 6 h
.iovFitData <- function() {
  .d <- nlmixr2data::theo_sd
  .d$occ <- ifelse(.d$TIME < 6, 1L, 2L)
  .d
}

# A quiet fit of .iovFitMod() with no covariance step; `control` adds to it
.iovFit <- function(est, control = list(), table = tableControl()) {
  suppressWarnings(suppressMessages(nlmixr2(
    .iovFitMod,
    .iovFitData(),
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
