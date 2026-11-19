# Whether the loaded rxode2 advances a linCmt() interval with that interval's
# covariates (rxode2 PR 1517: locf the start record's, midpoint their mean).
# Without it every interpolation advances with the end record's, and the
# carry reads each covariate at the record that ends the interval
# (end-record mode).  Remove with test-dependency-workarounds.R's
# interval-covariate block.

.carryIntervalEnv <- new.env(parent = emptyenv())

#' Does the loaded rxode2 advance a linCmt() interval with that interval's
#' covariates?  Probed once per session.
#' @noRd
.rxFoceiLinCmtIntervalCovs <- function() {
  if (is.null(.carryIntervalEnv$covs)) {
    .carryIntervalEnv$covs <- isTRUE(try(.rxFoceiLinCmtIntervalProbe(), silent = TRUE))
  }
  .carryIntervalEnv$covs
}

#' One linCmt() solve: the interval (0, 1] runs at k = 1 when it is advanced
#' with the start record's covariate under locf, at the end record's k = 2
#' when it is not
#' @noRd
.rxFoceiLinCmtIntervalProbe <- function() {
  .m <- suppressMessages(rxode2::rxode2("cl <- rxLcProbeCov\nv <- 1\ncp <- linCmt()"))
  .ev <- data.frame(
    id = 1L,
    time = c(0, 1),
    amt = c(1, 0),
    evid = c(1L, 0L),
    cmt = 1L,
    rxLcProbeCov = c(1, 2)
  )
  .r <- rxode2::rxSolve(.m, events = .ev, covsInterpolation = "locf", returnType = "data.frame")
  abs(.r$cp[1] - exp(-1)) < 1e-8
}
