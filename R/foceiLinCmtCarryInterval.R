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

#' Does the loaded rxode2 impute an alag() dose's covariates at its arrival?
#'
#' rxode2 PR 1635 (#1584): under `addlKeepsCov = FALSE` a dose moved by
#' alag() reads the covariates interpolated at its arrival over the data
#' records' times, where its own record's time stays a knot.  Probed once per
#' session.
#' @noRd
.rxFoceiLinCmtLagImputes <- function() {
  if (is.null(.carryIntervalEnv$lag)) {
    .carryIntervalEnv$lag <- isTRUE(try(.rxFoceiLinCmtLagProbe(), silent = TRUE))
  }
  .carryIntervalEnv$lag
}

#' f(depot) reads the covariate at the dose's arrival (2 under nocb) when
#' the lagged dose is imputed, its record's (1) when it is not
#' @noRd
.rxFoceiLinCmtLagProbe <- function() {
  .m <- suppressMessages(rxode2::rxode2(paste(
    "cl <- 0.1",
    "v <- 1",
    "ka <- 1",
    "f(depot) <- rxLcProbeCov",
    "alag(depot) <- 1",
    "cp <- linCmt()",
    sep = "\n"
  )))
  .ev <- data.frame(
    id = 1L,
    time = c(0, 10),
    amt = c(1, 0),
    evid = c(1L, 0L),
    cmt = 1L,
    rxLcProbeCov = c(1, 2)
  )
  .s <- function(keep) {
    .r <- rxode2::rxSolve(
      .m,
      events = .ev,
      covsInterpolation = "nocb",
      addlKeepsCov = keep,
      returnType = "data.frame"
    )
    .r$cp[nrow(.r)]
  }
  abs(.s(FALSE) / .s(TRUE) - 2) < 1e-6
}

#' Does the carry read each interval's covariates around lagged-dose knots?
#'
#' True for a linCmt() carry build of a model with alag() whose rxode2
#' imputes a lagged dose's covariates and whose control keeps
#' `addlKeepsCov = FALSE`.  The interval that ends at a lagged dose then
#' holds its record's data time as a knot (`.rxFoceiLinCmtCarryCovInterval()`).
#'
#' @param ui rxode2 ui with its control
#' @return logical
#' @noRd
.rxFoceiLinCmtCarryLagKnot <- function(ui) {
  .rx <- rxode2::rxGetControl(ui, "rxControl", NULL)
  !isTRUE(.rx$addlKeepsCov) &&
    .rxFoceiLinCmtCarryHasLag(ui) &&
    .rxFoceiLinCmtIntervalCovs() &&
    .rxFoceiLinCmtLagImputes()
}

#' Does the model declare an alag()/lag() on any compartment?
#' @noRd
.rxFoceiLinCmtCarryHasLag <- function(ui) {
  .txt <- vapply(ui$lstExpr, function(e) paste(deparse(e), collapse = " "), character(1))
  any(grepl("^\\s*(alag|lag)\\([^)]*\\)\\s*(<-|=)", .txt))
}
