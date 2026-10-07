#' Drop the non-eta columns from a `$eta`/`$ranef` data frame
#'
#' `$eta` carries `ID`, and for a mixture model `nmObjGet.ranef` also merges in
#' a `mixnum` column.  Neither is an eta, so neither may reach an `etaMat`:
#' `foceiSetup_` compares the column count against the model's `neta` and stops
#' with "The etaMat must have the same number of ETAs (cols) as the model."
#' Same strip as `R/resid.R`'s residual path.
#'
#' @param eta `$eta` / `$ranef` data frame
#' @return the same data frame with `ID`/`mixnum`/`MIXEST` removed
#' @noRd
#' @author Matthew L. Fidler
.nmDropNonEtaCols <- function(eta) {
  .w <- which(names(eta) %in% c("ID", "mixnum", "MIXEST"))
  if (length(.w) > 0L) {
    return(eta[, -.w, drop = FALSE])
  }
  eta
}

#' One occasion variable's etas as the "theta" IOV rewrite expands them
#'
#' The rewrite gives each occasion parameter one unit-variance eta per
#' occasion, `rx.<parameter>.<occasion>`, scaled by the parameter's standard
#' deviation; `$iov` holds the scaled occasion deviations.  The rewrite refuses a
#' correlated occasion block, so each parameter's own standard deviation is its
#' whole scale.
#'
#' @param n occasion variable
#' @param iov the fit's `$iov`
#' @param omega the fit's `$omega`, a list by level
#' @return data frame, one column per parameter and occasion
#' @noRd
.nmIovThetaEtas <- function(n, iov, omega) {
  .dt <- data.table::as.data.table(iov[[n]])
  .frm <- stats::as.formula(paste0("ID ~ ", n))
  .sd <- sqrt(diag(omega[[n]]))
  .ret <- NULL
  for (.nr in names(.dt)[-(1:2)]) {
    .df <- as.data.frame(data.table::dcast(.dt, formula = .frm, value.var = .nr)[, -1])
    names(.df) <- paste0("rx.", .nr, ".", names(.df))
    # with no occasion variance every eta gives the same (zero) deviation
    .df <- if (.sd[[.nr]] > 0) .df / .sd[[.nr]] else .df * 0
    .ret <- if (is.null(.ret)) .df else cbind(.ret, .df)
  }
  .ret
}

#' @export
nmObjGet.etaMat <- function(x, ...) {
  .ui <- x[[1]]
  if (is.null(.ui$eta)) {
    return(NULL)
  }
  .eta <- as.matrix(.nmDropNonEtaCols(.ui$eta))
  if (is.null(.ui$iov)) {
    return(.eta)
  }
  # $eta leaves the occasion etas out; when the IOV rewrite estimated them,
  # etaObf has every eta of the expanded model on the model's scale
  .eo <- .ui$etaObf
  if (is.null(.ui$iovNative) && is.data.frame(.eo)) {
    .eo <- as.matrix(.nmDropNonEtaCols(.eo[, names(.eo) != "OBJI", drop = FALSE]))
    if (ncol(.eo) > ncol(.eta) && all(colnames(.eta) %in% colnames(.eo))) {
      return(.eo)
    }
  }
  # saem's own IOV handling ("twoLevel", "collapsed") estimates the occasion
  # deviations directly: give a refit the etas its "theta" rewrite expands to
  as.matrix(do.call(
    cbind,
    c(list(.eta), lapply(names(.ui$iov), .nmIovThetaEtas, iov = .ui$iov, omega = .ui$omega))
  ))
}
