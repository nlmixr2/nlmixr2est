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

#' One occasion variable's etas as a refit expands them
#'
#' A refit expands an uncorrelated level with one unit-variance eta per
#' occasion, `rx.<parameter>.<occasion>`, scaled by the parameter's standard
#' deviation; `$iov` holds the scaled occasion deviations.  A correlated block,
#' or `iovMethod = "omega"`, expands the deviations unscaled.
#'
#' @param n occasion variable
#' @param iov the fit's `$iov`
#' @param omega the fit's `$omega`, a list by level
#' @param unscaled `TRUE` when the refit expands the level unscaled
#' @return data frame, one column per parameter and occasion
#' @noRd
.nmIovThetaEtas <- function(n, iov, omega, unscaled = FALSE) {
  .dt <- data.table::as.data.table(iov[[n]])
  .frm <- stats::as.formula(paste0("ID ~ ", n))
  .m <- if (is.list(omega)) omega[[n]] else NULL
  .scale <- !unscaled && is.matrix(.m) && all(.m[upper.tri(.m)] == 0)
  .ret <- NULL
  for (.nr in names(.dt)[-(1:2)]) {
    .df <- as.data.frame(data.table::dcast(.dt, formula = .frm, value.var = .nr)[, -1])
    names(.df) <- paste0("rx.", .nr, ".", names(.df))
    if (.scale && .nr %in% rownames(.m)) {
      .sd <- sqrt(.m[.nr, .nr])
      # with no occasion variance every eta gives the same (zero) deviation
      .df <- if (.sd > 0) .df / .sd else .df * 0
    }
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
  # $iov holds each occasion eta on the natural scale; rescale it the way a
  # refit expands the level
  .om <- tryCatch(.ui$ui$omega, error = function(e) NULL)
  .omegaMode <- identical(tryCatch(.ui$control$iovMethod, error = function(e) NULL), "omega")
  .ret <- as.matrix(do.call(
    cbind,
    c(list(.eta), lapply(names(.ui$iov), .nmIovThetaEtas, iov = .ui$iov, omega = .om, unscaled = .omegaMode))
  ))
  # a correlated occasion block is expanded occasion by occasion; etaObf
  # carries the expanded model's eta order
  .eo <- tryCatch(names(.ui$etaObf), error = function(e) NULL)
  if (all(colnames(.ret) %in% .eo)) {
    .ret <- .ret[, intersect(.eo, colnames(.ret)), drop = FALSE]
  }
  .ret
}
