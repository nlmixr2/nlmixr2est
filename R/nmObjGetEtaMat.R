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

#' @export
nmObjGet.etaMat <- function(x, ...) {
  .ui <- x[[1]]
  if (is.null(.ui$eta)) {
    return(NULL)
  }
  .eta <- as.matrix(.nmDropNonEtaCols(.ui$eta))
  if (is.null(.ui$iov)) {
    .eta
  } else {
    # $iov holds each occasion eta on the natural scale; a refit expands an
    # uncorrelated level with unit-variance etas, so divide by its SD there
    .om <- tryCatch(.ui$ui$omega, error = function(e) NULL)
    .omegaMode <- identical(tryCatch(.ui$control$iovMethod, error = function(e) NULL), "omega")
    .iovSd <- function(n, d) {
      .m <- if (is.list(.om)) .om[[n]] else NULL
      if (.omegaMode || !is.matrix(.m) || !(d %in% rownames(.m)) || any(.m[upper.tri(.m)] != 0)) {
        return(1)
      }
      sqrt(.m[d, d])
    }
    .n <- names(.ui$iov)
    as.matrix(do.call(
      `cbind`,
      c(
        list(.eta),
        lapply(.n, function(n) {
          .dt <- data.table::as.data.table(.ui$iov[[n]])
          .frm <- eval(str2lang(paste0("ID ~ ", n)))
          .nr <- names(.dt)[-(1:2)]
          do.call(
            `cbind`,
            lapply(.nr, function(nr) {
              .dt0 <- .dt[, c("ID", n, nr)]
              .df <- as.data.frame(data.table::dcast(.dt, formula = .frm, value.var = nr)[, -1])
              names(.df) <- paste0("rx.", nr, ".", names(.df))
              .df / .iovSd(n, nr)
            })
          )
        })
      )
    ))
  }
}
