#' Within-subject IWRES autocorrelation of a fit table
#'
#' Pooled lag-1 Pearson correlation and Durbin-Watson statistic of the
#' individual weighted residuals, pairing successive observations of the same
#' subject (and endpoint, when the table has several) ordered by time.
#' Dosing, censored and non-finite IWRES rows are dropped first; subjects with
#' fewer than 2 valid residuals are skipped.
#'
#' @param df fit table (data frame) with at least `ID` and `IWRES`
#' @return `list(lag1R=, durbinWatson=, nSubjects=, nPairs=)`, or `NULL` when
#'   the table has no IWRES
#' @noRd
.iwresAutocor <- function(df) {
  if (!is.data.frame(df) || !all(c("ID", "IWRES") %in% names(df))) {
    return(NULL)
  }
  .keep <- is.finite(df$IWRES)
  if ("EVID" %in% names(df)) {
    .keep <- .keep & !is.na(df$EVID) & df$EVID == 0
  }
  if ("CENS" %in% names(df)) {
    .keep <- .keep & !is.na(df$CENS) & df$CENS == 0
  }
  .id <- as.integer(factor(df$ID[.keep]))
  .grp <- if ("CMT" %in% names(df)) as.integer(factor(df$CMT[.keep])) else rep(1L, length(.id))
  .time <- if ("TIME" %in% names(df)) df$TIME[.keep] else seq_along(.id)
  .iwres <- as.double(df$IWRES[.keep])
  .o <- order(.id, .grp, .time)
  .Call(`_nlmixr2est_iwresAutocor_`, .iwres[.o], .id[.o], .grp[.o])
}

#' @export
nmObjGetData.iwresAutocor <- function(x, ...) {
  .iwresAutocor(as.data.frame(x[[1]]))
}
attr(nmObjGetData.iwresAutocor, "desc") <- "within-subject IWRES lag-1 correlation and Durbin-Watson statistic"
attr(nmObjGetData.iwresAutocor, "rstudio") <- list(lag1R = 0, durbinWatson = 2, nSubjects = 1L, nPairs = 1L)

#' Run note for within-subject IWRES autocorrelation
#'
#' Warns (so it reaches `$runInfo`) when the pooled Durbin-Watson statistic
#' is below 1.5 or above 2.5, unless `tableControl(iwresAutocor = FALSE)`.
#' @param fit the finalized fit
#' @param table the fit's `tableControl()`
#' @return the `fit$iwresAutocor` list (invisibly), or `NULL`
#' @noRd
.iwresAutocorNote <- function(fit, table) {
  if (isFALSE(table$iwresAutocor) || !is.data.frame(fit)) {
    return(invisible(NULL))
  }
  .ac <- tryCatch(.iwresAutocor(as.data.frame(fit)), error = function(e) NULL)
  .dw <- .ac$durbinWatson
  if (length(.dw) == 1L && is.finite(.dw) && (.dw < 1.5 || .dw > 2.5)) {
    warning(
      sprintf(
        "IWRES %s autocorrelated within subjects: Durbin-Watson %.2f",
        if (.dw < 1.5) "positively" else "negatively",
        .dw
      ),
      call. = FALSE
    )
  }
  invisible(.ac)
}
