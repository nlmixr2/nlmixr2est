#' Parameters that reach a delay() duration
#'
#' The symbolic sensitivities miss the lag's own term (-ydot(t - tau) * dtau/dp)
#' and the jump each dose makes when it reaches t_dose + tau, so these
#' parameters are finite-differenced, like a dosing-event parameter.
#'
#' @param s symengine environment
#' @return list with integer vectors `eta` and `theta` (parameter numbers)
#' @noRd
.foceiDelayLagParams <- function(s) {
  .ret <- list(eta = integer(0), theta = integer(0))
  .vars <- c(
    ls(s, pattern = "^rx__d_dt_.*__$", all.names = TRUE),
    intersect(c("rx_pred_", "rx_r_"), ls(s, all.names = TRUE))
  )
  .taus <- character(0)
  .walk <- function(x) {
    if (is.call(x)) {
      if (identical(x[[1]], quote(delay)) && length(x) == 3L) {
        .taus <<- c(.taus, deparse1(x[[3]]))
      }
      for (.i in seq_along(x)[-1L]) {
        .walk(x[[.i]])
      }
    }
  }
  for (.v in .vars) {
    .e <- get(.v, envir = s) # rxFromSE() captures its argument
    .txt <- rxode2::rxFromSE(.e)
    if (!grepl("delay(", .txt, fixed = TRUE)) {
      next
    }
    .walk(str2lang(.txt))
  }
  if (length(.taus) == 0L) {
    return(.ret)
  }
  .num <- function(what) {
    .m <- regmatches(.taus, gregexpr(paste0("\\b", what, "\\[[0-9]+\\]"), .taus))
    .m <- unlist(.m, use.names = FALSE)
    sort(unique(as.integer(sub(".*\\[([0-9]+)\\]$", "\\1", .m))))
  }
  .ret$eta <- .num("ETA")
  .ret$theta <- .num("THETA")
  .ret
}
