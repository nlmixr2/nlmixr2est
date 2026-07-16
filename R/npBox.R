# Support-point box for the nonparametric engines.  A support point is an eta
# vector, so the grid lives in eta space; this defines the per-eta box the Sobol
# initial grid covers.  Control-selectable via gridBounds:
#   "auto" (default) -- symmetric box +/- gridWidth * initial eta SD (from the
#      ini omega diagonal); falls back to +/- 1.2 for a degenerate/zero variance.
#   "ini"  -- use the ini lower/upper of the population parameter mu-referenced by
#      each eta when finite (interpreted on the eta / transformed scale), else auto.
#   "both" -- ini bounds when present, auto otherwise.

#' Eta-space box for the nonparametric support-point grid.
#' @param ui rxode2 ui (bounded-transformed)
#' @param control impmapControl-derived control (reads gridBounds, gridWidth)
#' @return list(lower, upper, names) ordered by eta index (neta1)
#' @noRd
.npEtaBox <- function(ui, control) {
  .iniDf <- ui$iniDf
  .eta <- .iniDf[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2, , drop = FALSE]
  .eta <- .eta[order(.eta$neta1), , drop = FALSE]
  if (nrow(.eta) == 0L) stop("npag/npb require at least one random effect (eta)", call. = FALSE)
  .w <- if (is.null(control$gridWidth)) 4 else as.numeric(control$gridWidth)
  .gb <- if (is.null(control$gridBounds)) "auto" else match.arg(control$gridBounds, c("auto", "ini", "both"))
  .sd <- sqrt(pmax(.eta$est, 0))
  .half <- .w * .sd
  .half[!is.finite(.half) | .half <= 0] <- 1.2   # degenerate/zero-variance fallback
  .lower <- -.half
  .upper <- .half
  if (.gb %in% c("ini", "both")) {
    # override with the mu-referenced parameter's ini bounds where finite
    .mr <- ui$muRefDataFrame
    for (i in seq_len(nrow(.eta))) {
      .th <- .mr$theta[match(.eta$name[i], .mr$eta)]
      if (is.na(.th)) next
      .row <- .iniDf[!is.na(.iniDf$ntheta) & .iniDf$name == .th, , drop = FALSE]
      if (nrow(.row) != 1L) next
      .lo <- .row$lower; .hi <- .row$upper; .est <- .row$est
      if (is.finite(.lo) && is.finite(.hi)) {
        # bounds are on the population (transformed) scale; center on the estimate
        .lower[i] <- .lo - .est
        .upper[i] <- .hi - .est
      } else if (identical(.gb, "ini")) {
        # ini requested but unbounded -> keep the auto half-width (already set)
      }
    }
  }
  list(lower = .lower, upper = .upper, names = .eta$name,
       fixed = !is.na(.eta$fix) & .eta$fix)
}
