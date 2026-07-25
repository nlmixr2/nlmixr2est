# omegaBlock.R -- shared helper: build the full ini omega block matrix (and the
# per-entry FIXED status) from an iniDf, in a given eta order.  Used by the
# est="vae" and est="emvi"/"fbvi" data preps, which estimate the full modeled omega
# (diagonal + declared off-diagonals); the estimation mask is the matrix's
# nonzero structure, matching rxSymInvCholCreate (a correlation initialized at
# exactly 0 is structurally absent, as in focei).

#' @return list(mat = neta x neta ini omega, fixMat = logical neta x neta)
#' @noRd
.omegaBlockFromIniDf <- function(idf, etaNames) {
  .neta <- length(etaNames)
  .mat <- matrix(0, .neta, .neta, dimnames = list(etaNames, etaNames))
  .fix <- matrix(FALSE, .neta, .neta, dimnames = list(etaNames, etaNames))
  .etaRows <- idf[!is.na(idf$neta1), , drop = FALSE]
  ## eta index (position in etaNames) keyed by the iniDf neta numbering
  .diagRows <- .etaRows[.etaRows$neta1 == .etaRows$neta2, , drop = FALSE]
  .idx <- setNames(match(.diagRows$name, etaNames), as.character(.diagRows$neta1))
  for (.r in seq_len(nrow(.etaRows))) {
    ## single-bracket lookup: an unmatched neta number gives NA, not an error
    .i <- .idx[as.character(.etaRows$neta1[.r])]
    .j <- .idx[as.character(.etaRows$neta2[.r])]
    if (is.na(.i) || is.na(.j)) next
    .v <- as.numeric(.etaRows$est[.r])
    .f <- isTRUE(as.logical(.etaRows$fix[.r]))
    .mat[.i, .j] <- .mat[.j, .i] <- .v
    .fix[.i, .j] <- .fix[.j, .i] <- .f
  }
  list(mat = .mat, fixMat = .fix)
}

#' Does `mat` carry any modeled (nonzero) off-diagonal?
#' @noRd
.omegaHasOffDiag <- function(mat) {
  any(mat[upper.tri(mat)] != 0)
}

#' Write a fitted omega matrix into a ui's ini(), block-aware.
#'
#' Uncorrelated etas are written `eta ~ v`; each correlated block is written
#' with the block syntax `e1 + e2 ~ c(v11, v21, v22)` (lower-tri row-major).
#' Blocks are the connected components of the model's DECLARED structure unioned
#' with the fitted matrix's nonzeros, so a covariance estimated at exactly 0 is
#' still written into its block rather than silently left at its ini value.
#' @noRd
.omegaWriteIni <- function(u, omegaMat) {
  .nm <- colnames(omegaMat)
  .n <- nrow(omegaMat)
  ## Block on the model's DECLARED structure, not on which fitted values happen
  ## to be non-zero.  An estimated covariance of exactly 0 (the correlation hold
  ## running to the end of a short fit) would otherwise disconnect the etas, emit
  ## them as separate singletons, and leave the existing covariance row sitting
  ## at its ini value -- so the reported omega would disagree with the omega the
  ## fit actually used.
  .decl <- tryCatch(.omegaBlockFromIniDf(rxode2::rxUiDecompress(u)$iniDf, .nm)$mat,
                    error = function(e) NULL)
  .adj <- if (is.null(.decl)) omegaMat != 0 else (.decl != 0 | omegaMat != 0)
  diag(.adj) <- TRUE
  .comp <- integer(.n)
  .c <- 0L
  for (.i in seq_len(.n)) {
    if (.comp[.i] != 0L) next
    .c <- .c + 1L
    .stack <- .i
    while (length(.stack)) {
      .v <- .stack[[1L]]
      .stack <- .stack[-1L]
      if (.comp[.v] != 0L) next
      .comp[.v] <- .c
      .stack <- c(.stack, which(.adj[.v, ] & .comp == 0L))
    }
  }
  for (.b in seq_len(.c)) {
    .idx <- which(.comp == .b)
    if (length(.idx) == 1L) {
      .expr <- paste0(.nm[.idx], " ~ ", signif(omegaMat[.idx, .idx], 12))
    } else {
      .vals <- character(0)
      for (.r in seq_along(.idx)) {
        for (.s in seq_len(.r)) {
          .vals <- c(.vals, as.character(signif(omegaMat[.idx[.r], .idx[.s]], 12)))
        }
      }
      .expr <- paste0(paste(.nm[.idx], collapse = " + "), " ~ c(",
                      paste(.vals, collapse = ", "), ")")
    }
    u <- do.call(rxode2::ini, list(u, str2lang(.expr)))
  }
  u
}

#' The fitted omega as a dimnamed matrix from a vae/vi fit list: the full
#' `omegaMat` when present, else the diagonal vector.
#' @noRd
.omegaFitMat <- function(fit, etaNames) {
  .om <- fit$omegaMat
  if (is.null(.om)) .om <- diag(as.numeric(fit$omega), length(etaNames))
  dimnames(.om) <- list(etaNames, etaNames)
  .om
}
