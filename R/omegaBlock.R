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
    if (is.na(.i) || is.na(.j)) {
      next
    }
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
  ## Block on the model's DECLARED structure, not on which fitted values happen
  ## to be non-zero.  An estimated covariance of exactly 0 (the correlation hold
  ## running to the end of a short fit) would otherwise disconnect the etas, emit
  ## them as separate singletons, and leave the existing covariance row sitting
  ## at its ini value -- so the reported omega would disagree with the omega the
  ## fit actually used.
  .decl <- tryCatch(.omegaBlockFromIniDf(rxode2::rxUiDecompress(u)$iniDf, .nm)$mat, error = function(e) NULL)
  .adj <- if (is.null(.decl)) omegaMat != 0 else (.decl != 0 | omegaMat != 0)
  .comp <- .omegaBlockIds(.adj)
  for (.b in unique(.comp)) {
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
      .expr <- paste0(paste(.nm[.idx], collapse = " + "), " ~ c(", paste(.vals, collapse = ", "), ")")
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
  if (is.null(.om)) {
    .om <- diag(as.numeric(fit$omega), length(etaNames))
  }
  dimnames(.om) <- list(etaNames, etaNames)
  .om
}

#' Connected-component block ids of a symmetric adjacency pattern.
#'
#' @param adj logical (or coercible) symmetric matrix; the diagonal is ignored
#' @return integer vector of block ids, one per row
#' @noRd
.omegaBlockIds <- function(adj) {
  .adj <- adj != 0
  .n <- nrow(.adj)
  diag(.adj) <- TRUE
  .comp <- integer(.n)
  .c <- 0L
  for (.i in seq_len(.n)) {
    if (.comp[.i] != 0L) {
      next
    }
    .c <- .c + 1L
    .stack <- .i
    while (length(.stack)) {
      .v <- .stack[[1L]]
      .stack <- .stack[-1L]
      if (.comp[.v] != 0L) {
        next
      }
      .comp[.v] <- .c
      .stack <- c(.stack, which(.adj[.v, ] & .comp == 0L))
    }
  }
  .comp
}

#' Zeros `rxSymInvCholCreate()` cannot hold at zero.
#'
#' It counts its parameters from omega's zero pattern but fills them from each
#' block's cholesky factor, and it takes a block to be the whole index SPAN of
#' a correlated group.  So the patterns it accepts are exactly those whose
#' connected components are contiguous index ranges, each one dense
#' (rxode2#1365); anything else makes the two counts disagree and the theta
#' setter refuses the matrix with "theta has to have N elements".
#'
#' Measured over every 4x4 pattern: "components are contiguous and dense"
#' matches which matrices the call accepts 64/64, where "dense components"
#' alone misses 7 of them.
#'
#' Closing each component up to its span can merge components (spans overlap),
#' so grow the pattern to a fixed point.
#'
#' @param mat symmetric matrix
#' @return two-column (row, col) matrix of upper-triangle positions that have
#'   to become nonzero, empty when the pattern is already acceptable
#' @noRd
.omegaBlockZeros <- function(mat) {
  .adj <- mat != 0
  .adj[is.na(.adj)] <- FALSE
  repeat {
    .comp <- .omegaBlockIds(.adj)
    .new <- .adj
    for (.c in unique(.comp)) {
      .idx <- which(.comp == .c)
      .span <- seq.int(min(.idx), max(.idx))
      .new[.span, .span] <- TRUE
    }
    if (identical(.new, .adj)) {
      break
    }
    .adj <- .new
  }
  which(upper.tri(mat) & .adj & mat == 0, arr.ind = TRUE)
}

#' Fill the block-internal zeros of `mat` with a negligible covariance.
#'
#' @param mat symmetric positive-definite matrix
#' @param cor correlation written into each filled cell
#' @return the filled matrix, or `NULL` when there is nothing to fill or the
#'   fill cannot be made (a non-positive diagonal, or a result that is no
#'   longer cholesky-able)
#' @noRd
.omegaFillBlockZeros <- function(mat, cor = 1e-10) {
  .idx <- .omegaBlockZeros(mat)
  if (nrow(.idx) == 0L) {
    return(NULL)
  }
  .d <- diag(mat)
  .ret <- mat
  for (.k in seq_len(nrow(.idx))) {
    .i <- .idx[.k, 1L]
    .j <- .idx[.k, 2L]
    .v <- cor * sqrt(.d[.i] * .d[.j])
    if (!is.finite(.v) || .v <= 0) {
      return(NULL)
    }
    .ret[.i, .j] <- .ret[.j, .i] <- .v
  }
  if (inherits(try(chol(.ret), silent = TRUE), "try-error")) {
    return(NULL)
  }
  .ret
}

#' The random effects named by `.omegaBlockZeros()` positions, comma separated
#' and truncated so the warning stays on one line.
#' @noRd
.omegaBlockZeroNames <- function(mat, idx, width = 35L) {
  .nm <- colnames(mat)
  if (is.null(.nm)) {
    .nm <- paste0("eta", seq_len(nrow(mat)))
  }
  .use <- unique(as.vector(idx))
  .txt <- paste(.nm[.use], collapse = ", ")
  if (nchar(.txt) > width) {
    .txt <- paste0(substr(.txt, 1L, width - 3L), "...")
  }
  .txt
}
