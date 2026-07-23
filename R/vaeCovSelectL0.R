# L0Learn-backed candidate generation for the VAE covariate M-step.  This is the
# ONLY file that names L0Learn: it produces CANDIDATE supports, which the C++
# kernel then re-scores with the exact objective RSS_S/omega + penalty*|S| and
# polishes.  L0Learn therefore cannot shift a selection through its own internal
# scaling or lambda grid -- it only decides which subsets get looked at.

#' Is the suggested L0Learn package usable?
#' @return `TRUE` when `L0Learn` can be loaded
#' @noRd
.vaeL0Available <- function() {
  requireNamespace("L0Learn", quietly = TRUE)
}

#' Distinct supports visited on one L0Learn regularization path.
#' @param x covariate matrix (no intercept column)
#' @param y response vector
#' @param penalty `"L0"` or `"L0L2"`
#' @return list of ascending 0-based integer vectors (empty list on failure)
#' @noRd
.vaeL0Path <- function(x, y, penalty) {
  .f <- tryCatch(suppressWarnings(
    L0Learn::L0Learn.fit(x, y, penalty = penalty, algorithm = "CDPSI",
                         maxSuppSize = ncol(x))),
    error = function(e) NULL)
  if (is.null(.f) || !length(.f$beta)) return(list())
  # beta is one p x nLambda sparse matrix per gamma; each column is one support
  unlist(lapply(.f$beta, function(b) {
    b <- as.matrix(b)
    lapply(seq_len(ncol(b)), function(j) which(b[, j] != 0) - 1L)
  }), recursive = FALSE)
}

#' Candidate supports for one latent dimension.
#'
#' Pools the `L0` and `L0L2` paths -- `L0L2` is the reliable one on correlated
#' designs, which is the pharmacometric norm.  Extra candidates can only improve
#' the answer because the caller re-scores every one of them.
#' @param x covariate matrix (no intercept column)
#' @param y response vector
#' @return list of ascending 0-based integer vectors, always including the empty
#'   support so the intercept-only model is never lost
#' @noRd
.vaeL0Supports <- function(x, y) {
  .empty <- list(integer(0))
  if (!is.matrix(x) || ncol(x) == 0L || nrow(x) < 3L) return(.empty)
  if (!all(is.finite(x)) || !all(is.finite(y))) return(.empty)
  .all <- c(.empty, .vaeL0Path(x, y, "L0"), .vaeL0Path(x, y, "L0L2"))
  .all[!duplicated(vapply(.all, paste, character(1), collapse = ","))]
}

#' Resolve the per-latent-dimension covariate-selection mode.
#'
#' Pure: it emits nothing.  The caller raises the returned messages so they land
#' in the fit's `$runInfo`, which is how run-time notes reach the user.
#' @param nCand integer vector of candidate covariate counts, one per latent
#'   dimension (after `pinCovariates` trimming -- the real search size)
#' @param control a `vaeControl()` list
#' @param avail is L0Learn usable (an argument so the fallback branches are
#'   testable without installing/uninstalling the package)
#' @return list with `mode` (0 = exact branch-and-bound, 1 = L0Learn per
#'   dimension), `used` (`"bnb"`, `"l0learn"` or `"mixed"`) and `msg`
#' @noRd
.vaeCovSelectModes <- function(nCand, control, avail = .vaeL0Available()) {
  nCand <- as.integer(nCand)
  .method <- control$covSelectMethod
  if (is.null(.method)) .method <- "auto"
  .avail <- avail
  .msg <- character(0)
  .mode <- integer(length(nCand))
  if (.method == "l0learn") {
    if (!.avail) {
      stop("covSelectMethod=\"l0learn\" needs the L0Learn package; install it",
           call. = FALSE)
    }
    .mode[nCand > 0L] <- 1L
  } else if (.method == "auto") {
    .big <- nCand >= as.integer(control$covSelectMaxExact)
    if (any(.big)) {
      if (.avail) {
        .mode[.big] <- 1L
      } else {
        .msg <- c(.msg, "exact covariate search is slow here; install L0Learn")
      }
    }
  }
  .used <- if (all(.mode == 1L)) "l0learn" else if (any(.mode == 1L)) "mixed" else "bnb"
  if (.used != "bnb") {
    .msg <- c(.msg, "covariate search used L0Learn, not the exact search")
  }
  list(mode = .mode, used = .used, msg = .msg)
}

#' Candidate generator handed to the C++ M-step.
#'
#' Called ONCE per training iteration on the main thread (an `Rcpp::Function`
#' cannot be called from inside the OpenMP loop that scores the candidates).
#' @param y response matrix, one column per latent dimension
#' @param covMat covariate design (no intercept column)
#' @param mode per-dimension mode from [.vaeCovSelectModes()]
#' @param allowed list of 0-based allowed covariate columns per dimension, or
#'   `NULL` when every covariate is a candidate
#' @return list of length `ncol(y)`; each element is a list of candidate supports
#'   as 0-based GLOBAL covariate indices, or `NULL` for a branch-and-bound
#'   dimension
#' @noRd
.vaeL0Candidates <- function(y, covMat, mode, allowed = NULL) {
  lapply(seq_len(ncol(y)), function(k) {
    if (mode[k] != 1L) return(NULL)
    .cols <- if (is.null(allowed)) seq_len(ncol(covMat)) - 1L else allowed[[k]]
    if (!length(.cols)) return(list(integer(0)))
    .sup <- .vaeL0Supports(covMat[, .cols + 1L, drop = FALSE], y[, k])
    # map the reduced-design indices back to global covariate columns
    lapply(.sup, function(s) as.integer(.cols[s + 1L]))
  })
}
