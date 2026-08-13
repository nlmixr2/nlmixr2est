# L0Learn-backed candidate generation for the VAE covariate M-step.  This is the
# ONLY file that names L0Learn (an Imports dependency, so it is always available):
# it produces CANDIDATE supports, which the C++
# kernel then re-scores with the exact objective RSS_S/omega + penalty*|S| and
# polishes.  L0Learn therefore cannot shift a selection through its own internal
# scaling or lambda grid -- it only decides which subsets get looked at.

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
#'   support so the intercept-only model is never lost.  Carries a
#'   `vaeL0Fail` attribute when L0Learn was asked for proposals on a usable
#'   design and returned none -- the caller counts those, because such a
#'   dimension silently degrades to a local search from the intercept-only
#'   model, which is a far weaker search than the one that was requested.
#' @noRd
.vaeL0Supports <- function(x, y) {
  .empty <- list(integer(0))
  if (!is.matrix(x) || ncol(x) == 0L || nrow(x) < 3L) return(.empty)
  if (!all(is.finite(x)) || !all(is.finite(y))) return(.empty)
  .p <- c(.vaeL0Path(x, y, "L0"), .vaeL0Path(x, y, "L0L2"))
  ## Both penalty paths empty on a design that passed the guards above means
  ## L0Learn itself failed (or proposed nothing), not that the input was
  ## unusable -- the two cases must not be conflated, since only this one is
  ## worth telling the user about.
  if (!length(.p)) return(structure(.empty, vaeL0Fail = TRUE))
  .all <- c(.empty, .p)
  .all[!duplicated(vapply(.all, paste, character(1), collapse = ","))]
}

#' Run-time note for latent dimensions whose L0Learn proposals failed.
#'
#' Pure: the caller raises it so it lands in the fit's `$runInfo`.
#' @param n number of dimension-iterations that fell back
#' @return character(0) or a single message
#' @noRd
.vaeL0FailMsg <- function(n) {
  if (!length(n) || is.na(n) || n <= 0L) return(character(0))
  "L0Learn proposed no supports; a weaker search was used"
}

#' Resolve the per-latent-dimension covariate-selection mode.
#'
#' Pure: it emits nothing.  The caller raises the returned messages so they land
#' in the fit's `$runInfo`, which is how run-time notes reach the user.
#' @param nCand search size per latent dimension, in BITS of feasible-support
#'   space: `sum over allowed groups of log2(1 + columns in that group)` after
#'   `pinCovariates` trimming.  One column per group costs exactly 1 bit, so a
#'   single-shape search reads as a plain candidate count (its historic meaning),
#'   while two shape families of one covariate cost `log2(3)` -- the exact search
#'   then gets the same worst-case node budget either way.
#' @param control a `vaeControl()` list
#' @return list with `mode` (0 = exact branch-and-bound, 1 = L0Learn per
#'   dimension), `used` (`"bnb"`, `"l0learn"` or `"mixed"`) and `msg`
#' @noRd
.vaeCovSelectModes <- function(nCand, control) {
  ## kept numeric: a bit budget is fractional once a covariate carries more than
  ## one shape, and truncating it would quietly widen the exact-search region
  nCand <- as.numeric(nCand)
  .method <- control$covSelectMethod
  if (is.null(.method)) .method <- "auto"
  .msg <- character(0)
  .mode <- integer(length(nCand))
  if (.method == "l0learn") {
    .mode[nCand > 0] <- 1L
  } else if (.method == "auto") {
    # covSelectMaxExact may be Inf (force the exact search everywhere), so compare
    # against the raw value rather than coercing it to integer (as.integer(Inf) is
    # NA).
    .mode[nCand >= control$covSelectMaxExact] <- 1L
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
#' @param failEnv optional environment with an integer `n`, incremented once per
#'   latent dimension whose L0Learn proposals failed.  An environment because
#'   this is called once per training iteration from C++, so the count has to
#'   survive across calls.
#' @return list of length `ncol(y)`; each element is a list of candidate supports
#'   indexed into that dimension's REDUCED design (the caller already maps those
#'   back to global covariate columns), or `NULL` for a branch-and-bound dimension
#' @noRd
.vaeL0Candidates <- function(y, covMat, mode, allowed = NULL, failEnv = NULL) {
  lapply(seq_len(ncol(y)), function(k) {
    if (mode[k] != 1L) return(NULL)
    .cols <- if (is.null(allowed)) seq_len(ncol(covMat)) - 1L else allowed[[k]]
    if (!length(.cols)) return(list(integer(0)))
    .s <- .vaeL0Supports(covMat[, .cols + 1L, drop = FALSE], y[, k])
    if (isTRUE(attr(.s, "vaeL0Fail"))) {
      if (!is.null(failEnv)) failEnv$n <- failEnv$n + 1L
      ## strip before it crosses into C++, which reads this as a plain list
      attr(.s, "vaeL0Fail") <- NULL
    }
    .s
  })
}
