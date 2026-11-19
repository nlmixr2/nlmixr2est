#' The omega-collapse warning threshold of a fit
#'
#' @param fit the finalized fit
#' @return the fraction of the initial omega below which a diagonal omega is
#'   called collapsed; `NA` when the warning is disabled
#' @noRd
.omegaCollapseWarnLimit <- function(fit) {
  .ctl <- tryCatch(fit$control, error = function(e) NULL)
  .v <- if (is.list(.ctl)) .ctl[["omegaCollapseWarn"]] else NULL
  if (is.null(.v)) {
    return(1e-3)
  }
  .v <- suppressWarnings(as.double(.v))
  if (length(.v) != 1L || !is.finite(.v) || .v <= 0) {
    return(NA_real_)
  }
  .v
}

#' The diagonal omega lower-bound divisor a fit's optimizer enforced
#'
#' @param fit the finalized fit
#' @return `diagOmegaBoundLower` from a FOCEi-family control, or `NA` when the
#'   fit has no such bound
#' @noRd
.omegaCollapseBoundLower <- function(fit) {
  .ctl <- tryCatch(fit$control, error = function(e) NULL)
  .v <- if (is.list(.ctl)) .ctl[["diagOmegaBoundLower"]] else NULL
  .v <- suppressWarnings(as.double(.v))
  if (length(.v) != 1L || !is.finite(.v) || .v <= 1) {
    return(NA_real_)
  }
  .v
}

#' The initial values of the free diagonal omegas
#'
#' @param iniDf a model's `iniDf` before estimation
#' @return named numeric vector of the positive, unfixed diagonal omegas
#' @noRd
.omegaCollapseIni <- function(iniDf) {
  if (!is.data.frame(iniDf) || !all(c("neta1", "neta2", "fix", "est", "name") %in% names(iniDf))) {
    return(setNames(numeric(0), character(0)))
  }
  .d <- iniDf[
    !is.na(iniDf$neta1) & iniDf$neta1 == iniDf$neta2 & !iniDf$fix & is.finite(iniDf$est) & iniDf$est > 0,
    ,
    drop = FALSE
  ]
  setNames(.d$est, .d$name)
}

#' Diagonal omegas that collapsed or ended at their lower bound
#'
#' @param ini named initial diagonal omegas (free only)
#' @param final named final diagonal omegas
#' @param limit fraction of the initial omega that counts as collapsed; `NA`
#'   names nothing
#' @param boundLower the `diagOmegaBoundLower` divisor the optimizer enforced,
#'   or `NA` when there was none
#' @param omega0,omega the full initial and final omega matrices; when given,
#'   the bound is checked on the parameters FOCEi bounds (see
#'   `.omegaCollapseAtBound()`), otherwise on the variance ratio
#' @return list with `collapsed` and `bound`, the offending eta names
#' @noRd
.omegaCollapsed <- function(ini, final, limit, boundLower = NA_real_, omega0 = NULL, omega = NULL) {
  .none <- list(collapsed = character(0), bound = character(0))
  if (length(limit) != 1L || !is.finite(limit)) {
    return(.none)
  }
  .nm <- intersect(names(ini), names(final))
  if (length(.nm) == 0L) {
    return(.none)
  }
  .ratio <- final[.nm] / ini[.nm]
  .ok <- is.finite(.ratio)
  .collapsed <- .nm[.ok & .ratio < limit]
  .bound <- character(0)
  if (length(boundLower) == 1L && is.finite(boundLower)) {
    # within 1% of start/boundLower: the optimizer stopped on the bound.  The
    # variance is the bounded parameter only for an eta outside any block.
    .ratioAt <- .nm[.ok & .ratio <= 1.01 / boundLower]
    .cmp <- .omegaCollapseAtBound(omega0, omega, boundLower)
    .at <- if (!is.null(.cmp)) {
      c(.cmp$at, intersect(.ratioAt, intersect(.cmp$unknown, .omegaCollapseUnblocked(omega0))))
    } else if (is.matrix(omega0) && !is.null(rownames(omega0))) {
      intersect(.ratioAt, .omegaCollapseUnblocked(omega0))
    } else {
      .ratioAt
    }
    .bound <- setdiff(intersect(.at, .nm), .collapsed)
  }
  list(collapsed = .collapsed, bound = .bound)
}

#' Etas with no covariance in an omega
#'
#' @param omega a named omega matrix
#' @return the names of the etas whose off-diagonal omega entries are all 0
#' @noRd
.omegaCollapseUnblocked <- function(omega) {
  .off <- omega
  diag(.off) <- 0
  rownames(omega)[rowSums(.off != 0) == 0]
}

#' Etas whose bounded omega parameter ended at its lower bound
#'
#' FOCEi bounds the diagonal of `chol(solve(omega))` (through `diagXform`, a
#' monotone map), with the lower-variance bound taken from
#' `omega0/boundLower`.  In a block that is not the marginal variance, so the
#' ratio rule can miss a pinned eta; for a diagonal omega the two agree.  Each
#' independent block is checked on its own (exact, since conditioning on
#' independent etas changes nothing), so one that cannot be inverted only
#' leaves its own etas unknown.
#' @param omega0,omega full initial and final omega matrices with matching
#'   dimnames
#' @param boundLower the `diagOmegaBoundLower` divisor
#' @return list with `at`, the etas within 1% of the bound, and `unknown`, the
#'   etas of blocks that could not be checked; `NULL` when the matrices do not
#'   match
#' @noRd
.omegaCollapseAtBound <- function(omega0, omega, boundLower) {
  if (
    !is.matrix(omega0) ||
      !is.matrix(omega) ||
      is.null(rownames(omega)) ||
      !identical(dimnames(omega0), dimnames(omega))
  ) {
    return(NULL)
  }
  # on the correlation scale, so a tiny variance cannot make solve() reject
  # the block
  .d <- function(m) {
    .s <- sqrt(diag(m))
    if (!all(is.finite(.s) & .s > 0)) {
      return(NULL)
    }
    tryCatch(diag(chol(solve(m / outer(.s, .s)))) / .s, error = function(e) NULL)
  }
  .at <- character(0)
  .unknown <- character(0)
  for (.b in .omegaBlocks(omega0 != 0 | omega != 0)) {
    .d0 <- .d(omega0[.b, .b, drop = FALSE])
    .d1 <- .d(omega[.b, .b, drop = FALSE])
    if (is.null(.d0) || is.null(.d1)) {
      .unknown <- c(.unknown, rownames(omega)[.b])
    } else {
      .at <- c(.at, rownames(omega)[.b][is.finite(.d1) & 1.01 * .d1^2 >= boundLower * .d0^2])
    }
  }
  list(at = .at, unknown = .unknown)
}

#' The independent blocks of an omega pattern
#'
#' @param nz logical matrix, `TRUE` where two etas covary
#' @return list of index vectors, one per connected block, in eta order
#' @noRd
.omegaBlocks <- function(nz) {
  .n <- nrow(nz)
  .grp <- seq_len(.n)
  repeat {
    .new <- vapply(seq_len(.n), function(i) min(.grp[nz[i, ] | seq_len(.n) == i]), integer(1))
    if (identical(.new, .grp)) {
      break
    }
    .grp <- .new
  }
  unname(split(seq_len(.n), .grp))
}

#' Warn when a fit's diagonal omega collapsed toward zero
#'
#' An outer optimizer can drive a diagonal omega to (or near) zero, where the
#' objective is flat in it, and report convergence while a much larger omega
#' gives a lower objective (#1394).  Emits concise warnings (collected into
#' `$runInfo`) naming the etas whose omega ended far below its initial
#' estimate or at its `diagOmegaBoundLower` bound.
#' @param fit the finalized fit
#' @param env the estimation environment; its `ui` holds the initial estimates
#' @return the list from `.omegaCollapsed()`, invisibly
#' @noRd
.omegaCollapseWarn <- function(fit, env) {
  .limit <- .omegaCollapseWarnLimit(fit)
  .ret <- list(collapsed = character(0), bound = character(0))
  if (is.na(.limit)) {
    return(invisible(.ret))
  }
  .ui <- if (is.environment(env)) get0("ui", envir = env, inherits = FALSE)
  .ini <- .omegaCollapseIni(tryCatch(.ui$iniDf, error = function(e) NULL))
  .om <- tryCatch(fit$omega, error = function(e) NULL)
  if (length(.ini) == 0L || !is.matrix(.om) || is.null(rownames(.om))) {
    return(invisible(.ret))
  }
  .ret <- .omegaCollapsed(
    .ini,
    setNames(diag(.om), rownames(.om)),
    .limit,
    .omegaCollapseBoundLower(fit),
    tryCatch(.ui$omega, error = function(e) NULL),
    .om
  )
  if (length(.ret$collapsed) > 0L) {
    warning(.omegaCollapseWarnMsg(.ret$collapsed, .limit), call. = FALSE)
  }
  if (length(.ret$bound) > 0L) {
    warning(.omegaCollapseWarnMsg(.ret$bound, .omegaCollapseBoundLower(fit), bound = TRUE), call. = FALSE)
  }
  invisible(.ret)
}

#' The omega-collapse warning, kept to one `$runInfo` line
#'
#' @param etas the offending eta names
#' @param value the collapse fraction, or the bound divisor when `bound`
#' @param bound `TRUE` for omegas stopped on their lower bound
#' @return the message
#' @noRd
.omegaCollapseWarnMsg <- function(etas, value, bound = FALSE) {
  .pre <- if (bound) {
    sprintf("omega at lower bound (initial/%s): ", format(value))
  } else {
    sprintf("omega collapsed (<%s x initial): ", format(value))
  }
  .runInfoNameMsg(.pre, etas)
}
