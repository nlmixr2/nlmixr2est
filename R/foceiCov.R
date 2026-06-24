# Full FOCEI/FOCE covariance (structural theta + residual sigma + Omega, diagonal or block)
# with a three-tier robustness ladder (Matthew's note that the augmented ODE is
# more likely to fail solving the higher its order):
#
#   1. analytic, EXACT 3rd-order sensitivities      (foceiCovAnalytic sens="exact3")
#   2. analytic, 2nd-order + Shi(2021) finite differences of the 2nd-order
#      sensitivities for the 3rd-order term         (foceiCovAnalytic sens="fd2")
#        -- a much lighter ODE (no O(neta^3) states) that solves where tier 1
#           does not, and reproduces it to ~1e-5 (it differences EXACT
#           derivatives, not the objective, so no catastrophic cancellation).
#   3. builtin: finite differences of nlmixr2's OWN focei objective over the
#      non-Cholesky parameters                       (.foceiCovBuiltinFD)
#        -- uses no sensitivity model at all (just the fitted model), so it works
#           when symengine / the augmented build is unavailable; tightened OFV
#           (sigdig) keeps the FD Hessian stable.  Returns Omega/residual SEs too.
#
# Every tier returns the SAME full covariance (theta + sigma + Omega, diagonal or
# block), so the fallback never silently drops the Omega/residual blocks.

#' Omega blocks (connected components of the random-effect covariance), so each
#' block can be set with one `ini()` formula.  Returns a list of eta-index vectors.
#' @noRd
.omegaBlocks <- function(Om, tol = 1e-10) {
  n <- nrow(Om); adj <- abs(Om) > tol; diag(adj) <- TRUE
  comp <- integer(n); k <- 0L
  for (i in seq_len(n)) if (comp[i] == 0L) {
    k <- k + 1L; q <- i
    while (length(q)) { v <- q[1]; q <- q[-1]
      if (comp[v] == 0L) { comp[v] <- k; q <- c(q, which(adj[v, ] & comp == 0L)) } }
  }
  unname(split(seq_len(n), comp))
}

#' Which rows of `pairs` (each `c(a, b)`) are FIXED Omega elements in `idf`.
#' @noRd
.omegaFixed <- function(idf, pairs) {
  if (nrow(pairs) == 0L) return(logical(0))
  vapply(seq_len(nrow(pairs)), function(k) {
    r <- which(!is.na(idf$neta1) & ((idf$neta1 == pairs[k, 1] & idf$neta2 == pairs[k, 2]) |
                                    (idf$neta1 == pairs[k, 2] & idf$neta2 == pairs[k, 1])))
    length(r) > 0L && isTRUE(idf$fix[r[1]])
  }, logical(1))
}

#' Set one Omega block on a UI by its `eta_i + eta_j + ... ~ c(...)` ini formula.
#' @noRd
.omBlockIni <- function(u, etas, vals) {
  lhs <- Reduce(function(x, y) call("+", x, as.name(y)), etas[-1], as.name(etas[1]))
  rhs <- if (length(vals) == 1L) vals[[1]] else as.call(c(quote(c), as.list(vals)))
  eval(bquote(rxode2::ini(u, .(call("~", lhs, rhs)))))
}

#' Central finite-difference Hessian of `fn` at `x` (NULL if any evaluation fails).
#' @noRd
.fdHessian <- function(fn, x, h) {
  np <- length(x); H <- matrix(0, np, np)
  for (i in seq_len(np)) for (j in i:np) {
    hi <- h * (abs(x[i]) + 1); hj <- h * (abs(x[j]) + 1)      # absolute floor h: no step->0 for ~0 params
    pp <- x; pp[i] <- pp[i]+hi; pp[j] <- pp[j]+hj
    pm <- x; pm[i] <- pm[i]+hi; pm[j] <- pm[j]-hj
    mp <- x; mp[i] <- mp[i]-hi; mp[j] <- mp[j]+hj
    mm <- x; mm[i] <- mm[i]-hi; mm[j] <- mm[j]-hj
    H[i, j] <- H[j, i] <- (fn(pp) - fn(pm) - fn(mp) + fn(mm)) / (4 * hi * hj)
  }
  if (anyNA(H)) NULL else H
}

#' Builtin finite-difference covariance from nlmixr2's native objective.
#'
#' Last-resort tier: re-evaluates nlmixr2's own focei objective at perturbed
#' parameters via a 0-outer-iteration fit (no sensitivity model needed) and
#' finite-differences it over the thetas and the FULL Omega lower-triangle
#' (variances AND covariances, non-Cholesky) -- so it covers block Omega.  Every
#' parameter is set through the same `ini()` interface used for the point fit:
#' thetas by named value, each Omega block by its `eta_i + eta_j ~ c(...)`
#' formula.  `sigdig` tightens the inner-EBE / ODE tolerance so the OFV is smooth
#' enough to difference.
#' @noRd
.foceiCovBuiltinFD <- function(fit, h = 1e-3, sigdig = 9) {
  ui <- fit$finalUi; idf <- ui$iniDf; muRef <- ui$muRefDataFrame
  Om <- fit$omega; neta <- nrow(Om)
  data <- tryCatch(getData(fit), error = function(e) NULL); if (is.null(data)) return(NULL)
  fixed <- function(v) !is.na(v) & v
  thRows <- which(!is.na(idf$ntheta) & !fixed(idf$fix)); thNames <- idf$name[thRows]  # estimated thetas only
  diagRows <- idf[!is.na(idf$neta1) & idf$neta1 == idf$neta2, , drop = FALSE]
  etaName <- diagRows$name[order(diagRows$neta1)]               # etaName[i] = name of eta i
  blocks <- .omegaBlocks(Om)
  pairs <- do.call(rbind, lapply(blocks, function(b)            # within-block lower-tri (i>=j)
    do.call(rbind, lapply(seq_along(b), function(a) cbind(b[a], b[seq_len(a)])))))
  pairs <- pairs[!.omegaFixed(idf, pairs), , drop = FALSE]      # estimated Omega elements only
  par0 <- unname(c(idf$est[thRows], Om[pairs])); nth <- length(thRows)
  inter <- tryCatch(as.integer(rxode2::rxGetControl(ui, "interaction", 1L)), error = function(e) 1L)
  ctl <- foceiControl(maxOuterIterations = 0L, maxInnerIterations = 1000L, covMethod = "",
                      print = 0L, sigdig = sigdig, interaction = inter)  # match FOCE/FOCEI of the fit
  objAt <- function(pv) {
    u <- tryCatch({
      thv <- setNames(pv[seq_len(nth)], thNames)             # via a variable: ini() uses substitute()
      u <- rxode2::ini(ui, thv)
      M <- Om; ov <- pv[-seq_len(nth)]                         # fixed Omega elements stay at the fit
      M[pairs] <- ov; M[pairs[, 2:1, drop = FALSE]] <- ov
      for (b in blocks) u <- .omBlockIni(u, etaName[b], unlist(lapply(seq_along(b), function(a) M[b[a], b[seq_len(a)]])))
      u
    }, error = function(e) NULL)
    if (is.null(u)) return(NA_real_)
    f <- try(suppressMessages(suppressWarnings(nlmixr2(u, data, "focei", ctl))), silent = TRUE)
    if (inherits(f, "try-error") || is.null(f$objf)) NA_real_ else f$objf
  }
  H <- .fdHessian(objAt, par0, h); if (is.null(H)) return(NULL)
  cov <- tryCatch(solve(0.5 * H), error = function(e) NULL); if (is.null(cov)) return(NULL)
  thOf <- function(e) { t <- muRef$theta[match(e, muRef$eta)]; if (is.na(t)) e else t }
  omNm <- apply(pairs, 1, function(p) if (p[1] == p[2]) paste0("om.", thOf(etaName[p[1]]))
                else paste0("cov.", thOf(etaName[p[1]]), ".", thOf(etaName[p[2]])))
  nm <- c(thNames, omNm); dimnames(cov) <- list(nm, nm)
  list(cov = cov, se = setNames(suppressWarnings(sqrt(diag(cov))), nm), params = nm, method = "fd-builtin")
}

#' Full covariance for a fitted nlmixr2 FOCEI/FOCE object
#'
#' Returns the full observed-information covariance over the structural thetas,
#' residual sigma, and Omega (diagonal or block) variances and covariances, via a
#' three-tier ladder: analytic exact 3rd-order, then 2nd-order with Shi finite differences,
#' then a finite difference of nlmixr2's native objective.  Each tier returns the
#' same full set of parameters; the analytic path is never partially applied.
#'
#' @param fit a fitted nlmixr2 focei object.
#' @return list with `cov`, `se`, `params`, and `method` (`"analytic"` /
#'   `"analytic-fd2"` / `"fd-builtin"`), or `NULL` if no tier succeeds.
#' @export
foceiCov <- function(fit) {
  # each tier guarded: an unexpected error (e.g. singular Omega solve) drops to the next
  r <- tryCatch(foceiCovAnalytic(fit, sens = "exact3"), error = function(e) NULL)
  if (!is.null(r)) { r$method <- "analytic"; return(r) }
  r <- tryCatch(foceiCovAnalytic(fit, sens = "fd2"), error = function(e) NULL)
  if (!is.null(r)) { r$method <- "analytic-fd2"; return(r) }
  tryCatch(.foceiCovBuiltinFD(fit), error = function(e) NULL)
}
