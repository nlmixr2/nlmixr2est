#' Snapshot a FOCEi problem's inputs for a later omega-collapse probe
#'
#' The fit fills the environment it is given, so a probe re-fit needs the
#' inputs as they were before the first fit.
#' @param env the FOCEi environment, before fitting
#' @param control the fit's control
#' @return a list of the environment's bindings, or `NULL` when no probe applies
#' @noRd
.foceiCollapseProbeSnapshot <- function(env, control) {
  if (
    !isTRUE(control$collapseProbe) ||
      isTRUE(control$maxOuterIterations == 0) ||
      !inherits(env$rxInv, "rxSymInvCholEnv") ||
      length(env$mixIdx) > 0L
  ) {
    return(NULL)
  }
  as.list(env, all.names = TRUE)
}

#' The starting omega of a collapse probe
#'
#' @param omega the fit's final omega
#' @param ini named initial diagonal omegas
#' @param etas the collapsed etas
#' @param omega0 the starting omega, whose block pattern the probe keeps; a
#'   covariance the fit left at exactly 0 gets a negligible correlation
#' @return `omega` with each collapsed eta's variance raised back to its initial
#'   estimate and its correlations halved.  Shrinking (not zeroing) them keeps
#'   the block structure, so the omega parameters stay the same, and keeps the
#'   matrix positive definite.
#' @noRd
.foceiCollapseProbeOmega <- function(omega, ini, etas, omega0 = NULL) {
  .sd <- sqrt(diag(omega))
  .cor <- omega / outer(.sd, .sd)
  .cor[!is.finite(.cor)] <- 0
  if (is.matrix(omega0) && identical(dimnames(omega0), dimnames(omega))) {
    .cor[omega0 != 0 & .cor == 0] <- 1e-6
  }
  diag(.cor) <- 1
  .s <- setNames(rep(1, nrow(omega)), rownames(omega))
  .s[etas] <- 0.5
  .cor <- .cor * outer(.s, .s)
  diag(.cor) <- 1
  .sd[etas] <- sqrt(ini[etas])
  .cor * outer(.sd, .sd)
}

#' Keep a probe start strictly inside the first fit's bounds
#'
#' @param theta the start
#' @param lower,upper the bounds (may be infinite)
#' @return `theta`, with a value on or past a bound moved just inside it (by at
#'   most a quarter of the interval, so a narrow one stays feasible); values
#'   already inside are left as they are
#' @noRd
.foceiCollapseProbeClamp <- function(theta, lower, upper) {
  .eps <- pmin((.Machine$double.eps)^(1 / 7), (upper - lower) / 4)
  ifelse(theta <= lower, lower + .eps, ifelse(theta >= upper, upper - .eps, theta))
}

#' Re-fit from a collapsed omega raised back to its initial estimate
#'
#' With `collapseProbe=TRUE`, a diagonal omega the outer optimizer drove toward
#' zero (or onto its lower bound) is raised back to its initial estimate, the
#' other parameters kept at their final values, and the problem is re-fit from
#' there; the fit with the lower objective is kept (#1394).
#' @param fit the finished FOCEi environment (or a `try-error`)
#' @param snap the inputs from `.foceiCollapseProbeSnapshot()`
#' @param ui the model before estimation
#' @param control the fit's control
#' @return the kept fit environment
#' @noRd
.foceiCollapseProbe <- function(fit, snap, ui, control) {
  if (is.null(snap) || !is.environment(fit) || !is.matrix(fit$omega)) {
    return(fit)
  }
  .ini <- .omegaCollapseIni(ui$iniDf)
  .om <- fit$omega
  .bad <- .omegaCollapsed(
    .ini,
    setNames(diag(.om), rownames(.om)),
    .omegaCollapseWarnLimit(list(control = control)),
    .omegaCollapseBoundLower(list(control = control)),
    tryCatch(ui$omega, error = function(e) NULL),
    .om
  )
  .etas <- c(.bad$collapsed, .bad$bound)
  if (length(.etas) == 0L) {
    return(fit)
  }
  .same <- ui$omegaSameMap
  if (length(.same) != nrow(.om)) {
    .same <- NULL
  }
  .inv <- tryCatch(
    .foceiSymInvCholCreate(
      .foceiCollapseProbeOmega(.om, .ini, .etas, tryCatch(ui$omega, error = function(e) NULL)),
      rxode2::rxGetControl(ui, "diagXform", "sqrt"),
      .same,
      warn = FALSE,
      fallback = FALSE
    )$rxInv,
    error = function(e) NULL
  )
  .nth <- length(snap$thetaIni)
  if (is.null(.inv) || length(.inv$theta) != length(snap$rxInv$theta) || length(fit$fullTheta) < .nth) {
    return(fit)
  }
  .probe <- list2env(snap, envir = new.env(parent = parent.env(fit)))
  .probe$thetaIni <- setNames(fit$fullTheta[seq_len(.nth)], names(snap$thetaIni))
  # keep the omega start inside the bounds the first fit was given
  .k <- .nth + seq_along(.inv$theta)
  .inv$theta <- .foceiCollapseProbeClamp(.inv$theta, snap$lower[.k], snap$upper[.k])
  .probe$rxInv <- .inv
  .minfo(paste0("probing collapsed omega: ", paste(.etas, collapse = ", ")))
  .fit2 <- try(.foceiFitInternal(.probe), silent = TRUE)
  .d <- if (is.environment(.fit2)) fit$objective - .fit2$objective else NA_real_
  if (isTRUE(.d > 0)) {
    warning(
      .runInfoNameMsg(sprintf("omega collapse probe lowered objective by %s: ", signif(.d, 3)), .etas),
      call. = FALSE
    )
    return(.fit2)
  }
  warning(.runInfoNameMsg("omega collapse probe kept the original fit: ", .etas), call. = FALSE)
  fit
}
