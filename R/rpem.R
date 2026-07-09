# RPEM estimation method (design/rpem/).  This file wires the validated C++ K=1
# engine (src/rpem.cpp) to a mu-referenced rxode2 UI: it classifies the model
# parameters and drives the E-step -> conjugate M-step iteration loop.
#
# M1 scope: single endpoint, diagonal omega, one additive residual (add.sd).
# Structural parameters without a between-subject eta are held fixed for now
# (their numeric M-step update is a follow-up); mixtures/IOV/censoring are later
# milestones.

#' Classify a mu-referenced UI into the pieces the RPEM engine needs.
#'
#' @param ui rxode2 UI object.
#' @return list with the engine inputs (see body).
#' @noRd
.rpemClassify <- function(ui) {
  .ini <- ui$iniDf
  .thetas <- .ini[!is.na(.ini$ntheta), , drop = FALSE]
  .thetas <- .thetas[order(.thetas$ntheta), , drop = FALSE]
  nTheta <- nrow(.thetas)
  .etas <- .ini[!is.na(.ini$neta1) & .ini$neta1 == .ini$neta2, , drop = FALSE]
  .etas <- .etas[order(.etas$neta1), , drop = FALSE]
  nEta <- nrow(.etas)
  if (nEta == 0L) stop("RPEM requires at least one between-subject random effect")
  if (any(.ini$fix)) stop("RPEM does not yet support fixed (fix()) parameters")
  # param vector order (rpemParams): THETA[1..nTheta], ETA[1..nEta] (+ DV from data)
  base <- c(.thetas$est, rep(0, nEta))
  etaIdx <- as.integer(nTheta + seq_len(nEta) - 1L)          # 0-based
  omega0 <- diag(.etas$est, nEta)
  # mu-referenced typical-value theta for each eta (in eta order)
  .mu <- ui$muRefDataFrame
  .muName <- .mu$theta[match(.etas$name, .mu$eta)]
  if (anyNA(.muName)) stop("RPEM requires every random effect to be mu-referenced")
  muIdx <- as.integer(match(.muName, .thetas$name) - 1L)     # 0-based theta positions
  # additive residual (M1)
  .res <- .thetas[!is.na(.thetas$err), , drop = FALSE]
  if (nrow(.res) != 1L || .res$err[1] != "add")
    stop("RPEM M1 currently supports exactly one additive residual (add.sd)")
  addSdIdx <- as.integer(match(.res$name, .thetas$name) - 1L)
  list(base = base, nTheta = nTheta, nEta = nEta, etaIdx = etaIdx, omega0 = omega0,
       muIdx = muIdx, mu0 = .thetas$est[muIdx + 1L],
       addSdIdx = addSdIdx, addSd0 = .res$est,
       thetaNames = .thetas$name, etaNames = .etas$name, muNames = .muName)
}

#' Fit a mu-referenced model with RPEM (K=1 core).
#'
#' @param ui rxode2 UI object.
#' @param data Data frame (with a DV column).
#' @param control `rpemControl()`.
#' @return list of estimates (`mu`, `omega`, `addSd`) plus per-iteration traces.
#' @noRd
.rpemFit <- function(ui, data, control = rpemControl()) {
  .cl <- .rpemClassify(ui)
  .m <- ui$rpemRxModel$predOnly
  .nm <- c(paste0("THETA[", seq_len(.cl$nTheta), "]"),
           paste0("ETA[", seq_len(.cl$nEta), "]"))
  .e <- new.env()
  .e$predOnly <- .m
  .e$rxControl <- rxode2::rxControl(atol = control$atol, rtol = control$rtol)
  .e$param <- stats::setNames(.cl$base, .nm)
  .e$data <- data

  base <- .cl$base; mu <- .cl$mu0; omega <- .cl$omega0; addSd <- .cl$addSd0
  niter <- control$niter
  muTr <- matrix(0, niter, .cl$nEta); omTr <- matrix(0, niter, .cl$nEta)
  sdTr <- numeric(niter); llTr <- numeric(niter)
  for (.it in seq_len(niter)) {
    base[.cl$muIdx + 1L] <- mu
    base[.cl$addSdIdx + 1L] <- addSd
    rxode2::rxSetSeed(control$seed + .it)
    .est <- rpemEstepK1Draw(.e, base, .cl$etaIdx, omega, control$nGauss, control$cores)
    .ms  <- rpemMstepK1(mu, addSd, control$nMH, control$mhBurn)
    mu <- .ms$mu; omega <- .ms$omega; addSd <- .ms$addSd
    muTr[.it, ] <- mu; omTr[.it, ] <- diag(omega); sdTr[.it] <- addSd; llTr[.it] <- .est$lnL
  }
  rpemFree()

  .k <- min(control$collect, niter)
  .w <- (niter - .k + 1L):niter
  list(mu = stats::setNames(colMeans(muTr[.w, , drop = FALSE]), .cl$muNames),
       omega = stats::setNames(colMeans(omTr[.w, , drop = FALSE]), .cl$etaNames),
       addSd = mean(sdTr[.w]),
       lnL = llTr, muTrace = muTr, omegaTrace = omTr, sdTrace = sdTr,
       classify = .cl)
}

#' Validate the RPEM control (est="rpem")
#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.rpem <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- rpemControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("rpemControl", .ctl)
  if (!inherits(.ctl, "rpemControl")) {
    .minfo("invalid control for `est=\"rpem\"`, using default")
    .ctl <- rpemControl()
  } else {
    .ctl <- do.call(rpemControl, .ctl)
  }
  .ctl
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.rpem <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (!inherits(.control, "rpemControl")) .control <- rpemControl()
  # M1: single-endpoint, mu-referenced, diagonal omega, additive residual.
  # .rpemClassify raises a clear error if the model is outside that scope.
  .fit <- .rpemFit(.ui, env$data, .control)
  .fit$ui <- .ui
  .fit$control <- .control
  class(.fit) <- "nlmixr2rpem"
  .fit
}

#' @export
print.nlmixr2rpem <- function(x, ...) {
  cat("RPEM fit (K=1 minimal core -- not yet a full nlmixr2 fit object)\n")
  cat("-- typical values (mu):\n"); print(round(x$mu, 4))
  cat("-- omega (between-subject variance):\n"); print(round(x$omega, 4))
  cat(sprintf("-- add.sd: %.4f\n", x$addSd))
  cat(sprintf("-- final lnL: %.3f over %d iterations\n",
              x$lnL[length(x$lnL)], length(x$lnL)))
  invisible(x)
}
