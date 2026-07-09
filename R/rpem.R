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

  # Final estimate = mean over the converged iterations.
  .k <- min(control$collect, niter)
  .w <- (niter - .k + 1L):niter
  muHat <- colMeans(muTr[.w, , drop = FALSE])
  omHat <- colMeans(omTr[.w, , drop = FALSE])
  sdHat <- mean(sdTr[.w])

  # One final E-step at the converged estimates to compute per-subject EBEs
  # (posterior-mean etas, Eq 53): EBE_i = sum_j eta_ij * w_ij with the
  # self-normalized importance weights w_ij = softmax_j(log p_ij).
  base[.cl$muIdx + 1L] <- muHat
  base[.cl$addSdIdx + 1L] <- sdHat
  omegaHat <- diag(omHat, .cl$nEta)
  rxode2::rxSetSeed(control$seed)
  .fe <- rpemEstepK1Draw(.e, base, .cl$etaIdx, omegaHat, control$nGauss, control$cores)
  rpemFree()
  .nG <- control$nGauss
  .nsub <- length(.fe$logp) / .nG
  .etaM <- matrix(.fe$eta, ncol = .cl$nEta)
  ebe <- matrix(0, .nsub, .cl$nEta, dimnames = list(NULL, .cl$etaNames))
  for (.i in seq_len(.nsub)) {
    .idx <- ((.i - 1L) * .nG + 1L):(.i * .nG)
    .lw <- .fe$logp[.idx]; .wt <- exp(.lw - max(.lw)); .wt <- .wt / sum(.wt)
    ebe[.i, ] <- colSums(.etaM[.idx, , drop = FALSE] * .wt)
  }

  list(mu = stats::setNames(muHat, .cl$muNames),
       omega = stats::setNames(omHat, .cl$etaNames),
       addSd = sdHat, ebe = ebe,
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

#' Assemble a full nlmixr2FitData from RPEM estimates via eval-only FOCEI.
#'
#' Mirrors SAEM's finalize: set the estimated theta/omega and the RPEM EBEs on
#' the UI, then let FOCEI evaluate (0 outer/inner iterations) at those fixed
#' values to compute EBEs/residuals/tables (design/rpem/09).
#' @noRd
.rpemBuildFit <- function(env, ui, control, rfit) {
  .ret <- new.env(parent = emptyenv())
  .ret$table <- env$table
  .rxControl <- rxode2::rxControl(atol = control$atol, rtol = control$rtol,
                                  method = "liblsoda")
  .foceiPreProcessData(env$data, .ret, ui, .rxControl)
  .ret$ui <- ui
  .cl <- rfit$classify
  # full theta (named over all theta names): mu-referenced -> RPEM mu, additive
  # residual -> RPEM add.sd, held structural -> ini values.
  .tn <- .cl$thetaNames
  .ft <- stats::setNames(.cl$base[seq_along(.tn)], .tn)
  .ft[.cl$muNames] <- rfit$mu
  .ft[.tn[.cl$addSdIdx + 1L]] <- rfit$addSd
  .ret$fullTheta <- .ft
  # omega (diagonal) named over etas
  .om <- diag(rfit$omega, .cl$nEta)
  dimnames(.om) <- list(.cl$etaNames, .cl$etaNames)
  .ret$omega <- .om
  # per-subject EBEs -> etaMat for the FOCEI eval
  .eb <- rfit$ebe
  colnames(.eb) <- .cl$etaNames
  .ret$.etaMat <- .eb
  .ret$.etaMatBase <- .eb
  .ret$etaObf <- data.frame(ID = seq_len(nrow(.eb)),
                            stats::setNames(as.data.frame(.eb), .cl$etaNames),
                            OBJI = NA)
  .nlmixr2FitUpdateParams(.ret)
  .foceiControl <- foceiControl(maxOuterIterations = 0L, maxInnerIterations = 0L,
                                covMethod = 0L, etaMat = .eb, scaleTo = 0,
                                calcTables = TRUE, interaction = 1L,
                                rxControl = .rxControl, est = "rpem")
  .ret$control <- .foceiControl
  .ret$est <- "rpem"
  .ret$ofvType <- "rpem"
  .ret$adjObf <- TRUE
  # Store the control on the fit env so nmObjGetControl (hence methodOde, etc.)
  # resolve later.  nmObjHandleControlObject may null env$control, so pass the
  # saved .foceiControl variable (not .ret$control) to nlmixr2CreateOutputFromUi.
  nmObjHandleControlObject(.foceiControl, .ret)
  # Provide the compiled FOCEI model so the residual-table getters resolve:
  # nmObjGet.innerModel -> foceiModel$inner, nmObjGetIpredModel.default ->
  # foceiModel$predOnly (est="rpem" has no dedicated ipred-model getter). dataSav
  # is recomputed from origData on demand.
  .ret$foceiModel <- ui$focei
  .out <- nlmixr2CreateOutputFromUi(.ret$ui, data = .ret$origData,
                                    control = .foceiControl, table = .ret$table,
                                    env = .ret, est = "rpem")
  # The eval-only FOCEI pass produces the full nlmixr2FitData (parFixedDf, objDf,
  # omega, EBEs, shrinkage, and the per-observation IPRED/PRED/CWRES table).  If
  # anything goes wrong, require the full data.frame and let nlmixr2Est.rpem fall
  # back to the lightweight estimates object.
  if (!inherits(.out, "nlmixr2FitData")) {
    stop("rpem residual-table step incomplete")
  }
  .out
}

#' Retrieve the (focei eval) control stored on an RPEM fit.
#' @rdname nmObjGetControl
#' @export
nmObjGetControl.rpem <- function(x, ...) {
  .env <- x[[1]]
  for (.name in c("foceiControl0", "control")) {
    if (exists(.name, .env)) {
      .control <- get(.name, .env)
      if (inherits(.control, "foceiControl")) return(.control)
    }
  }
  stop("cannot find rpem control object", call.=FALSE)
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
  # Try the full nlmixr2FitData via eval-only FOCEI; fall back to the lightweight
  # estimates object if that path errors, so est="rpem" always returns something.
  .full <- tryCatch(.rpemBuildFit(env, .ui, .control, .fit),
                    error = function(e) {
                      .minfo(paste0("rpem: full fit object unavailable (",
                                    conditionMessage(e), "); returning estimates"))
                      NULL
                    })
  if (!is.null(.full)) return(.full)
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
