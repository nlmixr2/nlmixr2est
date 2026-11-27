#' Control options for the two-stage Stieltjes estimation method
#'
#' `est="stieltjes"` follows NONMEM's `$ESTIMATION STIELTJES`: a conditional
#' method (`stage1`) gives tentative estimates, then the population
#' parameters are re-estimated by minimizing the integrated likelihood on a
#' radial Stieltjes grid placed at each subject's EBE and eta Hessian at
#' those estimates and kept fixed (`agqCenter="start"`), so the second stage
#' does not move with the conditional etas.  The second stage evaluates the
#' EBEs and Hessians itself (FOCEi-style unless the options say otherwise).
#' A covariance recomputed later with `setCov()` places the grid at the
#' final estimates.
#'
#' @inheritParams agqControl
#' @param stage1 Conditional method of the first stage: `"focei"`, `"foce"`
#'   or `"laplace"`.
#' @param agqCenter Grid placement of the second stage (see
#'   [foceiControl()]); `"start"` keeps the grid where the second stage
#'   starts, at the first-stage estimates.
#' @param ... Other [foceiControl()] options, used by both stages.
#' @return stieltjesControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' stieltjesControl()
#'
#' # NONMEM's GRID=(5,2,.5,.95) with a FOCE first stage
#' stieltjesControl(stage1 = "foce", agqGrid = c(5, 2, 0.5, 0.95))
#'
stieltjesControl <- function(
  sigdig = 3,
  ...,
  stage1 = c("focei", "foce", "laplace"),
  interaction = TRUE,
  nAGQ = 2,
  agqRule = "stieltjes",
  agqGrid = c(nr = 3, ns = 2, r0 = 0.6, r1 = 0.9),
  agqCenter = c("start", "mode", "prior")
) {
  stage1 <- match.arg(stage1)
  agqCenter <- match.arg(agqCenter)
  .control <- foceiControl(
    sigdig = sigdig,
    ...,
    nAGQ = nAGQ,
    interaction = interaction,
    agqRule = agqRule,
    agqGrid = agqGrid,
    agqCenter = agqCenter
  )
  if (.control$nAGQ <= 1L) {
    stop("'stieltjesControl' needs nAGQ > 1", call. = FALSE)
  }
  .control$stage1 <- stage1
  class(.control) <- "stieltjesControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.stieltjesControl <- function(control, env) {
  assign("stieltjesControl", control, envir = env)
}

#' Convert a focei-family control to a `stieltjesControl()`
#'
#' A plain `foceiControl()` carries its own quadrature defaults (`nAGQ = 0`,
#' Gauss-Hermite on the mode), which would replace the method's; drop every
#' field still at its `foceiControl()` default.
#'
#' @inheritParams .foceiFamilyControlAs
#' @return a `stieltjesControl()`
#' @noRd
.stieltjesControlAs <- function(ctl, target) {
  if (!identical(class(ctl)[1], "foceiControl")) {
    return(.foceiFamilyControlAs(ctl, target))
  }
  .l <- unclass(ctl)
  .src <- unclass(foceiControl())
  .keep <- vapply(
    names(.l),
    function(.n) !(.n %in% names(.src)) || !isTRUE(all.equal(.l[[.n]], .src[[.n]])),
    logical(1)
  )
  do.call(target, .l[.keep])
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.stieltjes <- function(control) {
  .getValidCtl(
    control,
    "stieltjesControl",
    convert = .foceiFamilyControlConvertible,
    convertFun = .stieltjesControlAs
  )
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.stieltjes <- function(x, ...) .nmObjGetControlByClass(x, "stieltjesControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.stieltjes <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "stieltjesControl", assign = FALSE)
}

#' @export
rxUiDeparse.stieltjesControl <- function(object, var) {
  .rxUiDeparseFoceiControl(object, var, type = "stieltjesControl")
}

#' The first-stage control of a Stieltjes fit
#'
#' @param control a `stieltjesControl()`
#' @return the `stage1` method's control, without a covariance step
#' @noRd
.stieltjesStage1Control <- function(control) {
  .type <- paste0(control$stage1, "Control")
  .ctl <- .foceiFamilyControlAs(control, .type)
  .ctl$covMethod <- 0L
  .ctl
}

#' Run the first (conditional) stage on a copy of the estimation environment
#'
#' @param env the estimation environment
#' @param control the `stieltjesControl()`
#' @return list with the stage-1 fit and the warnings it raised
#' @noRd
.stieltjesStage1 <- function(env, control) {
  .env1 <- new.env(parent = emptyenv())
  for (.n in ls(env, all.names = TRUE)) {
    assign(.n, get(.n, envir = env), envir = .env1)
  }
  .env1$control <- .stieltjesStage1Control(control)
  .env1$est <- control$stage1
  class(.env1) <- c(control$stage1, "nlmixr2Est")
  .acc <- new.env(parent = emptyenv())
  .acc$warn <- character(0)
  .fit <- withCallingHandlers(
    get(paste0("nlmixr2Est.", control$stage1), mode = "function")(.env1),
    warning = function(w) {
      .acc$warn <- c(.acc$warn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = .fit, warnings = .acc$warn)
}

#' Put the stage-1 estimates into the model the second stage starts from
#'
#' The stage-1 fit reports bounded thetas back-transformed; the model holds
#' them on the internal (`rxBoundedTr.*`) scale.
#'
#' @param ui the estimation model
#' @param fit the stage-1 fit
#' @param transforms `ui$boundedTransforms` before stage 1 ran
#' @return nothing; `ui`'s `iniDf` is updated in place
#' @noRd
.stieltjesStage1Ini <- function(ui, fit, transforms = NULL) {
  .ini <- ui$iniDf
  .theta <- fit$theta
  for (.tr in transforms) {
    if (!(.tr$name %in% names(.theta))) {
      next
    }
    .v <- .theta[[.tr$name]]
    # keep the value inside its bounds, as .getBoundedParams() does
    if (is.finite(.tr$lower) && is.finite(.tr$upper)) {
      .eps <- (.tr$upper - .tr$lower) * 1e-6
      .v <- max(.tr$lower + .eps, min(.tr$upper - .eps, .v))
    } else if (is.finite(.tr$lower)) {
      .v <- max(.tr$lower + 1e-6, .v)
    } else if (is.finite(.tr$upper)) {
      .v <- min(.tr$upper - 1e-6, .v)
    }
    .theta[[.tr$internalName]] <- switch(
      .tr$type,
      logit = rxode2::logit(.v, .tr$lower, .tr$upper),
      lower_exp = log(.v - .tr$lower),
      upper_exp = log(.tr$upper - .v),
      .v
    )
  }
  .w <- which(is.na(.ini$neta1) & .ini$name %in% names(.theta))
  .ini$est[.w] <- .theta[.ini$name[.w]]
  # by eta name: the stage-1 omega need not list every eta of the model
  .omega <- fit$omega
  .etaNames <- .ini$name[which(!is.na(.ini$neta1) & .ini$neta1 == .ini$neta2)]
  .etaNames <- .etaNames[order(.ini$neta1[match(.etaNames, .ini$name)])]
  .w <- which(!is.na(.ini$neta1))
  for (.i in .w) {
    .e1 <- .etaNames[.ini$neta1[.i]]
    .e2 <- .etaNames[.ini$neta2[.i]]
    if (all(c(.e1, .e2) %in% rownames(.omega))) {
      .ini$est[.i] <- .omega[.e1, .e2]
    }
  }
  assign("iniDf", .ini, envir = ui)
  invisible(NULL)
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.stieltjes <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'stieltjes'", .var.name = .ui$modelName)
  .control <- env$control
  if (!inherits(.control, "stieltjesControl")) {
    .control <- getValidNlmixrCtl.stieltjes(list(.control))
    env$control <- .control
  }
  # stage 1 finishes a whole fit, which resets the pre-process hooks' state
  .state <- .preProcessHookState()
  .tr <- env$ui$boundedTransforms
  .s1 <- .stieltjesStage1(env, .control)
  .preProcessHookStateRestore(.state)
  if (!is.null(.tr)) {
    env$ui$boundedTransforms <- .tr
  }
  .stieltjesStage1Ini(env$ui, .s1$fit, .tr)
  .foceiFamilyControl(env, ..., type = "stieltjesControl")
  .foceiFamilyControlToFoceiControl(env, "stieltjesControl")
  on.exit({
    if (exists("control", envir = .ui)) {
      rm("control", envir = .ui)
    }
  })
  env$stieltjesControl <- .control
  env$est <- "stieltjes"
  .ui <- env$ui
  .ret <- .foceiFamilyReturn(env, .ui, ..., est = "stieltjes")
  if (inherits(.ret, "nlmixr2FitCore")) {
    assign(
      "stieltjesStage1",
      list(
        est = .control$stage1,
        objective = .s1$fit$objective,
        theta = .s1$fit$theta,
        omega = .s1$fit$omega,
        warnings = .s1$warnings
      ),
      envir = .ret$env
    )
  }
  .ret
}
attr(nlmixr2Est.stieltjes, "nlmixr2LikPenalty") <- "full"
attr(nlmixr2Est.stieltjes, "nlmixr2Priors") <- "general"
attr(nlmixr2Est.stieltjes, "covPresent") <- TRUE
attr(nlmixr2Est.stieltjes, "unbounded") <- .foUnbounded
