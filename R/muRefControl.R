#' Control options for the mfocei estimation method
#'
#' Mu-referenced-FOCEI-family closed-form-regression (`"lin"`) variant of
#' FOCEI; see `foceiControl(muModel=)`.
#'
#' @section Difference from `focei`:
#' The `mfocei`/`ifocei` (and related) methods apply the mu2+ covariate
#' hooks, which expand algebraic mu-referenced covariate expressions (e.g.
#' `cl.wt*log(WT/70)`) into estimable mu-referenced parameters and split
#' covariates into non-time-varying (absorbed into the phi term) and
#' time-varying (kept as `beta` regressors).  Calling `focei` directly does NOT
#' apply these hooks, so these methods can estimate more mu-referenced models
#' than plain `focei` -- there is a genuine difference between calling e.g.
#' `est="mfocei"` and `est="focei"`.
#'
#' All mu-referenced population thetas -- with or without covariates -- are
#' profiled out of the outer optimizer by the in-C++ regression
#' (intercept-only for covariate-free pairs), so outer gradients are only
#' calculated for the non-mu-referenced parameters (residual errors, omegas,
#' non-mu thetas).  Bounded mu-referenced parameters are regression-updated
#' with the update clamped to the bounds (a clamp is reported once as a fit
#' note); user-fixed (`fix()`) mu thetas stay out of the regression.
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for `mfoceiControl()`
#'   this is always `"lin"` and cannot be changed -- use `ifoceiControl()`
#'   for the IRLS variant.
#' @return mfoceiControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' mfoceiControl()
mfoceiControl <- function(sigdig=3,
                           ...,
                           muModel=c("lin", "irls", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., muModel="lin")
  class(.control) <- "mfoceiControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.mfoceiControl <- function(control, env) {
  assign("mfoceiControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mfocei <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- mfoceiControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("mfoceiControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to mfoceiControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(mfoceiControl, .ctl)
  } else if (!inherits(.ctl, "mfoceiControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- mfoceiControl()
  } else {
    .ctl <- do.call(mfoceiControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mfocei <- function(x, ...) {
  .env <- x[[1]]
  if (exists("mfoceiControl", .env)) {
    .control <- get("mfoceiControl", .env)
    if (inherits(.control, "mfoceiControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "mfoceiControl")) return(.control)
  }
  stop("cannot find mfocei related control object", call.=FALSE)
}

.mfoceiControlToFoceiControl <- function(env, assign=TRUE) {
  .mfoceiControl <- env$mfoceiControl
  .n <- names(.mfoceiControl)
  .foceiControl <- setNames(lapply(.n, function(n) .mfoceiControl[[n]]), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mfocei <- function(x, ...) {
  .env <- x[[1]]
  .mfoceiControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the ifocei estimation method
#'
#' Mu-referenced-FOCEI-family reweighted-regression (`"irls"`) variant of
#' FOCEI; see `foceiControl(muModel=)`.
#'
#' @inheritSection mfoceiControl Difference from `focei`
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for `ifoceiControl()`
#'   this is always `"irls"` and cannot be changed -- use `mfoceiControl()`
#'   for the closed-form OLS variant.
#' @return ifoceiControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' ifoceiControl()
ifoceiControl <- function(sigdig=3,
                             ...,
                             muModel=c("irls", "lin", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., muModel="irls")
  class(.control) <- "ifoceiControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.ifoceiControl <- function(control, env) {
  assign("ifoceiControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.ifocei <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- ifoceiControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("ifoceiControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to ifoceiControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(ifoceiControl, .ctl)
  } else if (!inherits(.ctl, "ifoceiControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- ifoceiControl()
  } else {
    .ctl <- do.call(ifoceiControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.ifocei <- function(x, ...) {
  .env <- x[[1]]
  if (exists("ifoceiControl", .env)) {
    .control <- get("ifoceiControl", .env)
    if (inherits(.control, "ifoceiControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "ifoceiControl")) return(.control)
  }
  stop("cannot find ifocei related control object", call.=FALSE)
}

.ifoceiControlToFoceiControl <- function(env, assign=TRUE) {
  .ifoceiControl <- env$ifoceiControl
  .n <- names(.ifoceiControl)
  .foceiControl <- setNames(lapply(.n, function(n) .ifoceiControl[[n]]), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.ifocei <- function(x, ...) {
  .env <- x[[1]]
  .ifoceiControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the mfoce estimation method
#'
#' Mu-referenced-FOCEI-family closed-form-regression (`"lin"`) variant of
#' FOCE (no interaction); see `foceiControl(muModel=)`.
#'
#' @inheritSection mfoceiControl Difference from `focei`
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`; it cannot be changed, use `mfocei` instead
#' @param muModel Selects the regression variant; for `mfoceControl()`
#'   this is always `"lin"` and cannot be changed -- use `ifoceControl()`
#'   for the IRLS variant.
#' @return mfoceControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' mfoceControl()
mfoceControl <- function(sigdig=3,
                          ...,
                          interaction=FALSE,
                          muModel=c("lin", "irls", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., interaction=FALSE, muModel="lin")
  class(.control) <- "mfoceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.mfoceControl <- function(control, env) {
  assign("mfoceControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mfoce <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- mfoceControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("mfoceControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to mfoceControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(mfoceControl, .ctl)
  } else if (!inherits(.ctl, "mfoceControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- mfoceControl()
  } else {
    .ctl <- do.call(mfoceControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mfoce <- function(x, ...) {
  .env <- x[[1]]
  if (exists("mfoceControl", .env)) {
    .control <- get("mfoceControl", .env)
    if (inherits(.control, "mfoceControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "mfoceControl")) return(.control)
  }
  stop("cannot find mfoce related control object", call.=FALSE)
}

.mfoceControlToFoceiControl <- function(env, assign=TRUE) {
  .mfoceControl <- env$mfoceControl
  .n <- names(.mfoceControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.mfoceControl$interaction)
                                     }
                                     .mfoceControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mfoce <- function(x, ...) {
  .env <- x[[1]]
  .mfoceControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the ifoce estimation method
#'
#' Mu-referenced-FOCEI-family reweighted-regression (`"irls"`) variant of
#' FOCE (no interaction); see `foceiControl(muModel=)`.
#'
#' @inheritSection mfoceiControl Difference from `focei`
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`; it cannot be changed, use `ifocei` instead
#' @param muModel Selects the regression variant; for `ifoceControl()`
#'   this is always `"irls"` and cannot be changed -- use `mfoceControl()`
#'   for the closed-form OLS variant.
#' @return ifoceControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' ifoceControl()
ifoceControl <- function(sigdig=3,
                            ...,
                            interaction=FALSE,
                            muModel=c("irls", "lin", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., interaction=FALSE, muModel="irls")
  class(.control) <- "ifoceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.ifoceControl <- function(control, env) {
  assign("ifoceControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.ifoce <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- ifoceControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("ifoceControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to ifoceControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(ifoceControl, .ctl)
  } else if (!inherits(.ctl, "ifoceControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- ifoceControl()
  } else {
    .ctl <- do.call(ifoceControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.ifoce <- function(x, ...) {
  .env <- x[[1]]
  if (exists("ifoceControl", .env)) {
    .control <- get("ifoceControl", .env)
    if (inherits(.control, "ifoceControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "ifoceControl")) return(.control)
  }
  stop("cannot find ifoce related control object", call.=FALSE)
}

.ifoceControlToFoceiControl <- function(env, assign=TRUE) {
  .ifoceControl <- env$ifoceControl
  .n <- names(.ifoceControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.ifoceControl$interaction)
                                     }
                                     .ifoceControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.ifoce <- function(x, ...) {
  .env <- x[[1]]
  .ifoceControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the magq estimation method
#'
#' Mu-referenced-FOCEI-family closed-form-regression (`"lin"`) variant of
#' adaptive Gauss-Hermite quadrature; see `foceiControl(muModel=)`.
#'
#' @inheritParams agqControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for `magqControl()`
#'   this is always `"lin"` and cannot be changed -- use `iagqControl()`
#'   for the IRLS variant.
#' @return magqControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' magqControl()
magqControl <- function(sigdig=3, nAGQ=2, ..., interaction=TRUE,
                         agqLow=-Inf, agqHi=Inf,
                         muModel=c("lin", "irls", "none")) {
  .control <- foceiControl(sigdig=sigdig, ...,
                           nAGQ=nAGQ, interaction=interaction,
                           agqLow=agqLow, agqHi=agqHi, muModel="lin")
  class(.control) <- "magqControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.magqControl <- function(control, env) {
  assign("magqControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.magq <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- magqControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("magqControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to magqControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(magqControl, .ctl)
  } else if (!inherits(.ctl, "magqControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- magqControl()
  } else {
    .ctl <- do.call(magqControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.magq <- function(x, ...) {
  .env <- x[[1]]
  if (exists("magqControl", .env)) {
    .control <- get("magqControl", .env)
    if (inherits(.control, "magqControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "magqControl")) return(.control)
  }
  stop("cannot find magq related control object", call.=FALSE)
}

.magqControlToFoceiControl <- function(env, assign=TRUE) {
  .magqControl <- env$magqControl
  .n <- names(.magqControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.magqControl$interaction)
                                     }
                                     .magqControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.magq <- function(x, ...) {
  .env <- x[[1]]
  .magqControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the iagq estimation method
#'
#' Mu-referenced-FOCEI-family reweighted-regression (`"irls"`) variant of
#' adaptive Gauss-Hermite quadrature; see `foceiControl(muModel=)`.
#'
#' @inheritParams agqControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for `iagqControl()`
#'   this is always `"irls"` and cannot be changed -- use `magqControl()`
#'   for the closed-form OLS variant.
#' @return iagqControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' iagqControl()
iagqControl <- function(sigdig=3, nAGQ=2, ..., interaction=TRUE,
                           agqLow=-Inf, agqHi=Inf,
                           muModel=c("irls", "lin", "none")) {
  .control <- foceiControl(sigdig=sigdig, ...,
                           nAGQ=nAGQ, interaction=interaction,
                           agqLow=agqLow, agqHi=agqHi, muModel="irls")
  class(.control) <- "iagqControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.iagqControl <- function(control, env) {
  assign("iagqControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.iagq <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- iagqControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("iagqControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to iagqControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(iagqControl, .ctl)
  } else if (!inherits(.ctl, "iagqControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- iagqControl()
  } else {
    .ctl <- do.call(iagqControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.iagq <- function(x, ...) {
  .env <- x[[1]]
  if (exists("iagqControl", .env)) {
    .control <- get("iagqControl", .env)
    if (inherits(.control, "iagqControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "iagqControl")) return(.control)
  }
  stop("cannot find iagq related control object", call.=FALSE)
}

.iagqControlToFoceiControl <- function(env, assign=TRUE) {
  .iagqControl <- env$iagqControl
  .n <- names(.iagqControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.iagqControl$interaction)
                                     }
                                     .iagqControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.iagq <- function(x, ...) {
  .env <- x[[1]]
  .iagqControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the mlaplace estimation method
#'
#' Mu-referenced-FOCEI-family closed-form-regression (`"lin"`) variant of
#' the Laplace method (`nAGQ=1`); see `foceiControl(muModel=)`.
#'
#' @inheritParams laplaceControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for `mlaplaceControl()`
#'   this is always `"lin"` and cannot be changed -- use
#'   `ilaplaceControl()` for the IRLS variant.
#' @return mlaplaceControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' mlaplaceControl()
mlaplaceControl <- function(sigdig=3, ..., nAGQ=1,
                             muModel=c("lin", "irls", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., nAGQ=nAGQ, muModel="lin")
  class(.control) <- "mlaplaceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.mlaplaceControl <- function(control, env) {
  assign("mlaplaceControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mlaplace <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- mlaplaceControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("mlaplaceControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to mlaplaceControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(mlaplaceControl, .ctl)
  } else if (!inherits(.ctl, "mlaplaceControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- mlaplaceControl()
  } else {
    .ctl <- do.call(mlaplaceControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mlaplace <- function(x, ...) {
  .env <- x[[1]]
  if (exists("mlaplaceControl", .env)) {
    .control <- get("mlaplaceControl", .env)
    if (inherits(.control, "mlaplaceControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "mlaplaceControl")) return(.control)
  }
  stop("cannot find mlaplace related control object", call.=FALSE)
}

.mlaplaceControlToFoceiControl <- function(env, assign=TRUE) {
  .mlaplaceControl <- env$mlaplaceControl
  .n <- names(.mlaplaceControl)
  .foceiControl <- setNames(lapply(.n, function(n) .mlaplaceControl[[n]]), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mlaplace <- function(x, ...) {
  .env <- x[[1]]
  .mlaplaceControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the ilaplace estimation method
#'
#' Mu-referenced-FOCEI-family reweighted-regression (`"irls"`) variant of
#' the Laplace method (`nAGQ=1`); see `foceiControl(muModel=)`.
#'
#' @inheritParams laplaceControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for
#'   `ilaplaceControl()` this is always `"irls"` and cannot be changed
#'   -- use `mlaplaceControl()` for the closed-form OLS variant.
#' @return ilaplaceControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' ilaplaceControl()
ilaplaceControl <- function(sigdig=3, ..., nAGQ=1,
                               muModel=c("irls", "lin", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., nAGQ=nAGQ, muModel="irls")
  class(.control) <- "ilaplaceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.ilaplaceControl <- function(control, env) {
  assign("ilaplaceControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.ilaplace <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- ilaplaceControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("ilaplaceControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to ilaplaceControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(ilaplaceControl, .ctl)
  } else if (!inherits(.ctl, "ilaplaceControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- ilaplaceControl()
  } else {
    .ctl <- do.call(ilaplaceControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.ilaplace <- function(x, ...) {
  .env <- x[[1]]
  if (exists("ilaplaceControl", .env)) {
    .control <- get("ilaplaceControl", .env)
    if (inherits(.control, "ilaplaceControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "ilaplaceControl")) return(.control)
  }
  stop("cannot find ilaplace related control object", call.=FALSE)
}

.ilaplaceControlToFoceiControl <- function(env, assign=TRUE) {
  .ilaplaceControl <- env$ilaplaceControl
  .n <- names(.ilaplaceControl)
  .foceiControl <- setNames(lapply(.n, function(n) .ilaplaceControl[[n]]), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.ilaplace <- function(x, ...) {
  .env <- x[[1]]
  .ilaplaceControlToFoceiControl(.env, assign=FALSE)
}
