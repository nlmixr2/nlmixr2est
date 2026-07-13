#' Control options for the mufocei estimation method
#'
#' Mu-referenced-FOCEI-family closed-form-regression (`"lin"`) variant of
#' FOCEI; see `foceiControl(muModel=)`.
#'
#' @section Difference from `focei`:
#' The `mufocei`/`irlsfocei` (and related) methods apply the mu2+ covariate
#' hooks, which expand algebraic mu-referenced covariate expressions (e.g.
#' `cl.wt*log(WT/70)`) into estimable mu-referenced parameters and split
#' covariates into non-time-varying (absorbed into the phi term) and
#' time-varying (kept as `beta` regressors).  Calling `focei` directly does NOT
#' apply these hooks, so these methods can estimate more mu-referenced models
#' than plain `focei` -- there is a genuine difference between calling e.g.
#' `est="mufocei"` and `est="focei"`.
#'
#' All mu-referenced population thetas -- with or without covariates -- are
#' profiled out of the outer optimizer by the in-C++ regression
#' (intercept-only for covariate-free pairs), so outer gradients are only
#' calculated for the non-mu-referenced parameters (residual errors, omegas,
#' non-mu thetas).  Bounded or fixed mu-referenced thetas stay
#' outer-optimized.
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for `mufoceiControl()`
#'   this is always `"lin"` and cannot be changed -- use `irlsfoceiControl()`
#'   for the IRLS variant.
#' @return mufoceiControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' mufoceiControl()
mufoceiControl <- function(sigdig=3,
                           ...,
                           muModel=c("lin", "irls", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., muModel="lin")
  class(.control) <- "mufoceiControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.mufoceiControl <- function(control, env) {
  assign("mufoceiControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mufocei <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- mufoceiControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("mufoceiControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to mufoceiControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(mufoceiControl, .ctl)
  } else if (!inherits(.ctl, "mufoceiControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- mufoceiControl()
  } else {
    .ctl <- do.call(mufoceiControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mufocei <- function(x, ...) {
  .env <- x[[1]]
  if (exists("mufoceiControl", .env)) {
    .control <- get("mufoceiControl", .env)
    if (inherits(.control, "mufoceiControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "mufoceiControl")) return(.control)
  }
  stop("cannot find mufocei related control object", call.=FALSE)
}

.mufoceiControlToFoceiControl <- function(env, assign=TRUE) {
  .mufoceiControl <- env$mufoceiControl
  .n <- names(.mufoceiControl)
  .foceiControl <- setNames(lapply(.n, function(n) .mufoceiControl[[n]]), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mufocei <- function(x, ...) {
  .env <- x[[1]]
  .mufoceiControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the irlsfocei estimation method
#'
#' Mu-referenced-FOCEI-family reweighted-regression (`"irls"`) variant of
#' FOCEI; see `foceiControl(muModel=)`.
#'
#' @inheritSection mufoceiControl Difference from `focei`
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for `irlsfoceiControl()`
#'   this is always `"irls"` and cannot be changed -- use `mufoceiControl()`
#'   for the closed-form OLS variant.
#' @return irlsfoceiControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' irlsfoceiControl()
irlsfoceiControl <- function(sigdig=3,
                             ...,
                             muModel=c("irls", "lin", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., muModel="irls")
  class(.control) <- "irlsfoceiControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.irlsfoceiControl <- function(control, env) {
  assign("irlsfoceiControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.irlsfocei <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- irlsfoceiControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("irlsfoceiControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to irlsfoceiControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(irlsfoceiControl, .ctl)
  } else if (!inherits(.ctl, "irlsfoceiControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- irlsfoceiControl()
  } else {
    .ctl <- do.call(irlsfoceiControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.irlsfocei <- function(x, ...) {
  .env <- x[[1]]
  if (exists("irlsfoceiControl", .env)) {
    .control <- get("irlsfoceiControl", .env)
    if (inherits(.control, "irlsfoceiControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "irlsfoceiControl")) return(.control)
  }
  stop("cannot find irlsfocei related control object", call.=FALSE)
}

.irlsfoceiControlToFoceiControl <- function(env, assign=TRUE) {
  .irlsfoceiControl <- env$irlsfoceiControl
  .n <- names(.irlsfoceiControl)
  .foceiControl <- setNames(lapply(.n, function(n) .irlsfoceiControl[[n]]), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.irlsfocei <- function(x, ...) {
  .env <- x[[1]]
  .irlsfoceiControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the mufoce estimation method
#'
#' Mu-referenced-FOCEI-family closed-form-regression (`"lin"`) variant of
#' FOCE (no interaction); see `foceiControl(muModel=)`.
#'
#' @inheritSection mufoceiControl Difference from `focei`
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`; it cannot be changed, use `mufocei` instead
#' @param muModel Selects the regression variant; for `mufoceControl()`
#'   this is always `"lin"` and cannot be changed -- use `irlsfoceControl()`
#'   for the IRLS variant.
#' @return mufoceControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' mufoceControl()
mufoceControl <- function(sigdig=3,
                          ...,
                          interaction=FALSE,
                          muModel=c("lin", "irls", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., interaction=FALSE, muModel="lin")
  class(.control) <- "mufoceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.mufoceControl <- function(control, env) {
  assign("mufoceControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mufoce <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- mufoceControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("mufoceControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to mufoceControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(mufoceControl, .ctl)
  } else if (!inherits(.ctl, "mufoceControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- mufoceControl()
  } else {
    .ctl <- do.call(mufoceControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mufoce <- function(x, ...) {
  .env <- x[[1]]
  if (exists("mufoceControl", .env)) {
    .control <- get("mufoceControl", .env)
    if (inherits(.control, "mufoceControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "mufoceControl")) return(.control)
  }
  stop("cannot find mufoce related control object", call.=FALSE)
}

.mufoceControlToFoceiControl <- function(env, assign=TRUE) {
  .mufoceControl <- env$mufoceControl
  .n <- names(.mufoceControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.mufoceControl$interaction)
                                     }
                                     .mufoceControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mufoce <- function(x, ...) {
  .env <- x[[1]]
  .mufoceControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the irlsfoce estimation method
#'
#' Mu-referenced-FOCEI-family reweighted-regression (`"irls"`) variant of
#' FOCE (no interaction); see `foceiControl(muModel=)`.
#'
#' @inheritSection mufoceiControl Difference from `focei`
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction Interaction term for the model, in this case the
#'   default is `FALSE`; it cannot be changed, use `irlsfocei` instead
#' @param muModel Selects the regression variant; for `irlsfoceControl()`
#'   this is always `"irls"` and cannot be changed -- use `mufoceControl()`
#'   for the closed-form OLS variant.
#' @return irlsfoceControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' irlsfoceControl()
irlsfoceControl <- function(sigdig=3,
                            ...,
                            interaction=FALSE,
                            muModel=c("irls", "lin", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., interaction=FALSE, muModel="irls")
  class(.control) <- "irlsfoceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.irlsfoceControl <- function(control, env) {
  assign("irlsfoceControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.irlsfoce <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- irlsfoceControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("irlsfoceControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to irlsfoceControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(irlsfoceControl, .ctl)
  } else if (!inherits(.ctl, "irlsfoceControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- irlsfoceControl()
  } else {
    .ctl <- do.call(irlsfoceControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.irlsfoce <- function(x, ...) {
  .env <- x[[1]]
  if (exists("irlsfoceControl", .env)) {
    .control <- get("irlsfoceControl", .env)
    if (inherits(.control, "irlsfoceControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "irlsfoceControl")) return(.control)
  }
  stop("cannot find irlsfoce related control object", call.=FALSE)
}

.irlsfoceControlToFoceiControl <- function(env, assign=TRUE) {
  .irlsfoceControl <- env$irlsfoceControl
  .n <- names(.irlsfoceControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.irlsfoceControl$interaction)
                                     }
                                     .irlsfoceControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.irlsfoce <- function(x, ...) {
  .env <- x[[1]]
  .irlsfoceControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the muagq estimation method
#'
#' Mu-referenced-FOCEI-family closed-form-regression (`"lin"`) variant of
#' adaptive Gauss-Hermite quadrature; see `foceiControl(muModel=)`.
#'
#' @inheritParams agqControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for `muagqControl()`
#'   this is always `"lin"` and cannot be changed -- use `irlsagqControl()`
#'   for the IRLS variant.
#' @return muagqControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' muagqControl()
muagqControl <- function(sigdig=3, nAGQ=2, ..., interaction=TRUE,
                         agqLow=-Inf, agqHi=Inf,
                         muModel=c("lin", "irls", "none")) {
  .control <- foceiControl(sigdig=sigdig, ...,
                           nAGQ=nAGQ, interaction=interaction,
                           agqLow=agqLow, agqHi=agqHi, muModel="lin")
  class(.control) <- "muagqControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.muagqControl <- function(control, env) {
  assign("muagqControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.muagq <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- muagqControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("muagqControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to muagqControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(muagqControl, .ctl)
  } else if (!inherits(.ctl, "muagqControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- muagqControl()
  } else {
    .ctl <- do.call(muagqControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.muagq <- function(x, ...) {
  .env <- x[[1]]
  if (exists("muagqControl", .env)) {
    .control <- get("muagqControl", .env)
    if (inherits(.control, "muagqControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "muagqControl")) return(.control)
  }
  stop("cannot find muagq related control object", call.=FALSE)
}

.muagqControlToFoceiControl <- function(env, assign=TRUE) {
  .muagqControl <- env$muagqControl
  .n <- names(.muagqControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.muagqControl$interaction)
                                     }
                                     .muagqControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.muagq <- function(x, ...) {
  .env <- x[[1]]
  .muagqControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the irlsagq estimation method
#'
#' Mu-referenced-FOCEI-family reweighted-regression (`"irls"`) variant of
#' adaptive Gauss-Hermite quadrature; see `foceiControl(muModel=)`.
#'
#' @inheritParams agqControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for `irlsagqControl()`
#'   this is always `"irls"` and cannot be changed -- use `muagqControl()`
#'   for the closed-form OLS variant.
#' @return irlsagqControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' irlsagqControl()
irlsagqControl <- function(sigdig=3, nAGQ=2, ..., interaction=TRUE,
                           agqLow=-Inf, agqHi=Inf,
                           muModel=c("irls", "lin", "none")) {
  .control <- foceiControl(sigdig=sigdig, ...,
                           nAGQ=nAGQ, interaction=interaction,
                           agqLow=agqLow, agqHi=agqHi, muModel="irls")
  class(.control) <- "irlsagqControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.irlsagqControl <- function(control, env) {
  assign("irlsagqControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.irlsagq <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- irlsagqControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("irlsagqControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to irlsagqControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(irlsagqControl, .ctl)
  } else if (!inherits(.ctl, "irlsagqControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- irlsagqControl()
  } else {
    .ctl <- do.call(irlsagqControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.irlsagq <- function(x, ...) {
  .env <- x[[1]]
  if (exists("irlsagqControl", .env)) {
    .control <- get("irlsagqControl", .env)
    if (inherits(.control, "irlsagqControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "irlsagqControl")) return(.control)
  }
  stop("cannot find irlsagq related control object", call.=FALSE)
}

.irlsagqControlToFoceiControl <- function(env, assign=TRUE) {
  .irlsagqControl <- env$irlsagqControl
  .n <- names(.irlsagqControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.irlsagqControl$interaction)
                                     }
                                     .irlsagqControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.irlsagq <- function(x, ...) {
  .env <- x[[1]]
  .irlsagqControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the mulaplace estimation method
#'
#' Mu-referenced-FOCEI-family closed-form-regression (`"lin"`) variant of
#' the Laplace method (`nAGQ=1`); see `foceiControl(muModel=)`.
#'
#' @inheritParams laplaceControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for `mulaplaceControl()`
#'   this is always `"lin"` and cannot be changed -- use
#'   `irlslaplaceControl()` for the IRLS variant.
#' @return mulaplaceControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' mulaplaceControl()
mulaplaceControl <- function(sigdig=3, ..., nAGQ=1,
                             muModel=c("lin", "irls", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., nAGQ=nAGQ, muModel="lin")
  class(.control) <- "mulaplaceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.mulaplaceControl <- function(control, env) {
  assign("mulaplaceControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mulaplace <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- mulaplaceControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("mulaplaceControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to mulaplaceControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(mulaplaceControl, .ctl)
  } else if (!inherits(.ctl, "mulaplaceControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- mulaplaceControl()
  } else {
    .ctl <- do.call(mulaplaceControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mulaplace <- function(x, ...) {
  .env <- x[[1]]
  if (exists("mulaplaceControl", .env)) {
    .control <- get("mulaplaceControl", .env)
    if (inherits(.control, "mulaplaceControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "mulaplaceControl")) return(.control)
  }
  stop("cannot find mulaplace related control object", call.=FALSE)
}

.mulaplaceControlToFoceiControl <- function(env, assign=TRUE) {
  .mulaplaceControl <- env$mulaplaceControl
  .n <- names(.mulaplaceControl)
  .foceiControl <- setNames(lapply(.n, function(n) .mulaplaceControl[[n]]), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mulaplace <- function(x, ...) {
  .env <- x[[1]]
  .mulaplaceControlToFoceiControl(.env, assign=FALSE)
}

#' Control options for the irlslaplace estimation method
#'
#' Mu-referenced-FOCEI-family reweighted-regression (`"irls"`) variant of
#' the Laplace method (`nAGQ=1`); see `foceiControl(muModel=)`.
#'
#' @inheritParams laplaceControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param muModel Selects the regression variant; for
#'   `irlslaplaceControl()` this is always `"irls"` and cannot be changed
#'   -- use `mulaplaceControl()` for the closed-form OLS variant.
#' @return irlslaplaceControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' irlslaplaceControl()
irlslaplaceControl <- function(sigdig=3, ..., nAGQ=1,
                               muModel=c("irls", "lin", "none")) {
  .control <- foceiControl(sigdig=sigdig, ..., nAGQ=nAGQ, muModel="irls")
  class(.control) <- "irlslaplaceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.irlslaplaceControl <- function(control, env) {
  assign("irlslaplaceControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.irlslaplace <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- irlslaplaceControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("irlslaplaceControl", .ctl)
  if (inherits(.ctl, "foceiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to irlslaplaceControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(irlslaplaceControl, .ctl)
  } else if (!inherits(.ctl, "irlslaplaceControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- irlslaplaceControl()
  } else {
    .ctl <- do.call(irlslaplaceControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.irlslaplace <- function(x, ...) {
  .env <- x[[1]]
  if (exists("irlslaplaceControl", .env)) {
    .control <- get("irlslaplaceControl", .env)
    if (inherits(.control, "irlslaplaceControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "irlslaplaceControl")) return(.control)
  }
  stop("cannot find irlslaplace related control object", call.=FALSE)
}

.irlslaplaceControlToFoceiControl <- function(env, assign=TRUE) {
  .irlslaplaceControl <- env$irlslaplaceControl
  .n <- names(.irlslaplaceControl)
  .foceiControl <- setNames(lapply(.n, function(n) .irlslaplaceControl[[n]]), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.irlslaplace <- function(x, ...) {
  .env <- x[[1]]
  .irlslaplaceControlToFoceiControl(.env, assign=FALSE)
}
