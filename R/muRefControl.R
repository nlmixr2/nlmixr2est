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
mfoceiControl <- function(sigdig = 3, ..., muModel = c("lin", "irls", "none")) {
  .control <- foceiControl(sigdig = sigdig, ..., muModel = "lin")
  class(.control) <- "mfoceiControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.mfoceiControl <- function(control, env) assign("mfoceiControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mfocei <- function(control) .foceiFamilyValidCtl(control, "mfoceiControl")

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mfocei <- function(x, ...) .nmObjGetControlByClass(x, "mfoceiControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mfocei <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "mfoceiControl", assign = FALSE)
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
ifoceiControl <- function(sigdig = 3, ..., muModel = c("irls", "lin", "none")) {
  .control <- foceiControl(sigdig = sigdig, ..., muModel = "irls")
  class(.control) <- "ifoceiControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.ifoceiControl <- function(control, env) assign("ifoceiControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.ifocei <- function(control) .foceiFamilyValidCtl(control, "ifoceiControl")

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.ifocei <- function(x, ...) .nmObjGetControlByClass(x, "ifoceiControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.ifocei <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "ifoceiControl", assign = FALSE)
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
mfoceControl <- function(sigdig = 3, ..., interaction = FALSE, muModel = c("lin", "irls", "none")) {
  .control <- foceiControl(sigdig = sigdig, ..., interaction = FALSE, muModel = "lin")
  class(.control) <- "mfoceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.mfoceControl <- function(control, env) assign("mfoceControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mfoce <- function(control) .foceiFamilyValidCtl(control, "mfoceControl")

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mfoce <- function(x, ...) .nmObjGetControlByClass(x, "mfoceControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mfoce <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "mfoceControl", assign = FALSE)
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
ifoceControl <- function(sigdig = 3, ..., interaction = FALSE, muModel = c("irls", "lin", "none")) {
  .control <- foceiControl(sigdig = sigdig, ..., interaction = FALSE, muModel = "irls")
  class(.control) <- "ifoceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.ifoceControl <- function(control, env) assign("ifoceControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.ifoce <- function(control) .foceiFamilyValidCtl(control, "ifoceControl")

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.ifoce <- function(x, ...) .nmObjGetControlByClass(x, "ifoceControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.ifoce <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "ifoceControl", assign = FALSE)
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
magqControl <- function(
  sigdig = 3,
  nAGQ = 2,
  ...,
  interaction = TRUE,
  agqLow = -Inf,
  agqHi = Inf,
  muModel = c("lin", "irls", "none")
) {
  .control <- foceiControl(
    sigdig = sigdig,
    ...,
    nAGQ = nAGQ,
    interaction = interaction,
    agqLow = agqLow,
    agqHi = agqHi,
    muModel = "lin"
  )
  class(.control) <- "magqControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.magqControl <- function(control, env) assign("magqControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.magq <- function(control) .foceiFamilyValidCtl(control, "magqControl")

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.magq <- function(x, ...) .nmObjGetControlByClass(x, "magqControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.magq <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "magqControl", assign = FALSE)
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
iagqControl <- function(
  sigdig = 3,
  nAGQ = 2,
  ...,
  interaction = TRUE,
  agqLow = -Inf,
  agqHi = Inf,
  muModel = c("irls", "lin", "none")
) {
  .control <- foceiControl(
    sigdig = sigdig,
    ...,
    nAGQ = nAGQ,
    interaction = interaction,
    agqLow = agqLow,
    agqHi = agqHi,
    muModel = "irls"
  )
  class(.control) <- "iagqControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.iagqControl <- function(control, env) assign("iagqControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.iagq <- function(control) .foceiFamilyValidCtl(control, "iagqControl")

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.iagq <- function(x, ...) .nmObjGetControlByClass(x, "iagqControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.iagq <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "iagqControl", assign = FALSE)
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
mlaplaceControl <- function(sigdig = 3, ..., nAGQ = 1, muModel = c("lin", "irls", "none")) {
  .control <- foceiControl(sigdig = sigdig, ..., nAGQ = nAGQ, muModel = "lin")
  class(.control) <- "mlaplaceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.mlaplaceControl <- function(control, env) assign("mlaplaceControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.mlaplace <- function(control) .foceiFamilyValidCtl(control, "mlaplaceControl")

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.mlaplace <- function(x, ...) .nmObjGetControlByClass(x, "mlaplaceControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.mlaplace <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "mlaplaceControl", assign = FALSE)
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
ilaplaceControl <- function(sigdig = 3, ..., nAGQ = 1, muModel = c("irls", "lin", "none")) {
  .control <- foceiControl(sigdig = sigdig, ..., nAGQ = nAGQ, muModel = "irls")
  class(.control) <- "ilaplaceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.ilaplaceControl <- function(control, env) assign("ilaplaceControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.ilaplace <- function(control) .foceiFamilyValidCtl(control, "ilaplaceControl")

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.ilaplace <- function(x, ...) .nmObjGetControlByClass(x, "ilaplaceControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.ilaplace <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "ilaplaceControl", assign = FALSE)
}

#' focei-family controls that the mu-referenced methods convert
#'
#' Excludes `impmapControl`: its fields (`isample`, `nIter`, ...) are not in
#' `.foceiControlInternal`, so `foceiControl()` rejects them with "unused
#' argument".  `foControl`/`foiControl` are not converted either: a mu-referenced
#' method replaces them by its default control.
#' @noRd
.foceiFamilyControlConvertible <-
  c(
    "foceiControl",
    "foceControl",
    "focepControl",
    "agqControl",
    "laplaceControl",
    "mfoceiControl",
    "ifoceiControl",
    "mfoceControl",
    "ifoceControl",
    "mfocepControl",
    "ifocepControl",
    "magqControl",
    "iagqControl",
    "mlaplaceControl",
    "ilaplaceControl"
  )

#' `getValidNlmixrCtl()` of a mu-referenced focei-family method
#' @param control the list `getValidNlmixrControl()` dispatches on
#' @param ctl name of the method's control constructor
#' @return the validated control
#' @noRd
.foceiFamilyValidCtl <- function(control, ctl) {
  .getValidCtl(control, ctl, convert = .foceiFamilyControlConvertible, convertFun = .foceiFamilyControlAs)
}

#' Convert one focei-family control into another, keeping only what the caller set
#'
#' Each `*Control()` records its method's identity in ordinary fields (`interaction`,
#' `nAGQ`, `foce`, `fo`, `muModel`), so handing the whole object to another
#' constructor lets the SOURCE method's identity override the TARGET's -- a
#' `foceControl()` would run `est="mfocei"` as FOCE, a `foceControl()`'s `nAGQ=0`
#' would run `est="magq"` as plain FOCE.  Stripping a fixed list of fields instead
#' throws away deliberate overrides (`agqControl(nAGQ=5)`).
#'
#' Distinguish them by comparing against the SOURCE control's OWN defaults: a field
#' the caller changed is carried over, a field still at its default belongs to the
#' source method and is dropped so the target's value wins.
#'
#' The `posthoc` field of `foControl()`/`foiControl()` is an argument of those two
#' constructors only, so it is dropped for any other target.
#' @param ctl control object to convert
#' @param target name of the target `*Control()` function
#' @return a control of class `target`
#' @noRd
.foceiFamilyControlAs <- function(ctl, target) {
  .cls <- class(ctl)[1]
  .ctl <- unclass(ctl)
  if (!(target %in% c("foControl", "foiControl"))) {
    .ctl$posthoc <- NULL
  }
  if (identical(.cls, target)) {
    return(do.call(target, .ctl))
  }
  # foceiControl is the BASE control: its defaults are neutral and encode no method
  # identity of their own (interaction/nAGQ/foce are ordinary user settings there), so
  # every field passes through -- which is also exactly what shipped before this helper
  # existed.  Only the DERIVED family controls record their method in their defaults.
  if (identical(.cls, "foceiControl")) {
    return(do.call(target, .ctl))
  }
  .src <- tryCatch(unclass(do.call(.cls, list())), error = function(e) NULL)
  if (is.null(.src)) {
    return(do.call(target, .ctl))
  }
  .keep <- vapply(
    names(.ctl),
    function(.n) {
      if (!(.n %in% names(.src))) {
        return(TRUE)
      }
      !isTRUE(all.equal(.ctl[[.n]], .src[[.n]]))
    },
    logical(1)
  )
  do.call(target, .ctl[.keep])
}

#' The foceiControl that runs a focei-family method
#'
#' Every control of the family is a `foceiControl()` that carries its method's
#' settings under its own class, so the translation keeps each field and
#' relabels the class.  The `posthoc` field of fo/foi is not a `foceiControl()`
#' field: it is dropped, and `posthoc = FALSE` turns the inner problem off
#' (`maxInnerIterations = 0`).
#' @param env fit environment holding the method's control
#' @param ctl name of the control in `env` (e.g. `"foceControl"`)
#' @param assign when `TRUE`, also store the result as `env$control`
#' @param set named list of fields whose values are replaced
#' @return the `foceiControl` object
#' @noRd
.foceiFamilyControlToFoceiControl <- function(env, ctl, assign = TRUE, set = list()) {
  .ctl <- env[[ctl]]
  if (isFALSE(.ctl$posthoc)) {
    set$maxInnerIterations <- 0L
  }
  .n <- names(.ctl)
  .n <- .n[.n != "posthoc"]
  .ret <- setNames(lapply(.n, function(n) if (n %in% names(set)) set[[n]] else .ctl[[n]]), .n)
  class(.ret) <- "foceiControl"
  if (assign) {
    env$control <- .ret
  }
  .ret
}
