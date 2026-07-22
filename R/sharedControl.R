#' Get valid nlmixr control object
#'
#' @param control nlmixr control object
#'
#' @param est Estimation routine
#'
#' @return Valid control object based on estimation method run.
#'
#' @details
#'
#' This is based on running the S3 method `getValidNlmixrCtl()` the
#'  `control` object is put into a list and the class of this new list
#'  is `c(est, "getValidNlmixrControl")`
#'
#' @export
getValidNlmixrControl <- function(control, est) {
  .ret <- list(control)
  class(.ret) <- c(est, "getValidNlmixrControl")
  getValidNlmixrCtl(.ret)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl <- function(control) {
  UseMethod("getValidNlmixrCtl")
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.focei <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- foceiControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("foceiControl", .ctl)
  if (!inherits(.ctl, "foceiControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- foceiControl()
  } else {
    .ctl <- do.call(foceiControl, .ctl)
  }
  .ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.nlme <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- nlmeControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("nlmeControl", .ctl)
  if (!inherits(.ctl, "nlmeControl")) {
    .minfo("invalid control for `est=\"nlme\"`, using default")
    .ctl <- nlmeControl()
  } else {
    .ctl <- do.call(nlmeControl, .ctl)
  }
  .ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.saem <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- saemControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("saemControl", .ctl)
  if (!inherits(.ctl, "saemControl")) {
    .minfo("invalid control for `est=\"saem\"`, using default")
    .ctl <- saemControl()
  } else {
    .ctl <- do.call(saemControl, .ctl)
  }
  .ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.rxSolve <- function(control) {
  .ctl <- control
  class(.ctl) <- NULL
  .ctl <- .ctl[[1]]
  .env <- nlmixr2global$nlmixrEvalEnv$envir
  if (!is.environment(.env)) {
    .env <- parent.frame(1)
  }
  if (is.null(.ctl)) .ctl <- rxControl(envir=.env)
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) {
    .ctl <- do.call(rxode2::rxControl, .ctl)
    .ctl$envir <- .env
  }
  if (!inherits(.ctl, "rxControl")) {
    .ctl <- .ctl$rxControl
    if (!inherits(.ctl, "rxControl")) {
      .minfo(paste0("invalid control for `est=\"", class(control)[1], "\"`, using default"))
      .ctl <- rxode2::rxControl(envir=.env)
    } else {
      .ctl <- do.call(rxode2::rxControl, .ctl)
      .ctl$envir <- .env
    }
  } else {
    .ctl <- do.call(rxode2::rxControl, .ctl)
    .ctl$envir <- .env
  }
  .ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.simulate <- getValidNlmixrCtl.rxSolve

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.simulation <- getValidNlmixrCtl.rxSolve

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.predict <- getValidNlmixrCtl.rxSolve



#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.tableControl <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- tableControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call(tableControl, .ctl)
  if (!inherits(.ctl, "tableControl")) {
    .minfo("invalid control for table, using default")
    .ctl <- tableControl()
  } else {
    .ctl <- do.call(tableControl, .ctl)
  }
  .ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.default <- function(control) {
  .cls <- class(control)[1]
  # An unknown est= reaches here before nlmixr2Est dispatch; show the tagged,
  # category-grouped list of available methods (issue #750) when one exists.
  .lines <- .nlmixr2EstTypeLines(current=.cls)
  if (length(.lines) > 0L) {
    stop("nlmixr2 estimation `est=\"", .cls, "\"` is not supported; available methods:\n",
         paste(.lines, collapse="\n"),
         call.=FALSE)
  }
  stop("do not know how to validate control for `est=\"", .cls, "\"`, please add `getValidNlmixrCtl.", .cls, "` method",
       call.=FALSE)
}

#'  Get specified control structure from reference
#'
#' @param ref Reference control
#' @param ... Other arguments for new control
#' @return List for new control object
#' @author Matthew L. Fidler
#' @noRd
.getControlFromDots <- function(ref, ...){
  .in <- list(...)
  .out <- vector(mode="list")
  for (.n in names(ref)) {
    .w <- which(names(.in) == .n)
    if (length(.w) == 1L) {
      .out[[.n]] <- .in[[.w]]
      .in <- .in[-.w]
    }
  }
  return(list(ctl=.out, rest=.in))
}

#' Optimizer convergence tolerance derived from `sigdig`
#'
#' `10^(-sigdig)` -- the same exponent `sigdig` sets for the ODE `rtol`, so the
#' optimizer converges to exactly the precision the solve can support (no chasing
#' below the solver noise floor).  Shared so every estimation method's optimizer
#' ties to `sigdig` the same way.
#' @param sigdig optimization significant digits
#' @return the tolerance
#' @noRd
.sigdigOptTol <- function(sigdig) 10^(-sigdig)

#' L-BFGS `factr` derived from `sigdig` (relative-f tolerance `10^-sigdig`)
#' @param sigdig optimization significant digits
#' @return the `factr` value (`tol / .Machine$double.eps`)
#' @noRd
.sigdigFactr <- function(sigdig) 10^(-sigdig) / .Machine$double.eps

#' Scale a tuned default tolerance by `sigdig` around `sigdig = 4`
#'
#' Keeps the method's tuned default at `sigdig = 4` and tightens/loosens it by one
#' order of magnitude per significant digit (`default * 10^(4 - sigdig)`).  Used
#' where a method's optimizer default should be preserved at the default `sigdig`
#' rather than replaced by the FOCEi formula (nlm, nls, nlme).
#' @param default the tolerance at `sigdig = 4`
#' @param sigdig optimization significant digits
#' @return the scaled tolerance
#' @noRd
.sigdigScale <- function(default, sigdig) default * 10^(4 - sigdig)

#' Scale ODE solver tolerances from the optimization `sigdig`
#'
#' The optimization `sigdig` sets the optimizer tolerances directly (`10^-sigdig`,
#' see [.sigdigOptTol()]); here it sets the ODE solver tolerances with ONE formula
#' -- the same for every solver (stiff, non-stiff, auto-switch) so the story is
#' simple and consistent with the optimizer:
#'   `rtol = 10^-sigdig`, `atol = 10^(-sigdig-3)`
#' The `rtol` exponent IS `sigdig`, and `atol` sits three orders below.  Sensitivity
#' solves match the main solve (the gradient/covariance are built from them);
#' steady-state solves run one order looser.  `tighten` shifts every exponent
#' down by that many orders for a method that needs a tighter solve than the
#' optimizer target (e.g. `est="nls"`, whose LM step is sensitive to solver noise).
#' This mirrors the mapping moved into `rxode2::rxControl(sigdig=)` (rxode2 PR
#' #1150); until that rxode2 release is the minimum dependency, nlmixr2est applies
#' it to every auto-built `rxControl`.
#' @param rxControl an `rxode2::rxControl()` object (modified and returned)
#' @param sigdig optimization significant digits; `NULL` leaves `rxControl` as-is
#' @param skip character vector of tolerance field names the user set explicitly
#'   (e.g. from a `rxControl = list(atol = ...)`); these are left untouched so an
#'   explicit `atol`/`rtol` overrides the `sigdig`-derived value
#' @param tighten extra orders of magnitude to tighten every tolerance (default 0)
#' @return `rxControl` with ODE solver tolerances set from `sigdig`
#' @noRd
.rxControlScaleSigdig <- function(rxControl, sigdig, skip = character(0), tighten = 0) {
  if (is.null(sigdig) || is.null(rxControl)) return(rxControl)
  .rtol <- 10^(-(sigdig + tighten))
  .atol <- 10^(-(sigdig + 3 + tighten))
  # only set a tolerance the user did not pass explicitly (skip).  Sensitivity
  # solves match the main solve: the outer gradient and covariance are built from
  # them, so a looser sens tolerance degrades analytic gradient/covariance
  # accuracy (and leaves the optimizer's gradient less accurate than its
  # objective).  Steady-state solves stay one order looser than the main solve.
  .set <- function(field, value) {
    if (!(field %in% skip)) rxControl[[field]] <<- rep_len(value, length(rxControl[[field]]))
  }
  .set("rtol", .rtol)
  .set("atol", .atol)
  .set("rtolSens", .rtol)
  .set("atolSens", .atol)
  .set("ssRtol", 10 * .rtol)
  .set("ssAtol", 10 * .atol)
  if (!is.null(rxControl$ssRtolSens)) .set("ssRtolSens", 10 * .rtol)
  if (!is.null(rxControl$ssAtolSens)) .set("ssAtolSens", 10 * .atol)
  rxControl
}
