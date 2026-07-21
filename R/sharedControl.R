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
#' The FOCEi-family convention: a convergence tolerance is `10^(-sigdig-1)`.
#' Shared so every estimation method's optimizer ties to `sigdig` the same way
#' FOCEi's outer optimizer does.
#' @param sigdig optimization significant digits
#' @return the tolerance
#' @noRd
.sigdigOptTol <- function(sigdig) 10^(-sigdig - 1)

#' L-BFGS `factr` derived from `sigdig` (matches `foceiControl(lbfgsFactr=)`)
#' @param sigdig optimization significant digits
#' @return the `factr` value (`tol / .Machine$double.eps`)
#' @noRd
.sigdigFactr <- function(sigdig) 10^(-sigdig - 1) / .Machine$double.eps

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

#' Is an ODE solver method purely non-stiff?
#'
#' Uses `rxode2::rxIsNonStiff()` when the installed rxode2 exports it; otherwise
#' falls back to the same comparison (switcher and stiff-only method codes are NOT
#' non-stiff).  The fallback code list must stay in sync with `rxode2::rxIsNonStiff`.
#' @param method integer method code (as stored on an `rxControl`)
#' @return `TRUE` if the method is a purely non-stiff solver
#' @noRd
.rxIsNonStiffCode <- function(method) {
  .m <- as.integer(method)
  if (exists("rxIsNonStiff", envir = asNamespace("rxode2"), inherits = FALSE)) {
    return(isTRUE(rxode2::rxIsNonStiff(.m)))
  }
  # fallback (older rxode2): switchers (lsoda/liblsoda/indLin) and stiff-only
  # solvers are not non-stiff
  .switcher <- c(1L, 2L, 3L)
  .stiff <- c(13L, 14L, 21L, 31L, 32L, 33L, 34L, 35L, 36L, 37L, 38L, 107L,
              213L, 236L, 233L, 234L, 238L, 235L, 231L, 232L, 237L, 221L)
  !(.m %in% c(.switcher, .stiff))
}

#' Scale ODE solver tolerances from the optimization `sigdig`, split by stiffness
#'
#' The optimization `sigdig` sets the optimizer tolerances directly (unchanged);
#' here it additionally sets the ODE solver tolerances.  This mirrors the mapping
#' moved into `rxode2::rxControl(sigdig=)` itself (rxode2 PR #1150) -- until that
#' rxode2 release is the minimum dependency, nlmixr2est applies it to every
#' auto-built `rxControl` so all estimation methods share one convention.
#'
#' A stiff solver needs tighter tolerances than a purely non-stiff one, and `atol`
#' sits well below `rtol` in both.  Classification uses `rxode2::rxIsNonStiff()`
#' when available (else an equivalent fallback, see [.rxIsNonStiffCode()]): the
#' non-stiff (looser) tolerances apply only to a purely non-stiff solver
#' (`dop853`, `dop5`, the Runge-Kutta / Verner methods, ...) that is not part of an
#' auto-switch composite; every stiff-only solver (`ros4`, `cvode`, `bdf`, ...),
#' every auto-switching solver (`lsoda`/`liblsoda`), and every stiff-secondary
#' composite (`dop853+ros4`, ...) gets the tighter stiff tolerances.  Sensitivity
#' and steady-state solves run one order looser.
#'
#' At `sigdig = S`:
#'   non-stiff: rtol = 10^-S,     atol = 10^(-S-3)
#'   stiff:     rtol = 10^(-S-3), atol = 10^(-S-5)
#'   sensitivity / steady-state: 10x the corresponding main tolerance
#' @param rxControl an `rxode2::rxControl()` object (modified and returned)
#' @param sigdig optimization significant digits; `NULL` leaves `rxControl` as-is
#' @param skip character vector of tolerance field names the user set explicitly
#'   (e.g. from a `rxControl = list(atol = ...)`); these are left untouched so an
#'   explicit `atol`/`rtol` overrides the `sigdig`-derived value
#' @return `rxControl` with ODE solver tolerances set from `sigdig`
#' @noRd
.rxControlScaleSigdig <- function(rxControl, sigdig, skip = character(0)) {
  if (is.null(sigdig) || is.null(rxControl)) return(rxControl)
  # non-stiff only for a purely non-stiff base solver with no stiff auto-switch
  # secondary; switchers (lsoda/liblsoda), stiff-only solvers, and composites are
  # all stiff and need the tighter tolerances
  .nonStiff <- isTRUE(.rxIsNonStiffCode(rxControl$method)) &&
    (is.null(rxControl$stiff2) || identical(as.integer(rxControl$stiff2), 0L))
  .rtol <- 10^(-(if (.nonStiff) sigdig else sigdig + 3))
  .atol <- 10^(-(if (.nonStiff) sigdig + 3 else sigdig + 5))
  # only set a tolerance the user did not pass explicitly (skip); sensitivity +
  # steady-state solves run one order looser than the main solve
  .set <- function(field, value) {
    if (!(field %in% skip)) rxControl[[field]] <<- rep_len(value, length(rxControl[[field]]))
  }
  .set("rtol", .rtol)
  .set("atol", .atol)
  .set("rtolSens", 10 * .rtol)
  .set("atolSens", 10 * .atol)
  .set("ssRtol", 10 * .rtol)
  .set("ssAtol", 10 * .atol)
  if (!is.null(rxControl$ssRtolSens)) .set("ssRtolSens", 10 * .rtol)
  if (!is.null(rxControl$ssAtolSens)) .set("ssAtolSens", 10 * .atol)
  rxControl
}
