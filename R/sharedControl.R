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
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("foceiControl", .ctl)
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
getValidNlmixrCtl.foce <- getValidNlmixrCtl.focei

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.fo <- getValidNlmixrCtl.focei

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.foi <- getValidNlmixrCtl.focei

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.posthoc <- getValidNlmixrCtl.focei

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.foce <- getValidNlmixrCtl.focei

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
  if (is.null(.ctl)) .ctl <- rxControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) {
    .ctl <- do.call(rxode2::rxControl, .ctl)
  }
  if (!inherits(.ctl, "rxControl")) {
    .ctl <- .ctl$rxControl
    if (!inherits(.ctl, "rxControl")) {
      .minfo(paste0("invalid control for `est=\"", class(control)[1], "\"`, using default"))
      .ctl <- rxode2::rxControl()
    } else {
      .ctl <- do.call(rxode2::rxControl, .ctl)
    }
  } else {
    .ctl <- do.call(rxode2::rxControl, .ctl)
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
  stop("do not know how to validate control for `est=\"", .cls, "\"`, please add `getValidNlmixrCtl.", .cls, "` method",
       call.=FALSE)
}

#'  Get specified control structure from reference
#'
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
