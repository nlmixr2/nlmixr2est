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
  if (!inherits(.ctl, "foceiControl")) {
    .minfo("invalid control for `est=\"", .cls, "\"`, using default")
    .ctl <- foceiControl()
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
  if (!inherits(.ctl, "nlmeControl")) {
    .minfo("invalid control for `est=\"nlme\"`, using default")
    .ctl <- nlmeControl()
  }
  .ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.saem <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- saemControl()
  if (!inherits(.ctl, "saemControl")) {
    .minfo("invalid control for `est=\"saem\"`, using default")
    .ctl <- saemControl()
  }
  .ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.rxSolve <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- rxControl()
  if (!inherits(.ctl, "rxControl")) {
    .cls <- .ctl$rxControl
    if (!inherits(.ctl, "rxControl")) {
      .minfo("invalid control for `est=\"", .cls, "\"`, using default")
    }
  }
  .cls
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.simulation <- getValidNlmixrCtl.rxSolve

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.default <- function(control) {
  .cls <- class(control)[1]
  stop("do not know how to validate control for `est=\"", .cls, "\"`, please add `getValidNlmixrControl.", .cls, "` method",
       call.=FALSE)
}
