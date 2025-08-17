#' Control options for the Laplace estimation method
#'
#' This is the control options for the adaptive Gauss-Hermite
#' quadrature for the likelihood.  Note that nAGQ=1 is the same as the
#' Laplace method.
#'
#' This method can be made to more closely matches NONMEM-style Laplace
#' estimation by requesting the log-likelihood from STAN as well as
#' numerically calculated Hessian matrix. This is done with adding
#' `+dnorm()` to the model for any normal end-points.
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`

#' @return laplaceControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' \donttest{
#'
#' laplaceControl()
#'
#' # Use adaptive quadrature
#'
#' # x =  Litter size after 21 days, and the modeled value
#'
#' r <- rats
#' r$dv <- r$x
#'
#' # Time is not used in this model, but it is required in nlmixr2
#' # currently, add a dummy value
#'
#' r$time <- 0
#'
#' f <- function() {
#'   ini({
#'     t1 <- 1
#'     t2 <- 1
#'     t3 <- 1
#'     eta1 ~ 1
#'   })
#'   model({
#'     lp <- t1 * x1 + t2 * x2 + (x1 + x2*t3) * eta1
#'     p <- pnorm(lp)
#'     m1 <- m # need to add outside of model specification
#'     x ~ dbinom(m1, p)
#'   })
#' }
#'
#' fit <- nlmixr(f, r, est="laplace")
#'
#'
#' p <- pump
#'
#' p$dv <- p$y
#' p$time <- 0 # dummy time
#'
#' f <- function() {
#'   ini({
#'     t1 <- 1
#'     t2 <- 1
#'     t3 <- 1
#'     t4 <- 1
#'     eta1 ~ 1
#'   })
#'   model({
#'     if (group == 1) {
#'        lp <- t1 + t2 * logtstd
#'     } else {
#'        lp <- t3 + t4 * logtstd
#'     }
#'     lp <- lp + eta1
#'     lam <- exp(lp)
#'     y ~ dpois(lam)
#'   })
#' }
#'
#' fit <- nlmixr(f, p, est="laplace")
#'
#' }
#'
laplaceControl <- function(sigdig=3, ..., nAGQ=1) {
  # interaction forces the calculation of the hessian, which is needed
  # for the adaptive Gaussian quadrature
  .control <- foceiControl(sigdig=sigdig,
                           ...,
                           nAGQ=nAGQ)
  class(.control) <- "laplaceControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.laplaceControl <- function(control, env) {
  assign("laplaceControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.laplace <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- laplaceControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("laplaceControl", .ctl)
  if (inherits(.ctl, "foceiControl") ||
        inherits(.ctl, "foceControl") ||
        inherits(.ctl, "agqControl") ||
        inherits(.ctl, "foControl") ||
        inherits(.ctl, "foiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to laplaceControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(laplaceControl, .ctl)
  } else if (!inherits(.ctl, "laplaceControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- laplaceControl()
  } else {
    .ctl <- do.call(laplaceControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.laplace <- function(x, ...) {
  .env <- x[[1]]
  if (exists("laplaceControl", .env)) {
    .control <- get("laplaceControl", .env)
    if (inherits(.control, "laplaceControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "laplaceControl")) return(.control)
  }
  stop("cannot find laplace related control object", call.=FALSE)
}

.laplaceControlToFoceiControl <- function(env, assign=TRUE) {
  .laplaceControl <- env$laplaceControl
  .ui <- env$ui
  .n <- names(.laplaceControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.laplaceControl$interaction)
                                     }
                                     .laplaceControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.laplace <- function(x, ...) {
  .env <- x[[1]]
  .laplaceControlToFoceiControl(.env, assign=FALSE)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.laplace <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'laplace'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .uiApplyIov(env)
  .foceiFamilyControl(env, ..., type="laplaceControl")
  .laplaceControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$laplaceControl <- .control
  env$est <- "laplace"
  .ui <- env$ui
  .uiFinalizeIov(.foceiFamilyReturn(env, .ui, ..., est="laplace"))
}
attr(nlmixr2Est.laplace, "covPresent") <- TRUE
