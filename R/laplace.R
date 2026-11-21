#' Control options for the Laplace estimation method
#'
#' This is the control options for the adaptive Gauss-Hermite
#' quadrature for the likelihood.  Note that nAGQ=1 is the same as the
#' Laplace method.
#'
#' The objective is NONMEM's `METHOD=COND LAPLACE INTERACTION` one (Wang
#' 2007, J Pharmacokinet Pharmacodyn 34:575, Eq. 15).  The inner Hessian
#' and the Laplace log-determinant use the exact second derivative of the
#' conditional -2 log-likelihood in the ETAs, not the Gauss-Newton FOCEi
#' curvature.  Each Gaussian endpoint is evaluated through its log-density,
#' exactly as adding `+ dnorm()` to it does, so it works on every model.
#' The FOCEi objective is available as `est = "focei"`; `est = "flaplace"`
#' computes the same Laplace objective from the analytic conditional
#' Hessian (faster, with an analytic outer gradient, but needs `fast = TRUE`
#' and second-order sensitivities).
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
laplaceControl <- function(sigdig = 3, ..., nAGQ = 1) {
  # interaction forces the calculation of the hessian, which is needed
  # for the adaptive Gaussian quadrature
  .control <- .laplaceLlikSet(foceiControl(sigdig = sigdig, ..., nAGQ = nAGQ), list(...))
  class(.control) <- "laplaceControl"
  .control
}

#' Mark a Laplace-family control to evaluate Gaussian endpoints by log-density
#'
#' NONMEM's `METHOD=COND LAPLACE INTERACTION` objective (Wang 2007, J
#' Pharmacokinet Pharmacodyn 34:575, Eq. 15) uses the exact second
#' derivative of the conditional -2 log-likelihood in eta at the mode.  The
#' Gaussian FOCEi inner problem uses the Gauss-Newton curvature instead
#' (Wang Eq. 16), so a Gaussian endpoint is instead evaluated through its
#' log-density (`dnorm`, exactly as `+ dnorm()` in the model does), whose
#' inner Hessian is the exact one.  This needs only first-order
#' sensitivities, so it applies to every model (`linCmt()`, `lag()`,
#' censoring, ...).
#'
#' The full analytic conditional Hessian (`innerHessian`/`detHessian =
#' "conditional"`, `est = "flaplace"`) is the same Laplace objective by a
#' different route, so the flag is off there.
#'
#' @param control a `foceiControl()` object
#' @param xtra the `...` the constructor was called with; an explicit
#'   `laplaceLlik` there is kept
#' @return `control` with `laplaceLlik` set
#' @noRd
.laplaceLlikSet <- function(control, xtra) {
  if (is.null(xtra$laplaceLlik)) {
    control$laplaceLlik <- !(control$innerHessian == "conditional" ||
      control$detHessian == "conditional")
  }
  control
}

#' Does this control ask for log-density (exact Hessian) Gaussian endpoints?
#'
#' @param control a focei-family control (or `NULL`)
#' @return logical
#' @noRd
.laplaceLlikOn <- function(control) {
  isTRUE(control$laplaceLlik) &&
    !identical(control$innerHessian, "conditional") &&
    !identical(control$detHessian, "conditional")
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.laplaceControl <- function(control, env) assign("laplaceControl", control, envir = env)

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.laplace <- function(control) {
  .getValidCtl(
    control,
    "laplaceControl",
    convert = c("foceiControl", "foceControl", "agqControl", "foControl", "foiControl"),
    convertFun = .foceiFamilyControlAs
  )
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.laplace <- function(x, ...) .nmObjGetControlByClass(x, "laplaceControl")

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.laplace <- function(x, ...) {
  .foceiFamilyControlToFoceiControl(x[[1]], "laplaceControl", assign = FALSE)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.laplace <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'laplace'", .var.name = .ui$modelName)
  .control <- env$control
  .foceiFamilyControl(env, ..., type = "laplaceControl")
  .foceiFamilyControlToFoceiControl(env, "laplaceControl")
  on.exit({
    if (exists("control", envir = .ui)) {
      rm("control", envir = .ui)
    }
  })
  env$laplaceControl <- .control
  env$est <- "laplace"
  .ui <- env$ui
  .foceiFamilyReturn(env, .ui, ..., est = "laplace")
}
attr(nlmixr2Est.laplace, "nlmixr2Priors") <- "general"
attr(nlmixr2Est.laplace, "iov") <- TRUE
attr(nlmixr2Est.laplace, "covPresent") <- TRUE
attr(nlmixr2Est.laplace, "unbounded") <- .foUnbounded

#' @export
rxUiDeparse.laplaceControl <- function(object, var) .rxUiDeparseFoceiControl(object, var, type = "laplaceControl")
