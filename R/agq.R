#' Get the adaptive Gauss-Hermite quadrature points and weights
#'
#' @param neta number of eta parameters in the model
#' @param nAGQ number of adaptive quadrature points to use
#' @return A list with the following elements:
#' \item{x}{A matrix of quadrature points, one row per point.}
#' \item{w}{A matrix of quadrature weights, one row per point.}
#' \item{n}{The number of quadrature points.}
#' \item{neta}{The number of eta parameters.}
#' \item{nAQD}{The number of adaptive quadrature points.}
#' \item{first}{A logical indicating if the first point is zero.}
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
#' @examples
#' .agq(neta=2, nAGQ=3)
.agq <- function(neta=2, nAGQ=3) {
  .n <- .nlmixr2estAgq(NA_integer_)
  if (nAGQ > .n) {
    rxode2::rxReq("fastGHQuad") # conditionally require
    .gh <- fastGHQuad::gaussHermiteData(nAGQ)
    .x <- .gh$x
    .w <- .gh$w/sqrt(pi)
  } else {
    .gh <- .nlmixr2estAgq(as.integer(nAGQ))
    .x <- .gh$x
    .w <- .gh$w
  }
  .first <- FALSE
  if (nAGQ %% 2 == 1) {
    # If nAQD is odd, have the zero weight at the beginning so that
    # the F value is cached and it doesn't need to evaluate it twice
    .zero <- (nAGQ+1L)/2L
    .x <- c(0, .x[-.zero])
    .w <- c(.w[.zero], .w[-.zero])
    .first <- TRUE
  }
  .x <-   as.matrix(do.call("expand.grid", lapply(1:neta, function(x) .x)))
  .w <-   as.matrix(do.call("expand.grid", lapply(1:neta, function(x) .w)))
  list(
    x = .x,
    w = .w,
    n = nrow(.x),
    neta = neta,
    nAGQ = nAGQ,
    first = .first
  )
}

#' Control options for the agq estimation method
#'
#' This is the control options for the adaptive Gauss-Hermite
#' quadrature for the likelihood.  Note that nAGQ=1 is the same as the
#' Laplace method.
#'
#' @inheritParams foceiControl
#' @param ... Parameters used in the default `foceiControl()`
#' @param interaction boolean, Interaction term for the model, in this
#'   case the default is `TRUE`; For adaptive quadrature, with normal
#'   distribution the Hessian is calculated with the foce(i)
#'   approximation

#' @return agqControl object
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' \donttest{
#'
#' agqControl()
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
#' fit <- nlmixr(f, r, est="agq")
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
#' fit <- nlmixr(f, p, est="agq", control=agqControl(nAGQ=5))
#'
#' one.cmt <- function() {
#'   ini({
#'     ## You may label each parameter with a comment
#'     tka <- 0.45 # Log Ka
#'     tcl <- log(c(0, 2.7, 100)) # Log Cl
#'     ## This works with interactive models
#'     ## You may also label the preceding line with label("label text")
#'     tv <- 3.45; label("log V")
#'     ## the label("Label name") works with all models
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'     add.sd <- 0.7
#'  })
#'  model({
#'    ka <- exp(tka + eta.ka)
#'    cl <- exp(tcl + eta.cl)
#'    v <- exp(tv + eta.v)
#'    linCmt() ~ add(add.sd)
#'  })
#' }
#'
#' fit <- nlmixr(one.cmt, theo_sd, est="agq")
#'
#' }
#'
agqControl <- function(sigdig=3, nAGQ=2, ..., interaction=TRUE,
                       agqLow=-Inf,
                       agqHi=Inf) {
  # interaction forces the calculation of the hessian, which is needed
  # for the adaptive Gaussian quadrature
  .control <- foceiControl(sigdig=sigdig, ...,
                           nAGQ=nAGQ, interaction=interaction,
                           agqLow=agqLow,
                           agqHi=agqHi)
  class(.control) <- "agqControl"
  .control
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.agqControl <- function(control, env) {
  eval(rxode2::rxUiDeparse(control, "control"))
  assign("agqControl", control, envir=env)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.agq <- function(control) {
  .ctl <- control[[1]]
  .cls <- class(control)[1]
  if (is.null(.ctl)) .ctl <- agqControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list"))
    .ctl <- do.call("agqControl", .ctl)
  if (inherits(.ctl, "foceiControl") ||
        inherits(.ctl, "foceControl") ||
        inherits(.ctl, "foControl") ||
        inherits(.ctl, "foiControl")) {
    .minfo(paste0("converting ", class(.ctl)[1], " to agqControl"))
    class(.ctl) <- NULL
    .ctl <- do.call(agqControl, .ctl)
  } else if (!inherits(.ctl, "agqControl")) {
    .minfo(paste0("invalid control for `est=\"", .cls, "\"`, using default"))
    .ctl <- agqControl()
  } else {
    .ctl <- do.call(agqControl, .ctl)
  }
  .ctl
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.agq <- function(x, ...) {
  .env <- x[[1]]
  if (exists("agqControl", .env)) {
    .control <- get("agqControl", .env)
    if (inherits(.control, "agqControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "agqControl")) return(.control)
  }
  stop("cannot find agq related control object", call.=FALSE)
}

.agqControlToFoceiControl <- function(env, assign=TRUE) {
  .agqControl <- env$agqControl
  .ui <- env$ui
  .n <- names(.agqControl)
  .foceiControl <- setNames(lapply(.n,
                                   function(n) {
                                     if (n == "interaction") {
                                       return(.agqControl$interaction)
                                     }
                                     .agqControl[[n]]
                                   }), .n)
  class(.foceiControl) <- "foceiControl"
  if (assign) env$control <- .foceiControl
  .foceiControl
}

#' @rdname nmObjGetFoceiControl
#' @export
nmObjGetFoceiControl.agq <- function(x, ...) {
  .env <- x[[1]]
  .agqControlToFoceiControl(.env, assign=FALSE)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.agq <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiIovNoCor(.ui, " for the estimation routine 'agq'",
                             .var.name=.ui$modelName)
  .control <- env$control
  .uiApplyIov(env)
  .foceiFamilyControl(env, ..., type="agqControl")
  .agqControlToFoceiControl(env)
  on.exit({
    if (exists("control", envir=.ui)) {
      rm("control", envir=.ui)
    }
  })
  env$agqControl <- .control
  env$est <- "agq"
  .ui <- env$ui
  .uiFinalizeIov(.foceiFamilyReturn(env, .ui, ..., est="agq"))
}
attr(nlmixr2Est.agq, "covPresent") <- TRUE
