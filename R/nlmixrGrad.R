##' Get the optimal forward difference interval by Gill83 method
##'
##' @inheritParams base::do.call
##' @param which Which parameters to calculate the forward difference
##'     and optimal forward difference interval
##' @inheritParams foceiControl
##' @author Matthew Fidler
##'
##' @return
##'
##' A data frame with the following columns:
##' \itemize{
##' \item{info}{Gradient evaluation/forward difference information}
##' \item{hf}{Forward difference final estimate}
##' \item{df}{Derivative estimate}
##' \item{df2}{2nd Derivative Estimate}
##' \item{err}{Error of the final estimate derivative}
##' \item{aEps}{Absolute difference for forward numerical differences}
##' \item{rEps}{Relative Difference for backward numerical differences}
##' \item{aEpsC}{Absolute difference for central numerical differences}
##' \item{rEpsC}{Relative difference for central numerical differences}
##' }
##'
##' The \code{info} returns one of the following:
##' \itemize{
##' \item{Not Assessed}{Gradient wasn't assessed}
##' \item{Good}{Success in Estimating optimal forward difference interval}
##' \item{High Grad Error}{Large error; Derivative estimate error \code{fTol} or more of the derivative}
##' \item{Constant Grad}{Function constant or nearly constant for this parameter}
##' \item{Odd/Linear Grad}{Function odd or nearly linear, df = K, df2 ~ 0}
##' \item{Grad changes quickly}{df2 increases rapidly as h decreases}
##' }
##' @examples
##'
##' ## These are taken from the numDeriv::grad examples to show how
##' ## simple gradients are assessed with nlmixr2Gill83
##'
##' nlmixr2Gill83(sin, pi)
##'
##' nlmixr2Gill83(sin, (0:10)*2*pi/10)
##'
##' func0 <- function(x){ sum(sin(x))  }
##' nlmixr2Gill83(func0 , (0:10)*2*pi/10)
##'
##' func1 <- function(x){ sin(10*x) - exp(-x) }
##' curve(func1,from=0,to=5)
##'
##' x <- 2.04
##' numd1 <- nlmixr2Gill83(func1, x)
##' exact <- 10*cos(10*x) + exp(-x)
##' c(numd1$df, exact, (numd1$df - exact)/exact)
##'
##' x <- c(1:10)
##' numd1 <- nlmixr2Gill83(func1, x)
##' exact <- 10*cos(10*x) + exp(-x)
##' cbind(numd1=numd1$df, exact, err=(numd1$df - exact)/exact)
##'
##' sc2.f <- function(x){
##'   n <- length(x)
##'    sum((1:n) * (exp(x) - x)) / n
##' }
##'
##' sc2.g <- function(x){
##'   n <- length(x)
##'   (1:n) * (exp(x) - 1) / n
##' }
##'
##' x0 <- rnorm(100)
##' exact <- sc2.g(x0)
##'
##' g <- nlmixr2Gill83(sc2.f, x0)
##'
##' max(abs(exact - g$df)/(1 + abs(exact)))
##'
##' @export
nlmixr2Gill83 <- function(what, args, envir = parent.frame(),
                         which, gillRtol = sqrt(.Machine$double.eps), gillK = 10L, gillStep = 2, gillFtol = 0) {
  if (missing(which)) {
    which <- rep(TRUE, length(args))
  }
  return(nlmixr2Gill83_(what, args, envir, which,
    gillRtol = sqrt(.Machine$double.eps), gillK = 10L, gillStep = 2, gillFtol = 0
  ))
}

.nlmixr2GradInfo <- new.env(parent = emptyenv())
##' Create a gradient function based on gill numerical differences
##'
##' @param thetaNames Names for the theta parameters
##' @inheritParams nlmixr2Gill83
##' @inheritParams foceiControl
##' @param theta for the internal functions theta is the parameter
##'     values
##' @param md5 the md5 identifier for the internal gradient function
##'     information.
##'
##' @return A list with `eval`, `grad`, `hist` and `unscaled`
##'   functions.  This is an internal module used with dynmodel
##'
##' @keywords internal
##'
##' @examples
##'
##' func0 <- function(x){ sum(sin(x))  }
##'
##' ## This will printout every interation or when print=X
##' gf <- nlmixr2GradFun(func0)
##'
##' ## x
##' x <- (0:10)*2*pi/10;
##' gf$eval(x)
##' gf$grad(x)
##'
##' ## x2
##' x2 <- x+0.1
##' gf$eval(x2)
##' gf$grad(x2)
##'
##' ## Gives the parameter history as a data frame
##' gf$hist()
##'
##' @export
nlmixr2GradFun <- function(what, envir = parent.frame(), which, thetaNames,
                          gillRtol = sqrt(.Machine$double.eps), gillK = 10L, gillStep = 2, gillFtol = 0,
                          useColor = crayon::has_color(),
                          printNcol = floor((getOption("width") - 23) / 12),
                          print = 1) {
  .md5 <- digest::digest(list(what, gillRtol, gillK, gillStep, gillFtol))
  .nlmixr2GradInfo[["printNcol"]] <- printNcol
  .nlmixr2GradInfo[["useColor"]] <- useColor
  .nlmixr2GradInfo[["isRstudio"]] <- (Sys.getenv("RSTUDIO") == "1")
  .nlmixr2GradInfo[["print"]] <- print
  if (!missing(which)) {
    .nlmixr2GradInfo[[paste0(.md5, ".w")]] <- which
  }
  if (!missing(thetaNames)) {
    .nlmixr2GradInfo[["thetaNames"]] <- thetaNames
  }
  .nlmixr2GradInfo[[paste0(.md5, ".n")]] <- 0L
  .nlmixr2GradInfo[[paste0(.md5, ".f")]] <- what
  .nlmixr2GradInfo[[paste0(.md5, ".e")]] <- envir
  .nlmixr2GradInfo[[paste0(.md5, ".rtol")]] <- gillRtol
  .nlmixr2GradInfo[[paste0(.md5, ".k")]] <- gillK
  .nlmixr2GradInfo[[paste0(.md5, ".s")]] <- gillStep
  .nlmixr2GradInfo[[paste0(.md5, ".ftol")]] <- gillFtol
  .eval <- eval(parse(text = paste0("function(theta){
        nlmixr2Eval_(theta, \"", .md5, "\");
    }")))
  .grad <- eval(parse(text = paste0("function(theta){
        nlmixr2Grad_(theta, \"", .md5, "\");
    }")))
  .hist <- eval(parse(text = paste0("function(){
        nlmixr2ParHist_(md5=\"", .md5, "\");
    }")))
  .unscaled <- eval(parse(text = paste0("function(theta){
        nlmixr2Unscaled_(theta,md5=\"", .md5, "\");
    }")))
  return(list(eval = .eval, grad = .grad, hist = .hist, unscaled = .unscaled))
}

##' Calculate Hessian
##'
##' Unlike `stats::optimHess` which assumes the gradient is accurate,
##' nlmixr2Hess does not make as strong an assumption that the gradient
##' is accurate but takes more function evaluations to calculate the
##' Hessian.  In addition, this procedures optimizes the forward
##' difference interval by \code{\link{nlmixr2Gill83}}
##'
##' If you have an analytical gradient function, you should use
##' `stats::optimHess`
##'
##' @inheritParams stats::optimHess
##' @param ... Extra arguments sent to \code{\link{nlmixr2Gill83}}
##' @inheritParams base::do.call
##' @author Matthew Fidler
##' @return Hessian matrix based on Gill83
##' @export
##' @seealso \code{\link{nlmixr2Gill83}}, \code{\link{optimHess}}
##'
##' @examples
##'  func0 <- function(x){ sum(sin(x))  }
##'  x <- (0:10)*2*pi/10
##'  nlmixr2Hess(x, func0)
##'
##' fr <- function(x) {   ## Rosenbrock Banana function
##'     x1 <- x[1]
##'     x2 <- x[2]
##'     100 * (x2 - x1 * x1)^2 + (1 - x1)^2
##' }
##' grr <- function(x) { ## Gradient of 'fr'
##'     x1 <- x[1]
##'     x2 <- x[2]
##'     c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
##'        200 *      (x2 - x1 * x1))
##' }
##'
##' h1 <- optimHess(c(1.2,1.2), fr, grr)
##'
##' h2 <- optimHess(c(1.2,1.2), fr)
##'
##' ## in this case h3 is closer to h1 where the gradient is known
##'
##' h3 <- nlmixr2Hess(c(1.2,1.2), fr)
nlmixr2Hess <- function(par, fn, ..., envir = parent.frame()) {
  .gill <- nlmixr2Gill83(fn, par, envir = envir, ...)
  return(nlmixr2Hess_(par, fn, envir, .gill))
}
