#' Control Values for nlme Fit with extra options for nlmixr
#'
#' The values supplied in the function call replace the defaults and
#' a list with all possible arguments is returned.  The returned list
#' is used as the ‘control’ argument to the ‘nlme’ function.
#'
#' @inheritParams nlme::nlmeControl
#' @inheritParams foceiControl
#' @export
#' @examples
#' nlmixr2::nlmeControl()
nlmeControl <- function(maxIter = 50, pnlsMaxIter = 7, msMaxIter = 50, minScale = 0.001,
    tolerance = 1e-05, niterEM = 25, pnlsTol = 0.001, msTol = 1e-06,
    returnObject = FALSE, msVerbose = FALSE, msWarnNoConv = TRUE,
    gradHess = TRUE, apVar = TRUE, .relStep = .Machine$double.eps^(1/3),
    minAbsParApVar = 0.05, opt = c("nlminb", "nlm"), natural = TRUE,
    sigma = NULL, optExpression=TRUE, sumProd=FALSE,
    rxControl=rxode2::rxControl(atol=1e-4, rtol=1e-4), ...) {

  checkmate::assertLogical(optExpression, len=1, any.missing=FALSE)
  checkmate::assertLogical(sumProd, len=1, any.missing=FALSE)
  checkmate::assertLogical(returnObject, len=1, any.missing=FALSE)
  checkmate::assertLogical(msVerbose, len=1, any.missing=FALSE)
  checkmate::assertLogical(msWarnNoConv, len=1, any.missing=FALSE)
  checkmate::assertLogical(gradHess, len=1, any.missing=FALSE)
  checkmate::assertLogical(apVar, len=1, any.missing=FALSE)
  checkmate::assertLogical(natural, len=1, any.missing=FALSE)

  checkmate::assertIntegerish(pnlsMaxIter, len=1, any.missing=FALSE, lower=1)
  checkmate::assertIntegerish(msMaxIter, len=1, any.missing=FALSE, lower=1)
  checkmate::assertIntegerish(niterEM, len=1, any.missing=FALSE, lower=1)

  checkmate::assertNumeric(minScale, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(pnlsTol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(msTol, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(tolerance, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(.relStep, len=1, any.missing=FALSE, lower=0)
  checkmate::assertNumeric(minAbsParApVar len=1, any.missing=FALSE, lower=0)

  if (!inherits(rxControl, "rxControl")) rxControl <- do.call("rbind", rxControl)

  if (is.null(sigma))
    sigma <- 0
  else if (!is.finite(sigma) || length(sigma) != 1 || sigma <
             0)
    stop("Within-group std. dev. must be a positive numeric value")
  .ret <- list(maxIter = maxIter, pnlsMaxIter = pnlsMaxIter, msMaxIter = msMaxIter,
               minScale = minScale, tolerance = tolerance, niterEM = niterEM,
               pnlsTol = pnlsTol, msTol = msTol, returnObject = returnObject,
               msVerbose = msVerbose, msWarnNoConv = msWarnNoConv, gradHess = gradHess,
               apVar = apVar, .relStep = .relStep, minAbsParApVar = minAbsParApVar,
               opt = match.arg(opt), natural = natural, sigma = sigma,
               optExpression=optExpression, sumProd=sumProd,
               rxControl=rxControl,
               ...)
  class(.ret) <- "nlmeControl"
  .ret
}

.nlmeFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- nlmixr2::nlmeControl()
  }
  if (!inherits(.control, "nlmeControl")){
    .control <- do.call(nlmixr2::saemControl, .control)
  }
  assign("control", .control, envir=.ui)

}
nlmixr2Est.nlme <- function(env, ...) {
  .nlmeFamilyControl(env, ...)
  on.exit({rm("control", envir=.ui)})
  .nlmeFamilyFit(env,  ...)
}


