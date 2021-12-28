#' nlmixr2 fits population PK and PKPD non-linear mixed effects models.
#'
#' nlmixr2 is an R package for fitting population pharmacokinetic (PK)
#' and pharmacokinetic-pharmacodynamic (PKPD) models.
#'
#' The nlmixr2 generalized function allows common access to the nlmixr2
#' estimation routines.
#'
#' @template uif
#'
#' @param object Fitted object or function specifying the model.
#' @inheritParams nlmixr2_fit
#' @param ... Other parameters
#' @param save Boolean to save a nlmixr2 object in a rds file in the
#'     working directory.  If \code{NULL}, uses option "nlmixr2.save"
#' @return Either a nlmixr2 model or a nlmixr2 fit object
#' @author Matthew L. Fidler
#' @examples
#'
#' \donttest{
#'
#' f_ode <- function(){
#'     ini({
#'         lCl <- 1.6      #log Cl (Lhr)
#'         lVc <- log(80)   #log Vc (L)
#'         lKA <- 0.3      #log Ka (1/hr)
#'         prop.err <- c(0, 0.2, 1)
#'         eta.Cl ~ 0.3 ## BSV Cl
#'         eta.Vc ~ 0.2 ## BSV Vc
#'         eta.KA ~ 0.1 ## BSV Ka
#'     })
#'     model({
#'         ## First parameters are defined in terms of the initial estimates
#'         ## parameter names.
#'         Cl <- exp(lCl + eta.Cl)
#'         Vc = exp(lVc + eta.Vc)
#'         KA <- exp(lKA + eta.KA)
#'         ## After the differential equations are defined
#'         kel <- Cl / Vc;
#'         d/dt(depot)    = -KA*depot;
#'         d/dt(centr)  =  KA*depot-kel*centr;
#'         ## And the concentration is then calculated
#'         cp = centr / Vc;
#'         ## Last, nlmixr2 is told that the plasma concentration follows
#'         ## a proportional error (estimated by the parameter prop.err)
#'         cp ~ prop(prop.err)
#'     })
#' }
#' f_linCmt <- function(){
#'     ini({
#'         lCl <- 1.6      #log Cl (L/hr)
#'         lVc <- log(90)   #log Vc (L)
#'         lKA <- 0.1      #log Ka (1/hr)
#'         prop.err <- c(0, 0.2, 1)
#'         add.err <- c(0, 0.01)
#'         eta.Cl ~ 0.1 ## BSV Cl
#'         eta.Vc ~ 0.1 ## BSV Vc
#'         eta.KA ~ 0.1 ## BSV Ka
#'     })
#'     model({
#'         Cl <- exp(lCl + eta.Cl)
#'         Vc = exp(lVc + eta.Vc)
#'         KA <- exp(lKA + eta.KA)
#'         ## Instead of specifying the ODEs, you can use
#'         ## the linCmt() function to use the solved system.
#'         ##
#'         ## This function determines the type of PK solved system
#'         ## to use by the parameters that are defined.  In this case
#'         ## it knows that this is a one-compartment model with first-order
#'         ## absorption.
#'         linCmt() ~ add(add.err) + prop(prop.err)
#'     })
#' }
#'
#' # Use nlme algorithm
#' fit_linCmt_nlme <- try(nlmixr2(f_ode, Oral_1CPT, est="nlme",
#'                control=nlmeControl(maxstepsOde = 50000, pnlsTol=0.4)))
#' if (!inherits(fit_linCmt_nlme, "try-error")) print(fit_linCmt_nlme)
#'
#' # Use Focei algorithm
#' fit_linCmt_focei <- try(nlmixr2(f_linCmt, Oral_1CPT, est="focei"))
#' if (!inherits(fit_linCmt_focei, "try-error")) print(fit_linCmt_focei)
#'
#' # The ODE model can be fitted using the saem algorithm, more
#' # iterations should be used for real applications
#'
#' fit_ode_saem <- try(nlmixr2(f_ode, Oral_1CPT, est = "saem",
#'         control = saemControl(n.burn = 50, n.em = 100, print = 50)))
#' if (!inherits(fit_ode_saem, "try-error")) print(fit_ode_saem)
#'
#' }
#' @export
nlmixr2 <- function(object, data, est = NULL, control = list(),
                   table = tableControl(), ..., save = NULL,
                   envir = parent.frame()) {
  assignInMyNamespace(".nlmixr2Time", proc.time())
  rxode2::rxSolveFree()
  rxode2::.setWarnIdSort(FALSE)
  on.exit(rxode2::.setWarnIdSort(TRUE))
  force(est)
  ## verbose?
  ## https://tidymodels.github.io/model-implementation-principles/general-conventions.html
  UseMethod("nlmixr2")
}


#' @rdname nlmixr
#' @export
nlmixr <- nlmixr2


##' @rdname nlmixr2
##' @export
nlmixr2.function <- function(object, data, est = NULL, control = NULL, table = tableControl(), ...,
                            save = NULL, envir = parent.frame()) {
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .modName <- deparse(substitute(object))
  .uif <- rxode2::rxode(object)
  if (missing(data) && missing(est)) {
    return(.uif)
  } else if (missing(data)) {
    stop("need data", call.=FALSE)
  }
  .env <- new.env(parent=emptyenv())
  .env$ui <- .uif
  .env$data <- data
  .env$control <- control
  .env$table <- table
  class(.env) <- c(est, "nlmixr2Est")
  nlmixr2Est(.env)
}

##' @rdname nlmixr2
##' @export
nlmixr2.rxUi <- function(object, data, est = NULL, control = NULL, table = tableControl(), ...,
                         save = NULL, envir = parent.frame()) {
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .modName <- deparse(substitute(object))
  .uif <- object
  if (missing(data) && missing(est)) {
    return(.uif)
  } else if (missing(data)) {
    stop("need data", call.=FALSE)
  }
  .env <- new.env(parent=emptyenv())
  .env$ui <- .uif
  .env$data <- data
  .env$control <- control
  .env$table <- table
  class(.env) <- c(est, "nlmixr2Est")
  nlmixr2Est(.env)
}

#' @rdname nlmixr2
#' @export
nlmixr2.nlmixr2FitCore <- function(object, data, est = NULL, control = NULL, table = tableControl(), ...,
                                   save = NULL, envir = parent.frame()) {
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .modName <- deparse(substitute(object))
  .uif <- object
  if (missing(data)) {
    data <- object$origData
  }
  if (missing(est)) {
    est <- object$est
  }
  if (missing(control)) {
    control <- object$control
  }
  if (missing(table)) {
    table <- object$table
  }
  .env <- new.env(parent=emptyenv())
  .env$ui <- object$ui
  .env$data <- data
  .env$control <- control
  .env$table <- table
  class(.env) <- c(est, "nlmixr2Est")
  nlmixr2Est(.env)
}

#' @rdname nlmixr2
#' @export
nlmixr2.nlmixr2FitData <- nlmixr2.nlmixr2FitCore