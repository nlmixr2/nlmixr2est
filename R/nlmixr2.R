#' @importFrom stats predict logLik na.fail pchisq approxfun cov cov2cor dlnorm median na.omit qchisq
#' @importFrom n1qn1 n1qn1
#' @importFrom brew brew
#' @importFrom nlme nlme fixed.effects random.effects
#' @importFrom nlme groupedData
#' @importFrom nlme getData
#' @importFrom nlme pdDiag
#' @importFrom rxode2 rxode2
#' @importFrom graphics abline lines matplot plot points title
#' @importFrom stats as.formula nlminb optimHess rnorm terms predict anova optim sd var AIC BIC asOneSidedFormula coef end fitted resid setNames start simulate nobs qnorm quantile time
#' @importFrom utils assignInMyNamespace getFromNamespace head stack sessionInfo tail str getParseData .DollarNames
#' @importFrom parallel mclapply
#' @importFrom methods is
#' @importFrom Rcpp evalCpp
#' @importFrom lbfgsb3c lbfgsb3c
#' @importFrom ggplot2 ggplot aes geom_point facet_wrap geom_line geom_abline xlab geom_smooth aes_string
#' @importFrom rxode2 rxUiGet .malert .minfo .msuccess .mwarn
#' @useDynLib nlmixr2est, .registration=TRUE

rex::register_shortcuts("nlmixr2est")
## GGplot use and other issues...
utils::globalVariables(c("DV", "ID", "IPRED", "IRES", "PRED", "TIME", "grp", "initCondition", "values", "nlmixr2_pred", "iter", "val", "EVID"))

nlmixr2.logo <- "         _             _             \n        | | %9s (_) %s\n  _ __  | | _ __ ___   _ __  __ _ __\n | '_ \\ | || '_ ` _ \\ | |\\ \\/ /| '__|\n | | | || || | | | | || | >  < | |\n |_| |_||_||_| |_| |_||_|/_/\\_\\|_|\n"

#' Messages the nlmixr2 logo...
#'
#' @param str String to print
#' @param version Version information (by default use package version)
#' @return nothing; Called to display version information
#' @author Matthew L. Fidler
nlmixr2Logo <- function(str = "", version = sessionInfo()$otherPkgs$nlmixr2$Version) {
  message(sprintf(nlmixr2.logo, str, version))
}
#' Display nlmixr2's version
#'
#' @author Matthew L. Fidler
#' @return Nothing, called for its side effects
#' @export
nlmixr2Version <- function() {
  nlmixr2Logo()
}

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
#' @param data nlmixr data
#' @param est estimation method
#' @param control The estimation control object.  These are expected
#'   to be different for each type of estimation method
#' @param table The output table control object (like
#'   `tableControl()`)
#' @param ... Other parameters
#' @param save Boolean to save a nlmixr2 object in a rds file in the
#'   working directory.  If \code{NULL}, uses option "nlmixr2.save"
#' @param envir Environment where the nlmixr object/function is
#'   evaluated before running the estimation routine.
#' @return Either a nlmixr2 model or a nlmixr2 fit object
#' @author Matthew L. Fidler
#' @examples
#'
#' \donttest{
#'
#' one.cmt <- function() {
#'  ini({
#'    ## You may label each parameter with a comment
#'    tka <- 0.45 # Ka
#'    tcl <- log(c(0, 2.7, 100)) # Log Cl
#'    ## This works with interactive models
#'    ## You may also label the preceding line with label("label text")
#'    tv <- 3.45; label("log V")
#'    ## the label("Label name") works with all models
#'    eta.ka ~ 0.6
#'    eta.cl ~ 0.3
#'    eta.v ~ 0.1
#'    add.sd <- 0.7
#'    prop.sd <- 0.01
#'  })
#'  model({
#'    ka <- exp(tka + eta.ka)
#'    cl <- exp(tcl + eta.cl)
#'    v <- exp(tv + eta.v)
#'    linCmt() ~ add(add.sd) + prop(prop.sd)
#'  })
#' }
#'
#' fitF <- nlmixr(one.cmt, theo_sd, "focei")
#'
#' fitS <- nlmixr(one.cmt, theo_sd, "saem")
#'
#' }
#'
#' @export
nlmixr2 <- function(object, data, est = NULL, control = list(),
                    table = tableControl(), ..., save = NULL,
                    envir = parent.frame()) {
  rxode2::rxUnloadAll()
  assignInMyNamespace(".nlmixr2Time", proc.time())
  on.exit(.finalizeOverallTiming(), add=TRUE)
  nmSuppressMsg()
  rxode2::rxSuppressMsg()
  rxode2::rxSolveFree()
  rxode2::.setWarnIdSort(FALSE)
  on.exit(rxode2::.setWarnIdSort(TRUE), add=TRUE)
  force(est)
  ## verbose?
  ## https://tidymodels.github.io/model-implementation-principles/general-conventions.html
  UseMethod("nlmixr2")
}

#' @rdname nlmixr2
#' @export
nlmixr <- nlmixr2

#' @rdname nlmixr2
#' @export
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
  .env$missingTable <- missing(table)
  .env$missingControl <- missing(control)
  .env$missingEst <- missing(est)
  class(.env) <- c(est, "nlmixr2Est")
  nlmixr2Est0(.env)
}

#' @rdname nlmixr2
#' @export
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
  .env$missingTable <- missing(table)
  .env$missingControl <- missing(control)
  .env$missingEst <- missing(est)
  class(.env) <- c(est, "nlmixr2Est")
  nlmixr2Est0(.env)
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
  nlmixr2Est0(.env)
}

#' @rdname nlmixr2
#' @export
nlmixr2.nlmixr2FitData <- nlmixr2.nlmixr2FitCore
