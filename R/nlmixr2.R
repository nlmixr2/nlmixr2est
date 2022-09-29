#' @importFrom stats predict logLik na.fail pchisq approxfun cov cov2cor dlnorm median na.omit qchisq
#' @importFrom symengine S
#' @importFrom n1qn1 n1qn1
#' @importFrom nlme nlme fixed.effects random.effects
#' @importFrom nlme groupedData
#' @importFrom nlme getData
#' @importFrom nlme pdDiag
#' @importFrom rxode2 rxode2
#' @importFrom graphics abline lines matplot plot points title
#' @importFrom stats as.formula nlminb optimHess rnorm terms predict anova optim sd var AIC BIC asOneSidedFormula coef end fitted resid setNames start simulate nobs qnorm quantile time
#' @importFrom utils assignInMyNamespace getFromNamespace head stack sessionInfo tail str getParseData .DollarNames
#' @importFrom methods is
#' @importFrom Rcpp evalCpp
#' @importFrom lbfgsb3c lbfgsb3c
#' @importFrom rxode2 rxUiGet .malert .minfo .msuccess .mwarn
#' @import nlmixr2data
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
#' @export
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

.nlmixr2objectName <- NULL

#' Allows external methods (like those in nlmixr2) to assign object name
#'
#' @param x String or null for assigning a nlmixr object name
#' @return nothing called for side effects
#' @author Matthew L. Fidler
#' @export
#' @keywords internal
.nlmixr2objectNameAssign <- function(x) {
  assignInMyNamespace(".nlmixr2objectName",x)
  invisible()
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
#' @param est estimation method (all methods are shown by
#'   `nlmixr2AllEst()`). Methods can be added for other tools
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
  assignInMyNamespace(".finalUiCompressed", FALSE)
  .objectName <- try(as.character(substitute(object)), silent=TRUE)
  if (inherits(.objectName, "try-error")) .objectName <- "object"
  if (!identical(.objectName, "object")) {
    assignInMyNamespace(".nlmixr2objectName", .objectName)
  }
  on.exit(.finalizeOverallTiming(), add=TRUE)
  nmSuppressMsg()
  rxode2::rxSuppressMsg()
  rxode2::rxSolveFree()
  force(est)
  ## verbose?
  ## https://tidymodels.github.io/model-implementation-principles/general-conventions.html
  UseMethod("nlmixr2")
}

#' @rdname nlmixr2
#' @export
nlmixr <- nlmixr2

.nlmixr2pipeData <- NULL
.nlmixr2pipeControl <- NULL
.nlmixr2pipeTable <- NULL
.nlmixr2pipeEst <- NULL

.nlmixr2savePipe <- function(x) {
  assignInMyNamespace(".nlmixr2pipeData", x$origData)
  rxode2::rxSetCovariateNamesForPiping(names(x$origData))
  assignInMyNamespace(".nlmixr2pipeControl", x$control)
  assignInMyNamespace(".nlmixr2pipeTable", x$table)
  assignInMyNamespace(".nlmixr2pipeEst", x$est)
}

.nlmixr2clearPipe <- function(x) {
  assignInMyNamespace(".nlmixr2pipeData", NULL)
  assignInMyNamespace(".nlmixr2pipeControl", NULL)
  assignInMyNamespace(".nlmixr2pipeTable", NULL)
  assignInMyNamespace(".nlmixr2pipeEst", NULL)
  assignInMyNamespace(".finalUiCompressed", TRUE)
  rxode2::rxSetCovariateNamesForPiping(NULL)
}

#' @rdname nlmixr2
#' @export
nlmixr2.function <- function(object, data=NULL, est = NULL, control = NULL, table = tableControl(), ...,
                             save = NULL, envir = parent.frame()) {
  on.exit(.nlmixr2clearPipe())
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .uif <- rxode2::rxode2(object)
  .uif <- rxode2::rxUiDecompress(.uif)
  if (!is.null(.nlmixr2objectName)) {
    if (!identical(.nlmixr2objectName, "object")) {
      assign("modelName", .nlmixr2objectName, envir=.uif)
    }
  }
  .missingData <- FALSE
  if (is.null(data)) {
    .missingData <- TRUE
  }

  if (.missingData && missing(est)) {
    return(.uif)
  }
  .env <- new.env(parent=emptyenv())
  .env$ui <- .nlmixrPreprocessUi(.uif)
  if (is.null(data) && !is.null(.nlmixr2pipeData)) {
    .env$data <- .nlmixr2pipeData
    .minfo("use {.code data} from pipeline")
  } else if (.missingData) {
    stop("need data", call.=FALSE)
  } else {
    .env$data <- data
  }
  if (is.null(control) && !is.null(.nlmixr2pipeControl)) {
    .minfo("use {.code control} from pipeline")
    .env$control <- getValidNlmixrControl(.nlmixr2pipeControl, est)
  } else {
    .env$control <- getValidNlmixrControl(control, est)
  }
  if (is.null(table) && !is.null(.nlmixr2pipeTable)) {
    .env$table <- getValidNlmixrControl(.nlmixr2pipeTable, "tableControl")
    .minfo("use {.code table} from pipeline")
  } else {
    .env$table <- getValidNlmixrControl(table, "tableControl")
  }
  if (is.null(est) && !is.null(.nlmixr2pipeEst)) {
    .minfo("use {.code est} from pipeline")
    est <- .nlmixr2pipeEst
  }
  .env$missingTable <- missing(table)
  .env$missingControl <- missing(control)
  .env$missingEst <- missing(est)
  class(.env) <- c(est, "nlmixr2Est")
  nlmixr2Est0(.env)
}

#' @rdname nlmixr2
#' @export
nlmixr2.rxUi <- function(object, data=NULL, est = NULL, control = NULL, table = tableControl(), ...,
                         save = NULL, envir = parent.frame()) {
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .modelName <- try(as.character(substitute(object)), silent=TRUE)
  if (inherits(.modelName, "try-error")) .modelName <- NULL
  .uif <- object
  .uif <- rxode2::rxUiDecompress(.uif)
  if (is.null(.uif$modelName)) assign("modelName", .modelName, envir=.uif)
  if (is.null(data) && missing(est)) {
    return(.uif)
  }
  .env <- new.env(parent=emptyenv())
  .env$ui <- .nlmixrPreprocessUi(.uif)
  .missingData <- FALSE
  if (is.null(data)) {
    data <- NULL
  }
  if (is.null(data) && !is.null(.nlmixr2pipeData)) {
    .env$data <- .nlmixr2pipeData
    .minfo("use {.code data} from pipeline")
  } else if (.missingData) {
    stop("need data", call.=FALSE)
  } else {
    .env$data <- data
  }
  if (is.null(control) && !is.null(.nlmixr2pipeControl)) {
    .env$control <- getValidNlmixrControl(.nlmixr2pipeControl, est)
    .minfo("use {.code control} from pipeline")
  } else {
    .env$control <- getValidNlmixrControl(control, est)
  }
  if (is.null(table) && !is.null(.nlmixr2pipeTable)) {
    .env$table <- getValidNlmixrControl(.nlmixr2pipeTable, "tableControl")
    .minfo("use {.code table} from pipeline")
  } else {
    .env$table <- getValidNlmixrControl(table, "tableControl")
  }
  if (is.null(est) && !is.null(.nlmixr2pipeEst)) {
    est <- .nlmixr2pipeEst
    .minfo("use {.code est} from pipeline")
  }
  .env$missingTable <- missing(table)
  .env$missingControl <- missing(control)
  .env$missingEst <- missing(est)
  class(.env) <- c(est, "nlmixr2Est")
  nlmixr2Est0(.env)
}

.nlmixr2SimInfo <- NULL
#' @rdname nlmixr2
#' @export
nlmixr2.nlmixr2FitCore <- function(object, data=NULL, est = NULL, control = NULL, table = tableControl(), ...,
                                   save = NULL, envir = parent.frame()) {
  on.exit({
    .nlmixr2clearPipe()
    assignInMyNamespace(".nlmixr2SimInfo", NULL)
  })
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .modName <- deparse(substitute(object))
  assignInMyNamespace(".nlmixr2SimInfo", .simInfo(object))
  if (is.null(data) && !is.null(.nlmixr2pipeData)) {
    data <- .nlmixr2pipeData
    .minfo("use {.code data} from pipeline")
  }  else if (missing(data)) {
    data <- object$origData
    .minfo("use {.code data} from prior fit")
  }
  if (is.null(est) && !is.null(.nlmixr2pipeEst)) {
    est <- .nlmixr2pipeEst
    .minfo("use {.code est} from pipeline")
  } else if (missing(est)) {
    .minfo("use {.code est} from prior fit")
    est <- object$est
  }
  if (is.null(control) && !is.null(.nlmixr2pipeControl)) {
    .minfo("use {.code control} from pipeline")
    control <- getValidNlmixrControl(.nlmixr2pipeControl, est)
  } else if (is.null(control)) {
    .minfo("use {.code control} from prior fit")
    control <- getValidNlmixrControl(object$control, est)
  } else {
    control <- getValidNlmixrControl(control, est)
  }
  if (is.null(table) && !is.null(.nlmixr2pipeTable)) {
    table <- getValidNlmixrControl(.nlmixr2pipeTable, "tableControl")
    .minfo("use {.code table} from pipeline")
  } else if (is.null(table)) {
    .minfo("use {.code table} from prior fit")
    table <- getValidNlmixrControl(object$table, "tableControl")
  } else {
    table <- getValidNlmixrControl(table, "tableControl")
  }
  .env <- new.env(parent=emptyenv())
  .ui <- rxode2::rxUiDecompress(object$ui)
  .env$ui <- .nlmixrPreprocessUi(.ui)
  .env$data <- data
  .env$control <- control
  .env$table <- table
  class(.env) <- c(est, "nlmixr2Est")
  nlmixr2Est0(.env)
}

#' @rdname nlmixr2
#' @export
nlmixr2.nlmixr2FitData <- nlmixr2.nlmixr2FitCore
