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
#' @importFrom utils getFromNamespace head stack sessionInfo tail str getParseData .DollarNames
#' @importFrom methods is
#' @importFrom Rcpp evalCpp
#' @importFrom lbfgsb3c lbfgsb3c
#' @importFrom rxode2 rxUiGet .malert .minfo .msuccess .mwarn
#' @eval .nlmixr2estbuild()
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

#' Allows external methods (like those in nlmixr2) to assign object name
#'
#' @param x String or null for assigning a nlmixr object name
#' @return nothing called for side effects
#' @author Matthew L. Fidler
#' @export
#' @keywords internal
.nlmixr2objectNameAssign <- function(x) {
  nlmixr2global$nlmixr2objectName <- x
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
#'
#' @includeRmd man/rmdhunks/nlmixr2Keywords.Rmd
#'
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
#' # fitF <- nlmixr(one.cmt, theo_sd, "focei")
#'
#' fitS <- nlmixr(one.cmt, theo_sd, "saem")
#'
#' }
#'
#' @export
nlmixr2 <- function(object, data, est = NULL, control = list(),
                    table = tableControl(), ..., save = NULL,
                    envir = parent.frame()) {
  ## rxode2::rxUnloadAll() # don't unload everything anymore
  .nlmixr2globalReset()
  nlmixr2global$nlmixr2Time <- proc.time()
  nlmixr2global$finalUiCompressed <- FALSE
  nlmixr2global$nlmixrEvalEnv$envir <- envir
  if (inherits(object, "nlmixr2FitCore") &&
        !is.null(object$eta)) {
    nlmixr2global$etaMat <- object
  } else {
    nlmixr2global$etaMat <- NULL
  }
  .objectName <- try(as.character(substitute(object)), silent=TRUE)
  if (inherits(.objectName, "try-error")) .objectName <- "object"
  if (!identical(.objectName, "object")) {
    nlmixr2global$nlmixr2objectName <- .objectName
  }
  on.exit(.finalizeOverallTiming(), add=TRUE)
  nmSuppressMsg()
  rxode2::rxSuppressMsg()
  rxode2::rxSolveFree() # rxSolveFree unlocks evaluation environment
  # Add UDF environment for querying nlmixr2/rxode2 r-based user defined functions
  rxode2::.udfEnvSet(nlmixr2global$nlmixrEvalEnv$envir)
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
  nlmixr2global$nlmixr2pipeData <- x$origData
  rxode2::rxSetCovariateNamesForPiping(names(x$origData))
  nlmixr2global$nlmixr2pipeControl <- x$control
  nlmixr2global$nlmixr2pipeTable <- x$table
  nlmixr2global$nlmixr2pipeEst <- x$est
}

.nlmixr2clearPipe <- function(x) {
  nlmixr2global$nlmixr2pipeData <- NULL
  nlmixr2global$nlmixr2pipeControl <- NULL
  nlmixr2global$nlmixr2pipeTable <- NULL
  nlmixr2global$nlmixr2pipeEst <- NULL
  nlmixr2global$finalUiCompressed <- TRUE
  rxode2::rxSetCovariateNamesForPiping(NULL)
}
#' Infer missing estimation routine
#'
#' @param env prep environment for nlmixr2
#' @param est estimation routine, could actually be a control
#' @return actual estimation routine (could be inferred)
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2inferEst <- function(env, est) {
  if (!env$missingEst) {
    .cls <- class(est)
    if (length(.cls) == 1L && grepl("^.*?Control$", .cls)) {
      .est <- sub("^(.*?)Control$", "\\1", .cls)
      if (env$missingControl) {
        env$control <- getValidNlmixrControl(est, .est)
        est <- .est
        .minfo(paste0("infer estimation {.code ", est, "} from control"))
        env$missingControl <- FALSE
      }
    }
  }
  est
}

#' @rdname nlmixr2
#' @export
nlmixr2.function <- function(object, data=NULL, est = NULL, control = NULL, table = tableControl(), ...,
                             save = NULL, envir = parent.frame()) {
  on.exit(.nlmixr2clearPipe())
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .uif <- rxode2::rxode2(object)
  .uif <- rxode2::rxUiDecompress(.uif)
  if (!is.null(nlmixr2global$nlmixr2objectName)) {
    if (!identical(nlmixr2global$nlmixr2objectName, "object")) {
      assign("modelName", nlmixr2global$nlmixr2objectName, envir=.uif)
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
  .env$ui <- .uif
  if (is.null(data) && !is.null(.nlmixr2pipeData)) {
    .env$data <- .nlmixr2pipeData
    .minfo("use {.code data} from pipeline")
  } else if (.missingData) {
    stop("need data", call.=FALSE)
  } else {
    .env$data <- data
  }
  .env$missingTable <- missing(table)
  .env$missingControl <- missing(control)
  .env$missingEst <- missing(est)
  .env$control <- control
  est <- .nlmixr2inferEst(.env, est)
  control <- .env$control
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
  est <- .preProcessHooksRun(.env, est)
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
  .env$ui <- .uif
  .missingData <- FALSE
  if (is.null(data)) {
    data <- NULL
  }
  if (is.null(data) && !is.null(.nlmixr2pipeData)) {
    .env$data <- .nlmixr2pipeData
    .minfo("use {.code data} from pipeline")
  } else if (!is.null(.uif$nonmemData)) {
    .env$data <- .uif$nonmemData
    .minfo("use {.code data} from $nonmemData")
  } else if (.missingData) {
    stop("need data", call.=FALSE)
  } else {
    .env$data <- data
  }
  .env$missingTable <- missing(table)
  .env$missingControl <- missing(control)
  .env$missingEst <- missing(est)
  .env$control <- control
  est <- .nlmixr2inferEst(.env, est)
  control <- .env$control
  if (is.null(control) && !is.null(.nlmixr2pipeControl)) {
    .env$control <- getValidNlmixrControl(.nlmixr2pipeControl, est)
    .minfo("use {.code control} from pipeline")
  } else {
    .ctl <- getValidNlmixrControl(control, est)
    # The isTRUE() in the if below allows for the fact that
    # babelmixr2::pkncaControl() does not have a genRxControl element, so it is
    # NULL.
    if (isTRUE(.ctl$genRxControl)) {
      if (is.numeric(.uif$atol)) {
        .minfo("use rxControl(atol=) from $atol")
        .ctl$rxControl$atol <- .uif$atol
      }
      if (is.numeric(.uif$rtol)) {
        .minfo("use rxControl(rtol=) from $rtol")
        .ctl$rxControl$atol <- .uif$rtol
      }
      if (is.numeric(.uif$ssAtol)) {
        .minfo("use rxControl(ssAtol=) from $ssAtol")
        .ctl$rxControl$ssAtol <- .uif$ssAtol
      }
      if (is.numeric(.uif$ssRtol)) {
        .minfo("use rxControl(ssRtol=) from $ssRtol")
        .ctl$rxControl$ssRtol <- .uif$ssRtol
      }
      if (inherits(.uif, "nonmem2rx")) {
        .minfo("use rxControl(covsInterpolation=\"nocb\") since this model was import from NONMEM")
        .ctl$rxControl$covsInterpolation <- 2L
      }
    }
    .env$control <- .ctl
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
  est <- .preProcessHooksRun(.env, est)
  class(.env) <- c(est, "nlmixr2Est")
  nlmixr2Est0(.env)
}

#' @rdname nlmixr2
#' @export
nlmixr2.nlmixr2FitCore <- function(object, data=NULL, est = NULL, control = NULL, table = tableControl(), ...,
                                   save = NULL, envir = parent.frame()) {
  on.exit({
    .nlmixr2clearPipe()
    nlmixr2global$nlmixr2SimInfo <- NULL
  })
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .modName <- deparse(substitute(object))
  nlmixr2global$nlmixr2SimInfo <- .simInfo(object)
  .cls <- class(est)
  if (length(.cls) == 1L && grepl("^.*?Control$", .cls)) {
    .est <- sub("^(.*?)Control$", "\\1", .cls)
    if (is.null(control)) {
      control <- getValidNlmixrControl(est, .est)
      est <- .est
      .minfo(paste0("infer estimation {.code ", est, "} from control"))
    }
  }
  if (is.character(data) && length(data) == 1 &&
        data %fin% nlmixr2AllEst() &&
        is.null(est)) {
    est <- data
    data <- NULL
  } else {
    .cls <- class(data)
    if (length(.cls) == 1L && grepl("^.*?Control$", .cls)) {
      .est <- sub("^(.*?)Control$", "\\1", .cls)
      if (is.null(control)) {
        control <- getValidNlmixrControl(data, .est)
        est <- .est
        .minfo(paste0("infer estimation {.code ", est, "} from control"))
        data <- NULL
      }
    }
  }
  if (is.null(data) && !is.null(.nlmixr2pipeData)) {
    data <- .nlmixr2pipeData
    .minfo("use {.code data} from pipeline")
  } else if (missing(data)) {
    data <- object$origData
   .minfo("use {.code data} from prior/supplied fit")
  }
  if (!inherits(data, "data.frame")) {
    data <- object$origData
    .minfo("use {.code data} from prior/supplied fit")
  }
  if (is.null(est) && !is.null(.nlmixr2pipeEst)) {
    est <- .nlmixr2pipeEst
    .minfo("use {.code est} from pipeline")
  } else if (missing(est)) {
    .minfo("use {.code est} from prior/supplied fit")
    est <- object$est
  }
  .env <- new.env(parent=emptyenv())
  .env$control <- control
  if (!is.null(est)) {
    .env$missingEst <- FALSE
  } else {
    .env$missingEst <- missing(est)
  }
  est <- .nlmixr2inferEst(.env, est)
  control <- .env$control
  if (is.null(control) && !is.null(.nlmixr2pipeControl)) {
    .minfo("use {.code control} from pipeline")
    control <- getValidNlmixrControl(.nlmixr2pipeControl, est)
  } else if (is.null(control)) {
    .minfo("use/adapt {.code control} from prior/supplied fit")
    control <- getValidNlmixrControl(object$control, est)
  } else {
    control <- getValidNlmixrControl(control, est)
  }
  if (is.null(table) && !is.null(.nlmixr2pipeTable)) {
    table <- getValidNlmixrControl(.nlmixr2pipeTable, "tableControl")
    .minfo("use {.code table} from pipeline")
  } else if (is.null(table)) {
    .minfo("use {.code table} from prior/supplied fit")
    table <- getValidNlmixrControl(object$table, "tableControl")
  } else {
    table <- getValidNlmixrControl(table, "tableControl")
  }
  .ui <- rxode2::rxUiDecompress(object$ui)
  .env$ui <- .ui
  .env$data <- data
  .env$control <- control
  .env$table <- table
  .env$control <- control
  .env$missingEst <- missing(est)
  est <- .nlmixr2inferEst(.env, est)
  control <- .env$control
  est <- .preProcessHooksRun(.env, est)
  class(.env) <- c(est, "nlmixr2Est")
  nlmixr2Est0(.env)
}

#' @rdname nlmixr2
#' @export
nlmixr2.nlmixr2FitData <- nlmixr2.nlmixr2FitCore
