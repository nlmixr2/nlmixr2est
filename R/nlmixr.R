

##' @importFrom stats predict logLik na.fail pchisq approxfun cov cov2cor dlnorm median na.omit qchisq
##' @importFrom n1qn1 n1qn1
##' @importFrom brew brew
##' @importFrom nlme nlme fixed.effects random.effects
##' @importFrom nlme groupedData
##' @importFrom nlme getData
##' @importFrom nlme pdDiag
##' @importFrom rxode2 rxode2
##' @importFrom graphics abline lines matplot plot points title
##' @importFrom stats as.formula nlminb optimHess rnorm terms predict anova optim sd var AIC BIC asOneSidedFormula coef end fitted resid setNames start simulate nobs qnorm quantile time
##' @importFrom utils assignInMyNamespace getFromNamespace head stack sessionInfo tail str getParseData
##' @importFrom parallel mclapply
##' @importFrom methods is
##' @importFrom Rcpp evalCpp
##' @importFrom lbfgsb3c lbfgsb3c
##' @importFrom ggplot2 ggplot aes geom_point facet_wrap geom_line geom_abline xlab geom_smooth aes_string
##' @importFrom rxode2 rxUiGet .malert .minfo .msuccess .mwarn
##' @useDynLib nlmixr2, .registration=TRUE


rex::register_shortcuts("nlmixr2")
## GGplot use and other issues...
utils::globalVariables(c("DV", "ID", "IPRED", "IRES", "PRED", "TIME", "grp", "initCondition", "values", "nlmixr2_pred", "iter", "val", "EVID"))

nlmixr2.logo <- "         _             _             \n        | | %9s (_) %s\n  _ __  | | _ __ ___   _ __  __ _ __\n | '_ \\ | || '_ ` _ \\ | |\\ \\/ /| '__|\n | | | || || | | | | || | >  < | |\n |_| |_||_||_| |_| |_||_|/_/\\_\\|_|\n"

##' Messages the nlmixr2 logo...
##'
##' @param str String to print
##' @param version Version information (by default use package version)
##' @return nothing; Called to display version information
##' @author Matthew L. Fidler
nlmixr2Logo <- function(str = "", version = sessionInfo()$otherPkgs$nlmixr2$Version) {
  message(sprintf(nlmixr2.logo, str, version))
}
##' Display nlmixr2's version
##'
##' @author Matthew L. Fidler
##' @return Nothing, called for its side effects
##' @export
nlmixr2Version <- function() {
  nlmixr2Logo()
}

.nlmixr2Time <- NULL

##' @rdname nlmixr2
##' @export
nlmixr2000.nlmixr2FitCore <- function(object, data, est = NULL, control = list(), table = tableControl(), ...,
                                 save = NULL, envir = parent.frame()) {
  .uif <- .getUif(object)
  if (missing(data)) {
    data <- getData(object)
  }
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .args$data <- data
  .args$est <- est
  .args <- c(list(uif = .uif), .args[-1])
  return(do.call(nlmixr2_fit, .args, envir = envir))
}

##' @rdname nlmixr2
##' @export
nlmixr20000.nlmixr2UI <- function(object, data, est = NULL, control = list(), ...,
                            save = NULL, envir = parent.frame()) {
  .args <- as.list(match.call(expand.dots = TRUE))[-1]
  .uif <- object
  if (missing(data) && missing(est)) {
    return(.uif)
  } else {
    .args <- c(list(uif = .uif), .args[-1])
    if (missing(data) && !is.null(.getPipedData())) {
      data <- .getPipedData()
      .args$data <- data
      .args$est <- est
    } else {
      .uif$nmodel$data.name <- deparse(substitute(data))
      .args$data <- data
      .args$est <- est
    }
    return(do.call(nlmixr2_fit, .args, envir = envir))
  }
}

##' Convert/Format the data appropriately for nlmixr2
##'
##' @param data is the name of the data to convert.  Can be a csv file
##'     as well.
##' @param model This is the rxode2 model to use to translate against
##'     when parsing the data.
##' @return Appropriately formatted data
##' @author Matthew L. Fidler
##' @keywords internal
##' @export
nlmixr2Data <- function(data, model = NULL) {
  UseMethod("nlmixr2Data")
}
##' @export
##' @rdname nlmixr2Data
nlmixr2Data.character <- function(data, model = NULL) {
  if (!file.exists(data)) {
    stop(sprintf("%s does not exist.", data))
  }
  if (regexpr(rex::rex(".csv", end), data) != -1) {
    return(nlmixr2Data.default(utils::read.csv(data, na.strings = c(".", "NA", "na", ""))))
  } else {
    stop(sprintf("Do not know how to read in %s", data))
  }
}
##' @export
##' @rdname nlmixr2Data
nlmixr2Data.default <- function(data, model = NULL) {
  if (!is.null(model)) {
    dat <- rxode2::etTrans(data, model, addCmt = TRUE, dropUnits = TRUE, allTimeVar = TRUE)
  } else {
    dat <- .as.data.frame(data)
  }
  return(dat)
}

nlmixr2_fit0 <- function(uif, data, est = NULL, control = list(), ...,
                        keep=NULL, drop=NULL,
                        sum.prod = FALSE, table = tableControl(),
                        envir = parent.frame()) {
  if (is.null(est)) {
    stop("Estimation type must be specified by est=''")
  }
  .clearPipedData()
  .tmp <- deparse(body(uif$theta.pars))[-1]
  .tmp <- .tmp[-length(.tmp)]
  .origData <- data
  .meta <- uif$meta
  .drop <- NULL
  .keep <- NULL
  args <- as.list(match.call(expand.dots = TRUE))[-1]
  if (exists("drop", envir = .meta)) {
    .drop <- .meta$drop
    checkmate::assertCharacter(.drop, min.len=1, min.chars=1, .var.name="drop")
  }
  if (exists("keep", envir = .meta)) {
    .keep <- .meta$keep
    checkmate::assertCharacter(.keep, min.len=1, min.chars=1, .var.name="keep")
  }
  if (!is.null(control$keep)) {
    if (inherits(.keep, "character")) {
      warning("control(keep=) overwrites keep in model")
    }
    .keep <- control$keep
    checkmate::assertCharacter(.keep, min.len=1, min.chars=1, .var.name="keep")
  }
  if (!is.null(control$drop)) {
    if (inherits(.drop, "character")) {
      warning("control(drop=) overwrites drop in model")
    }
    .drop <- control$drop
    checkmate::assertCharacter(.drop, min.len=1, min.chars=1, .var.name="drop")
  }
  if (!is.null(keep)) {
    if (inherits(.keep, "character")) {
      warning("keep= overwrites other ways of specifying keep")
    }
    .keep <- keep
    checkmate::assertCharacter(.keep, min.len=1, min.chars=1, .var.name="keep")
  }
  if (!is.null(drop)) {
    if (inherits(.drop, "character")) {
      warning("drop= overwrites other ways of specifying drop")
    }
    .drop <- drop
    checkmate::assertCharacter(.drop, min.len=1, min.chars=1, .var.name="drop")

  }
  .extra <- ""
  if (inherits(.keep, "character")) {
    .extra <- paste("nlmixr2Extra~", paste(.keep, collapse="+"), "\n")
  }
  data <- rxode2::etTrans(inData=data, obj=paste(paste(.tmp, collapse = "\n"), "\n", uif$rxode, "\n", .extra),
                         addCmt=TRUE, dropUnits=TRUE, allTimeVar=TRUE, keepDosingOnly=FALSE)

  .nTv <- attr(class(data), ".rxode2.lst")$nTv
  .lab <- attr(class(data), ".rxode2.lst")$idLvl
  .nid <- attr(class(data), ".rxode2.lst")$nid
  .modelId <-
    digest::digest(
      list(
        sessionInfo()$otherPkgs$nlmixr2$Version,
        uif, data, est, control, sum.prod, table, ...
      )
    )
  .missingEst <- is.null(est)
  if (.missingEst & exists("est", envir = .meta)) {
    est <- .meta$est
  }
  if (.missingEst & missing(control) & exists("control", envir = .meta)) {
    control <- .meta$control
    if (is(control, "foceiControl")) {
      est <- "focei"
      if (.missingEst && est != "focei") {
        warning(sprintf("Using focei instead of %s since focei controls were specified.", est))
      }
    } else if (is(control, "saemControl")) {
      est <- "saem"
      if (.missingEst && est != "saem") {
        warning(sprintf("Using saem instead of %s since saem controls were specified.", est))
      }
    }
  }
  if (missing(table) && exists("table", envir = .meta)) {
    table <- .meta$table
  }
  start.time <- Sys.time()
  if (!is(table, "tableControl")) {
    if (is(table, "list")) {
      table <- do.call(tableControl, table, envir = envir)
    } else {
      table <- tableControl()
    }
  }

  dat <- nlmixr2Data(data)
  nobs2 <- sum(dat$EVID == 0)
  up.covs <- toupper(uif$all.covs)
  up.names <- toupper(names(dat))
  for (i in seq_along(up.covs)) {
    w <- which(up.covs[i] == up.names)
    if (length(w) == 1) {
      names(dat)[w] <- uif$all.covs[i]
    }
  }
  ## backSort <- attr(dat, "backSort")
  ## backSort2 <- attr(dat, "backSort2")
  ## attr(dat, "backSort") <- NULL
  ## attr(dat, "backSort2") <- NULL
  if (!is.null(uif$nmodel$lin.solved)) {
    uif$env$infusion <- max(dat$EVID) > 10000
  }
  bad.focei <- "Problem calculating residuals, returning fit without residuals."
  calc.resid <- table$cwres
  if (is.null(calc.resid)) {
    if (est == "saem") {
      calc.resid <- table$saemCWRES
    } else if (est == "nlme") {
      calc.resid <- table$nlmeCWRES
    }
  }
  .cur <- environment()
  class(.cur) <- c(est, "nlmixr2Est000")
  .ret <- nlmixr2Est000(.cur, ...)
  return(.ret)
}

##' Fit a nlmixr2 model
##'
##' @param data Dataset to estimate.  Needs to be rxode2 compatible (see
##'   \url{https://nlmixr2development.github.io/rxode2/articles/rxode2-event-types.html}
##'   for detailed dataset requirements).
##' @param uif Parsed nlmixr2 model (by \code{nlmixr2(mod.fn)}).
##' @param est Estimation method
##' @param control Estimation control options.  They could be
##'   \code{\link[nlme]{nlmeControl}}, \code{\link{saemControl}} or
##'   \code{\link{foceiControl}}
##' @param ... Parameters passed to estimation method.

##' @param sum.prod Take the rxode2 model and use more precise products/sums.
##'   Increases solving accuracy and solving time.
##' @param table A list controlling the table options (i.e. CWRES, NPDE etc).
##'   See \code{\link{tableControl}}.
##' @param save This option determines if the fit will be saved to be reloaded
##'   if already run.  If NULL, get the option from
##'   \code{options("nlmixr2.save")};
##' @param envir Environment that nlmixr2 is evaluated in.
##' @inheritParams foceiFit
##' @return nlmixr2 fit object
##' @author Matthew L. Fidler
##' @export
nlmixr2_fit <- function(uif, data, est = NULL, control = list(), ...,
                       sum.prod = FALSE, table = tableControl(),
                       keep=NULL, drop=NULL,
                       save = NULL, envir = parent.frame()) {
  rxode2::.setWarnIdSort(FALSE)
  on.exit(rxode2::.setWarnIdSort(TRUE))
  if (is.null(save)) {
    save <- getOption("nlmixr2.save", FALSE)
  }
  if (save) {
    .modName <- ifelse(is.null(uif$model.name), "", paste0(uif$model.name, "-"))
    if (.modName == ".-") .modName <- ""
    .dataName <- ifelse(is.null(uif$data.name), "", paste0(uif$data.name, "-"))
    if (.dataName == ".-") .dataName <- ""
    .digest <-
      digest::digest(
        list(
          gsub("<-", "=", gsub(" +", "", uif$fun.txt)),
          .as.data.frame(uif$ini),
          data,
          est,
          control,
          sum.prod,
          table,
          ...,
          as.character(utils::packageVersion("nlmixr2")),
          as.character(utils::packageVersion("rxode2"))
        )
      )
    .saveFile <- file.path(
      getOption("nlmixr2.save.dir", getwd()),
      paste0("nlmixr2-", .modName, .dataName, est, "-", .digest, ".rds")
    )
    if (file.exists(.saveFile)) {
      message(sprintf("Loading model already run (%s)", .saveFile))
      .ret <- readRDS(.saveFile)
      if (!is.null(.ret$warnings)) {
        sapply(.ret$warnings, warning)
      }
      return(.ret)
    }
  }
  .ret <-
    .collectWarnings(
      nlmixr2_fit0(
        uif = uif, data = data, est = est, control = control, ...,
        sum.prod = sum.prod, table = table, envir = envir,
        keep=keep, drop=drop
      ),
      TRUE
    )
  .ws <- .ret[[2]]
  .ret <- .ret[[1]]
  if (inherits(.ret, "nlmixr2FitCore")) {
    .env <- .ret$env
    .env$warnings <- .ws
  }
  for (.i in seq_along(.ws)) {
    warning(.ws[.i])
  }
  if (inherits(.ret, "nlmixr2FitCore")) {
    if (save) {
      AIC(.ret) # Calculate SAEM AIC when saving...
      .env <- .ret$env
      .extra <- (proc.time() - .nlmixr2Time)["elapsed"] - sum(.env$time)
      .env$time <- .data.frame(.env$time, "other" = .extra, check.names = FALSE)
      saveRDS(.ret, file = .saveFile)
    } else {
      .env <- .ret$env
      .extra <- (proc.time() - .nlmixr2Time)["elapsed"] - sum(.env$time)
      .env$time <- .data.frame(.env$time, "other" = .extra, check.names = FALSE)
    }
  }
  return(.ret)
}

##' Control Options for SAEM
##'
##' @param seed Random Seed for SAEM step.  (Needs to be set for
##'     reproducibility.)  By default this is 99.
##'
##' @param nBurn Number of iterations in the Stochastic Approximation
##'     (SA), or burn-in step. This is equivalent to Monolix's \code{K_0} or \code{K_b}.
##'
##' @param nEm Number of iterations in the Expectation-Maximization
##'     (EM) Step. This is equivalent to Monolix's \code{K_1}.
##'
##' @param nmc Number of Markov Chains. By default this is 3.  When
##'     you increase the number of chains the numerical integration by
##'     MC method will be more accurate at the cost of more
##'     computation.  In Monolix this is equivalent to \code{L}
##'
##' @param nu This is a vector of 3 integers. They represent the
##'     numbers of transitions of the three different kernels used in
##'     the Hasting-Metropolis algorithm.  The default value is \code{c(2,2,2)},
##'     representing 40 for each transition initially (each value is
##'     multiplied by 20).
##'
##'     The first value represents the initial number of multi-variate
##'     Gibbs samples are taken from a normal distribution.
##'
##'     The second value represents the number of uni-variate, or multi-
##'     dimensional random walk Gibbs samples are taken.
##'
##'     The third value represents the number of bootstrap/reshuffling or
##'     uni-dimensional random samples are taken.
##'
##' @param print The number it iterations that are completed before
##'     anything is printed to the console.  By default, this is 1.
##'
##' @param covMethod  Method for calculating covariance.  In this
##'     discussion, R is the Hessian matrix of the objective
##'     function. The S matrix is the sum of each individual's
##'     gradient cross-product (evaluated at the individual empirical
##'     Bayes estimates).
##'
##'  "\code{linFim}" Use the Linearized Fisher Information Matrix to calculate the covariance.
##'
##'  "\code{fim}" Use the SAEM-calculated Fisher Information Matrix to calculate the covariance.
##'
##'  "\code{r,s}" Uses the sandwich matrix to calculate the covariance, that is: \eqn{R^-1 \times S \times R^-1}
##'
##'  "\code{r}" Uses the Hessian matrix to calculate the covariance as \eqn{2\times R^-1}
##'
##'  "\code{s}" Uses the crossproduct matrix to calculate the covariance as \eqn{4\times S^-1}
##'
##'  "" Does not calculate the covariance step.
##'
##' @param logLik boolean indicating that log-likelihood should be
##'     calculate by Gaussian quadrature.
##'
##' @param trace An integer indicating if you want to trace(1) the
##'     SAEM algorithm process.  Useful for debugging, but not for
##'     typical fitting.
##'
##' @param nnodes.gq number of nodes to use for the Gaussian
##'     quadrature when computing the likelihood with this method
##'     (defaults to 1, equivalent to the Laplaclian likelihood)
##'
##' @param nsd.gq span (in SD) over which to integrate when computing
##'     the likelihood by Gaussian quadrature. Defaults to 3 (eg 3
##'     times the SD)
##'
##' @param adjObf is a boolean to indicate if the objective function
##'     should be adjusted to be closer to NONMEM's default objective
##'     function.  By default this is \code{TRUE}
##'
##' @param tol This is the tolerance for the regression models used
##'   for complex residual errors (ie add+prop etc)
##'
##' @param itmax This is the maximum number of iterations for the
##'   regression models used for complex residual errors.  The number
##'   of iterations is itmax*number of parameters
##'
##' @param ... Other arguments to control SAEM.
##'
##' @inheritParams rxode2::rxSolve
##' @inheritParams foceiControl
##' @inheritParams configsaem
##' @inheritParams nlmixr2_fit
##' @inheritParams rxode2::rxSEinner
##' @inheritParams rxode2::rxGenSaem
##' @return List of options to be used in \code{\link{nlmixr2}} fit for
##'     SAEM.
##' @author Wenping Wang & Matthew L. Fidler
##' @export
saemControl <- function(seed = 99,
                        nBurn = 200, nEm = 300,
                        nmc = 3,
                        nu = c(2, 2, 2),
                        atol = 1e-06,
                        rtol = 1e-04,
                        method = "liblsoda",
                        transitAbs = FALSE,
                        print = 1,
                        trace = 0,
                        covMethod = c("linFim", "fim", "r,s", "r", "s", ""),
                        calcTables = TRUE,
                        logLik = FALSE,
                        nnodes.gq = 3,
                        nsd.gq = 1.6,
                        optExpression = FALSE,
                        maxsteps = 100000L,
                        adjObf = TRUE,
                        sum.prod = FALSE,
                        addProp = c("combined2", "combined1"),
                        singleOde = TRUE,
                        tol = 1e-6,
                        itmax = 30,
                        type = c("nelder-mead", "newuoa"),
                        powRange = 10,
                        lambdaRange = 3,
                        loadSymengine=FALSE,
                        odeRecalcFactor=10^(0.5),
                        maxOdeRecalc=5L,
                        ...) {
  type <- match.arg(type)
  .xtra <- list(...)
  .rm <- c()
  if (missing(transitAbs) && !is.null(.xtra$transit_abs)) {
    transitAbs <- .xtra$transit_abs
    .rm <- c(.rm, "transit_abs")
  }
  if (missing(nBurn) && !is.null(.xtra$n.burn)) {
    nBurn <- .xtra$n.burn
    .rm <- c(.rm, "n.burn")
  }
  if (missing(nEm) && !is.null(.xtra$n.em)) {
    nEm <- .xtra$n.em
    .rm <- c(.rm, "n.em")
  }
  if (inherits(addProp, "numeric")) {
    if (addProp == 1) {
      addProp <- "combined1"
    } else if (addProp == 2) {
      addProp <- "combined2"
    } else {
      stop("addProp must be 1, 2, \"combined1\" or \"combined2\"", call.=FALSE)
    }
  } else {
    addProp <- match.arg(addProp)
  }
  .ret <- list(
    mcmc = list(niter = c(nBurn, nEm), nmc = nmc, nu = nu),
    ODEopt = rxode2::rxControl(
      atol = atol, rtol = rtol, method = method,
      transitAbs = transitAbs, maxsteps = maxsteps, ...
    ),
    seed = seed,
    print = print,
    DEBUG = trace,
    optExpression = optExpression,
    sum.prod = sum.prod,
    nnodes.gq = nnodes.gq,
    nsd.gq = nsd.gq,
    adjObf = adjObf,
    addProp = addProp,
    singleOde = singleOde,
    itmax = itmax,
    tol = tol,
    type = type,
    powRange = powRange,
    lambdaRange = lambdaRange,
    loadSymengine=loadSymengine,
    odeRecalcFactor=odeRecalcFactor,
    maxOdeRecalc=maxOdeRecalc,
    ...
  )
  if (length(.rm) > 0) {
    .ret <- .ret[!(names(.ret) %in% .rm)]
  }
  .ret[["covMethod"]] <- match.arg(covMethod)
  .ret[["logLik"]] <- logLik
  .ret[["calcTables"]] <- calcTables
  class(.ret) <- "saemControl"
  .ret
}


##'@rdname  nlmixr2Est000
##'@export
nlmixr2Est000.saem <- function(env, ...) {
  with(env, {
    if (.nid <= 1) stop("SAEM is for mixed effects models, try 'focei', which downgrades to nonlinear regression")
    pt <- proc.time()
    uif$env <- new.env(parent = emptyenv())
    .tv <- NULL
    if (.nTv != 0) {
      .tv <- names(data)[-seq(1, 6)]
    }
    uif$env$.curTv <- .tv
    if (length(uif$noMuEtas) > 0) {
      stop(sprintf("Cannot run SAEM since some of the parameters are not mu-referenced (%s)", paste(uif$noMuEtas, collapse = ", ")))
    }
    default <- saemControl()
    .getOpt <- function(arg, envir = parent.frame(1)) {
      if (arg %in% names(args)) {
        assign(paste0(".", arg), args[[arg]], envir = envir)
      } else if (arg %in% names(control)) {
        assign(paste0(".", arg), control[[arg]], envir = envir)
      } else if (arg %in% names(default)) {
        assign(paste0(".", arg), default[[arg]], envir = envir)
      }
    }
    for (a in c(
      "mcmc", "ODEopt", "seed", "print",
      "DEBUG", "covMethod", "calcTables",
      "logLik", "nnodes.gq",
      "nsd.gq", "nsd.gq", "adjObf",
      "optExpression", "addProp",
      "singleOde", "type", "tol", "itmax",
      "lambdaRange", "powRange", "loadSymengine"
    )) {
      .getOpt(a)
    }
    uif$env$optExpression <- .optExpression
    uif$env$singleOde <- .singleOde
    .addCov <- .covMethod == "linFim"

    if (uif$saemErr != "") {
      stop(paste0("For SAEM:\n", uif$saemErr))
    }
    if (is.null(uif$nlme.fun.mu)) {
      stop("SAEM requires all ETAS to be mu-referenced")
    }
    .err <- uif$ini$err
    .low <- uif$ini$lower
    .low <- .low[!is.na(.low) & is.na(.err)]
    .up <- uif$ini$upper
    .up <- .up[!is.na(.up) & is.na(.err)]
    if (any(.low != -Inf) | any(.up != Inf)) {
      warning("Bounds are ignored in SAEM", call. = FALSE)
    }
    uif$env$mcmc <- .mcmc
    uif$env$ODEopt <- .ODEopt
    uif$env$sum.prod <- sum.prod
    uif$env$covMethod <- .covMethod
    .dist <- uif$saem.distribution
    model <- uif$saem.model
    inits <- uif$saem.init
    if (length(uif$saem.fixed) > 0) {
      nphi <- attr(model$saem_mod, "nrhs")
      m <- cumsum(!is.na(matrix(inits$theta, byrow = TRUE, ncol = nphi)))
      fixid <- match(uif$saem.fixed, t(matrix(m, ncol = nphi)))
      names(inits$theta) <- rep("", length(inits$theta))
      names(inits$theta)[fixid] <- "FIXED"
    }
    .cfg <- configsaem(
      model = model, data = dat, inits = inits,
      mcmc = .mcmc, ODEopt = .ODEopt, seed = .seed,
      distribution = .dist, DEBUG = .DEBUG,
      addProp = .addProp, tol = .tol, itmax = .itmax, type = .type,
      powRange = .powRange, lambdaRange = .lambdaRange
    )
    if (is(.print, "numeric")) {
      .cfg$print <- as.integer(.print)
    }
    .cfg$cres <- uif$saem.cres
    .cfg$yj <- uif$saem.yj
    .cfg$lres <- uif$saem.lambda
    .cfg$low <- uif$saem.low
    .cfg$hi <- uif$saem.hi
    .cfg$propT <- uif$saem.propT
    .fit <- model$saem_mod(.cfg)
    .ret <-
      try(
        as.focei.saemFit(
          .fit, uif, pt,
          data = dat, calcResid = calc.resid, obf = .logLik,
          nnodes.gq = .nnodes.gq, nsd.gq = .nsd.gq, adjObf = .adjObf,
          calcCov = .addCov, calcTables = .calcTables,
          keep=.keep, drop=.drop, IDlabel=.lab, table=table
        ),
        silent = TRUE
      )
    if (inherits(.ret, "try-error")) {
      warning("Error converting to nlmixr2 UI object, returning saem object")
      return(.fit)
    }
    if (inherits(.ret, "nlmixr2FitCore")) {
      .ret <- nlmixr2FitUpdateParams(.ret, origData = .origData)
    }
    if (inherits(.ret, "nlmixr2FitCore")) {
      .env <- .ret$env
      .env$adjObj <- .adjObf
      .env$nnodes.gq <- .nnodes.gq
      .env$nsd.gq <- .nsd.gq
      assign("startTime", start.time, .env)
      assign("est", est, .env)
      assign("stopTime", Sys.time(), .env)
      assign("origControl", control, .env)
      assign("modelId", .modelId, .env)
    }
    return(.ret)
  })
}


##' @rdname nlmixr2Est000
##'@export
nlmixr2Est000.nlme <- function(env, ...) {
  with(env, {
    if (.nid <= 1) stop("nlme is for mixed effects models, try 'dynmodel' (need more than 1 individual)")
    if (.nTv != 0) stop("nlme does not support time-varying covariates (yet)")
    data <- .as.data.frame(data)
    if (length(uif$predDf$cond) > 1) stop("nlmixr2 nlme does not support multiple endpoints.")
    pt <- proc.time()
    est.type <- est
    if (est == "nlme.free") {
      fun <- uif$nlme.fun
      specs <- uif$nlme.specs
    } else if (est == "nlme.mu") {
      fun <- uif$nlme.fun.mu
      specs <- uif$nlme.specs.mu
    } else if (est == "nlme.mu.cov") {
      fun <- uif$nlme.fun.mu.cov
      specs <- uif$nlme.specs.mu.cov
    } else {
      if (!is.null(uif$nlme.fun.mu.cov)) {
        est.type <- "nlme.mu.cov"
        fun <- uif$nlme.fun.mu.cov
        specs <- uif$nlme.specs.mu.cov
      } else if (!is.null(uif$nlme.fun.mu)) {
        est.type <- "nlme.mu"
        fun <- uif$nlme.fun.mu
        specs <- uif$nlme.specs.mu
      } else {
        est.type <- "nlme.free"
        fun <- uif$nlme.fun
        specs <- uif$nlme.fun.specs
      }
    }
    grp.fn <- uif$grp.fn
    dat$nlmixr2.grp <-
      factor(apply(dat, 1, function(x) {
        cur <- x
        names(cur) <- names(dat)
        with(as.list(cur), {
          return(grp.fn())
        })
      }))
    dat$nlmixr2.num <- seq_along(dat$nlmixr2.grp)
    .addProp <- "combined2"
    if (!is.null(control$addProp)) .addProp <- control$addProp
    if (!any(.addProp == c("combined2", "combined1"))) stop("addProp needs to either be 'combined1' and 'combined2'")
    uif$env$.addProp <- .addProp
    weight <- uif$nlme.var
    if (sum.prod) {
      rxode <- rxode2::rxSumProdModel(uif$rxode.pred)
    } else {
      rxode <- uif$rxode.pred
    }
    .atol <- 1e-8
    if (!is.null(control$atol)) .atol <- control$atol
    .rtol <- 1e-8
    if (!is.null(control$rtol)) .rtol <- control$rtol
    .maxsteps <- 5000
    if (!is.null(control$maxstepsOde)) .maxsteps <- control$maxstepsOde
    if (is(weight, "varConstProp")) {
      control$sigma <- 1
    }
    fit <- nlme_ode(dat,
                    model = rxode,
                    par_model = specs,
                    par_trans = fun,
                    response = "nlmixr2_pred",
                    weight = weight,
                    verbose = TRUE,
                    control = control,
                    atol = .atol,
                    rtol = .rtol,
                    maxsteps = .maxsteps,
                    ...
                    )
    class(fit) <- c(est.type, class(fit))
    .ret <- try({
      as.focei.nlmixr2Nlme(fit, uif, pt, data = dat, calcResid = calc.resid, nobs2 = nobs2,
                          keep=.keep, drop=.drop, IDlabel=.lab, table=table)
    })
    if (inherits(.ret, "try-error")) {
      warning("Error converting to nlmixr2 UI object, returning nlme object")
      return(fit)
    }
    if (inherits(.ret, "nlmixr2FitCore")) {
      .ret <- nlmixr2FitUpdateParams(.ret, origData = .origData)
    }
    if (inherits(.ret, "nlmixr2FitCore")) {
      .env <- .ret$env
      assign("startTime", start.time, .env)
      assign("est", est, .env)
      assign("stopTime", Sys.time(), .env)
      assign("origControl", control, .env)
      assign("modelId", .modelId, .env)
    }
    return(.ret)
  })
}

##'@rdname nlmixr2Est000
##'@export
nlmixr2Est000.nlme.mu <- nlmixr2Est000.nlme

##'@rdname nlmixr2Est000
##'@export
nlmixr2Est000.nlme.mu.cov <- nlmixr2Est000.nlme

##'@rdname nlmixr2Est000
##'@export
nlmixr2Est000.nlme.free <- nlmixr2Est000.nlme

##'@rdname nlmixr2Est000
##'@export
nlmixr2Est000.posthoc <- function(env, ...) {
  with(env, {
    if (class(control) != "foceiControl") control <- do.call(nlmixr2::foceiControl, control)
    if (any(est == c("foce", "fo"))) {
      control$interaction <- FALSE
    }
    env <- new.env(parent = emptyenv())
    env$table <- table
    env$IDlabel <- .lab
    env$uif <- uif
    if (any(est == c("fo", "foi"))) {
      control$maxInnerIterations <- 0
      control$fo <- TRUE
      control$boundTol <- 0
      env$skipTable <- TRUE
    }
    uif$env$singleOde <- control$singleOde
    if (control$singleOde) {
      .mod <- uif$focei.rx1
      .pars <- NULL
    } else {
      .mod <- uif$rxode.pred
      .pars <- uif$theta.pars
    }
    fit <- foceiFit(dat,
      inits = uif$focei.inits,
      PKpars = .pars,
      ## par_trans=fun,
      model = .mod,
      pred = function() {
        return(nlmixr2_pred)
      },
      err = uif$error,
      lower = uif$focei.lower,
      upper = uif$focei.upper,
      fixed = uif$focei.fixed,
      thetaNames = uif$focei.names,
      etaNames = uif$eta.names,
      control = control,
      env = env,
      keep=.keep,
      drop=.drop,
      ...
    )
    if (any(est == c("fo", "foi"))) {
      ## Add posthoc.
      .default <- foceiControl()
      control$maxInnerIterations <- .default$maxInnerIterations
      control$maxOuterIterations <- 0L
      control$covMethod <- 0L
      control$fo <- 0L
      .uif <- fit$uif
      .thetas <- fit$theta
      for (.n in names(.thetas)) {
        .uif$ini$est[.uif$ini$name == .n] <- .thetas[.n]
      }
      .omega <- fit$omega
      for (.i in seq_along(.uif$ini$neta1)) {
        if (!is.na(.uif$ini$neta1[.i])) {
          .uif$ini$est[.i] <- .omega[.uif$ini$neta1[.i], .uif$ini$neta2[.i]]
        }
      }
      .time <- fit$time
      .objDf <- fit$objDf
      .message <- fit$env$message
      .time <- fit$time
      env <- new.env(parent = emptyenv())
      env$table <- table
      env$IDlabel <- .lab
      for (.w in c("cov", "covR", "covS", "covMethod")) {
        if (exists(.w, fit$env)) {
          assign(.w, get(.w, envir = fit$env), envir = env)
        }
      }
      env$time2 <- time
      env$uif <- .uif
      env$method <- "FO"
      if (control$singleOde) {
        .mod <- uif$focei.rx1
        .pars <- NULL
      } else {
        .mod <- uif$rxode.pred
        .pars <- uif$theta.pars
      }
      fit0 <-
        try(
          foceiFit(
            dat,
            inits = .uif$focei.inits,
            PKpars = .pars,
            ## par_trans=fun,
            model = .mod,
            pred = function() {
              return(nlmixr2_pred)
            },
            err = .uif$error,
            lower = .uif$focei.lower,
            upper = .uif$focei.upper,
            fixed = .uif$focei.fixed,
            thetaNames = .uif$focei.names,
            etaNames = .uif$eta.names,
            control = control,
            env = env,
            keep=.keep,
            drop=.drop,
            ...
          ),
          silent = TRUE
        )
      if (inherits(fit0, "try-error")) {
      } else {
        fit <- fit0
        assign("message2", fit$env$message, env)
        assign("message", .message, env)
        .tmp1 <- env$objDf
        if (any(names(.objDf) == "Condition Number")) {
          .tmp1 <- .data.frame(.tmp1, "Condition Number" = NA, check.names = FALSE)
        }
        if (any(names(.tmp1) == "Condition Number")) {
          .objDf <- .data.frame(.objDf, "Condition Number" = NA, check.names = FALSE)
        }
        env$objDf <- rbind(.tmp1, .objDf)
        row.names(env$objDf) <- c(ifelse(est == "fo", "FOCE", "FOCEi"), "FO")
        .tmp1 <- env$time
        .tmp1$optimize <- .time$optimize
        .tmp1$covariance <- .time$covariance
        assign("time", .tmp1, envir = env)
        setOfv(env, "fo")
      }
    }
    fit <- nlmixr2FitUpdateParams(fit, origData = .origData)
    assign("start.time", start.time, env)
    assign("est", est, env)
    assign("stop.time", Sys.time(), env)
    assign("origControl", control, env)
    assign("modelId", .modelId, env)
    return(fit)
  })
}

##'@rdname nlmixr2Est000
##'@export
nlmixr2Est000.focei <- function(env, ...) {
  with(env, {
    if (class(control) != "foceiControl") control <- do.call(nlmixr2::foceiControl, control)
    if (any(est == c("foce", "fo"))) {
      control$interaction <- FALSE
    }
    env <- new.env(parent = emptyenv())
    env$table <- table
    env$IDlabel <- .lab
    env$uif <- uif
    if (any(est == c("fo", "foi"))) {
      control$maxInnerIterations <- 0
      control$fo <- TRUE
      control$boundTol <- 0
      env$skipTable <- TRUE
    }
    uif$env$singleOde <- control$singleOde
    if (control$singleOde) {
      .mod <- uif$focei.rx1
      .pars <- NULL
    } else {
      .mod <- uif$rxode.pred
      .pars <- uif$theta.pars
    }
    .FoceiInits <- uif$focei.inits
    if (.nid == 1) {
      if (length(.FoceiInits$OMGA) > 0) {
        stop("a mixed effect model requires more than one subject/id\na population estimate requires no etas", call.=FALSE)
      }
    }
    fit <- foceiFit(dat,
                    inits = .FoceiInits,
                    PKpars = .pars,
                    ## par_trans=fun,
                    model = .mod,
                    pred = function() {
                      return(nlmixr2_pred)
                    },
                    err = uif$error,
                    lower = uif$focei.lower,
                    upper = uif$focei.upper,
                    fixed = uif$focei.fixed,
                    thetaNames = uif$focei.names,
                    etaNames = uif$eta.names,
                    control = control,
                    env = env,
                    keep=.keep,
                    drop=.drop,
                    ...
                    )
    if (any(est == c("fo", "foi"))) {
      ## Add posthoc.
      .default <- foceiControl()
      control$maxInnerIterations <- .default$maxInnerIterations
      control$maxOuterIterations <- 0L
      control$covMethod <- 0L
      control$fo <- 0L
      .uif <- fit$uif
      .thetas <- fit$theta
      for (.n in names(.thetas)) {
        .uif$ini$est[.uif$ini$name == .n] <- .thetas[.n]
      }
      .omega <- fit$omega
      for (.i in seq_along(.uif$ini$neta1)) {
        if (!is.na(.uif$ini$neta1[.i])) {
          .uif$ini$est[.i] <- .omega[.uif$ini$neta1[.i], .uif$ini$neta2[.i]]
        }
      }
      .time <- fit$time
      .objDf <- fit$objDf
      .message <- fit$env$message
      .time <- fit$time
      env <- new.env(parent = emptyenv())
      env$table <- table
      env$IDlabel <- .lab
      for (.w in c("cov", "covR", "covS", "covMethod")) {
        if (exists(.w, fit$env)) {
          assign(.w, get(.w, envir = fit$env), envir = env)
        }
      }
      env$time2 <- time
      env$uif <- .uif
      env$method <- "FO"
      if (control$singleOde) {
        .mod <- uif$focei.rx1
        .pars <- NULL
      } else {
        .mod <- uif$rxode.pred
        .pars <- uif$theta.pars
      }
      fit0 <-
        try(
          foceiFit(
            dat,
            inits = .uif$focei.inits,
            PKpars = .pars,
            ## par_trans=fun,
            model = .mod,
            pred = function() {
              return(nlmixr2_pred)
            },
            err = .uif$error,
            lower = .uif$focei.lower,
            upper = .uif$focei.upper,
            fixed = .uif$focei.fixed,
            thetaNames = .uif$focei.names,
            etaNames = .uif$eta.names,
            control = control,
            env = env,
            keep=.keep,
            drop=.drop,
            ...
          ),
          silent = TRUE
        )
      if (inherits(fit0, "try-error")) {
      } else {
        fit <- fit0
        assign("message2", fit$env$message, env)
        assign("message", .message, env)
        .tmp1 <- env$objDf
        if (any(names(.objDf) == "Condition Number")) {
          .tmp1 <- .data.frame(.tmp1, "Condition Number" = NA, check.names = FALSE)
        }
        if (any(names(.tmp1) == "Condition Number")) {
          .objDf <- .data.frame(.objDf, "Condition Number" = NA, check.names = FALSE)
        }
        env$objDf <- rbind(.tmp1, .objDf)
        row.names(env$objDf) <- c(ifelse(est == "fo", "FOCE", "FOCEi"), "FO")
        .tmp1 <- env$time
        .tmp1$optimize <- .time$optimize
        .tmp1$covariance <- .time$covariance
        assign("time", .tmp1, envir = env)
        setOfv(env, "fo")
      }
    }
    fit <- nlmixr2FitUpdateParams(fit, origData = .origData)
    assign("start.time", start.time, env)
    assign("est", est, env)
    assign("stop.time", Sys.time(), env)
    assign("origControl", control, env)
    assign("modelId", .modelId, env)
    return(fit)
  })
}

##'@rdname nlmixr2Est000
##'@export
nlmixr2Est000.foce <- nlmixr2Est000.focei


##'@rdname nlmixr2Est000000
##'@export
nlmixr2Est000.fo <- nlmixr2Est000.focei

##'@rdname nlmixr2Est000
##'@export
nlmixr2Est000.foi <- nlmixr2Est000.focei

 'y'##'@rdname nlmixr2Est000
##'@export
nlmixr2Est000.posthoc <- function(env, ...){
  with(env, {
    if (.nid <= 1) stop("'posthoc' estimation is for mixed effects models, try 'dynmodel' (need more than 1 individual)")
    if (class(control) != "foceiControl") control <- do.call(nlmixr2::foceiControl, control)
    control$covMethod <- 0L
    control$maxOuterIterations <- 0L
    .env <- new.env(parent = emptyenv())
    .env$table <- table
    .env$IDlabel <- .lab
    .env$uif <- uif
    if (control$singleOde) {
      .mod <- uif$focei.rx1
      .pars <- NULL
    } else {
      .mod <- uif$rxode.pred
      .pars <- uif$theta.pars
    }
    fit <- foceiFit(dat,
                    inits = uif$focei.inits,
                    PKpars = .pars,
                    ## par_trans=fun,
                    model = .mod,
                    pred = function() {
                      return(nlmixr2_pred)
                    },
                    err = uif$error,
                    lower = uif$focei.lower,
                    upper = uif$focei.upper,
                    thetaNames = uif$focei.names,
                    etaNames = uif$eta.names,
                    control = control,
                    env = .env,
                    keep=.keep,
                    drop=.drop,
                    ...
                    )
    if (inherits(fit, "nlmixr2FitData")) {
      .cls <- class(fit)
      .env <- attr(.cls, ".foceiEnv")
      .cls[1] <- "nlmixr2Posthoc"
      class(fit) <- .cls
    }
    assign("uif", .syncUif(fit, fit$popDf, fit$omega), fit$env)
    fit <- nlmixr2FitUpdateParams(fit, origData = .origData)
    assign("origControl", control, fit$env)
    assign("modelId", .modelId, fit$env)
    return(fit)
  })
}

##'@rdname nlmixr2Est000
##'@export
nlmixr2Est000.nlmixr2Est000 <- function(env, ...){
  with(env, stop("unknown estimation method est=\"", est, "\""))
}
