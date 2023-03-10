#' Get an item from a nlmixr core object
#'
#' @param x A specialized list with:
#' - First argument is a nlmixrFitCore environment
#' - Second argument is if the exact argument is requested
#' - The class would be the requested argument name followed by the class "nmObjGet"
#' @param ... Other arguments
#' @return Value of argument or NULL
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
nmObjGet <- function(x, ...) {
  if (!inherits(x, "nmObjGet")) {
    stop("'", as.character(substitute(x)), "' is wrong type for 'nmObjGet'", call.=FALSE)
  }
  .arg <- class(x)[1]
  if (any(.arg == c(
    "logLik", "value", "obf", "ofv",
    "objf", "OBJF", "objective", "AIC",
    "BIC"))) {
    .nmObjEnsureObjective(x[[1]])
  }
  UseMethod("nmObjGet")
}

#' @export
nmObjGet.iniUi <- function(x, ...) {
  .env <- x[[1]]
  .ui <- .cloneEnv(rxode2::rxUiDecompress(get("ui", .env)))
  .iniDf <- get("iniDf0", envir=.env)
  if (is.null(.iniDf)) return(NULL)
  assign("iniDf", .iniDf, envir=.ui)
  rxode2::rxUiCompress(.ui)
}
attr(nmObjGet.iniUi, "desc") <- "The initial ui used to run the model"


.finalUiCompressed <- TRUE

#' Set if the nlmixr2 object will return a compressed ui
#'
#' 
#' @param type is a boolean indicating if the compressed ui will be
#'   returned (`TRUE`) or not be returned (`FALSE`)
#' @return invisible logical type
#' @author Matthew L. Fidler
#' @export
#' @examples
#' 
#' nmObjUiSetCompressed(FALSE) # now the $ui will return an environment
#' nmObjUiSetCompressed(TRUE) # now the $ui will return a compressed value
#' 
nmObjUiSetCompressed <- function(type) {
  checkmate::assertLogical(type,len=1, any.missing=FALSE)
  assignInMyNamespace(".finalUiCompressed", type)
  invisible(type)
}

#' @export
nmObjGet.finalUi <- function(x, ...) {
  .env <- x[[1]]
  .ui <- .cloneEnv(rxode2::rxUiDecompress(get("ui", .env)))
  if (.finalUiCompressed) {
    rxode2::rxUiCompress(.ui)
  } else {
    rxode2::rxUiDecompress(.ui)
  }
}
attr(nmObjGet.finalUi, "desc") <- "The final ui used to run the model"

#' @export
nmObjGet.finalUiEnv <- function(x, ...) {
  .env <- x[[1]]
  .cloneEnv(rxode2::rxUiDecompress(get("ui", .env)))
}

#' @export
nmObjGet.ui <- nmObjGet.finalUi

#' Get an item from a nlmixr2FitData object
#'
#' @param x A specialized list with:
#' - First argument is a nlmixr2FitData object
#' - Second argument is if the exact argument is requested
#' - The class would be the requested argument name followed by the class "nmObjGet"
#' @param ... Other arguments
#' @return Value of argument or NULL
#' @author Matthew L. Fidler
#' @keywords internal
#' @export
nmObjGetData <- function(x, ...) {
  if (!inherits(x, "nmObjGetData")) {
    stop("'x' is wrong type for 'nmObjGetData'", call.=FALSE)
  }
  UseMethod("nmObjGetData")
}

#' @export
nmObjGetData.dataLloq <- function(x, ...) {
  .fit <- x[[1]]
  .df <- as.data.frame(.fit)
  if (!any(names(.df) == "CENS")) return(NULL)
  if (!any(names(.df) == "upperLim")) return(NULL)
  .w <- which(.df$CENS == 1)
  if (length(.w) == 0) return(NULL)
  mean(.df$upperLim[.w])
}

#' @export
nmObjGetData.dataUloq <- function(x, ...) {
  .fit <- x[[1]]
  .df <- as.data.frame(.fit) 
  if (!any(names(.df) == "CENS")) return(NULL)
  if (!any(names(.df) == "lowerLim")) return(NULL)
  .w <- which(.df$CENS == -1)
  if (length(.w) == 0) return(NULL)
  mean(.df$lowerLim[.w])
}

#' @export
nmObjGet.dataNormInfo <- function(x, ...) {
  .fit <- x[[1]]
  .ui <- .fit$ui
  .datSav <- .fit$dataSav
  .predDf <-.ui$predDf
  if (all(.predDf$dist %in% c("norm", "dnorm","t", "cauchy"))) {
    return(list(filter=rep(TRUE, length(.datSav[,1])),
                 nnorm=length(.datSav[,1]),
                 nlik=0,
                 nother=0,
                 nlmixrRowNums=.datSav$nlmixrRowNums))
  }
  .ret <- .Call(`_nlmixr2est_filterNormalLikeAndDoses`,
                .datSav$CMT, .predDf$distribution, .predDf$cmt)
  .ret$nlmixrRowNums <- .datSav[.ret$filter, "nlmixrRowNums"]
  .ret
}

#' @export
nmObjGet.warnings <-function(x, ...) {
  get("runInfo", x[[1]])
}

##' @export
nmObjGetData.default <- function(x, ...) {
  NULL
}

#' @rdname nmObjGet
#' @export
nmObjGet.default <- function(x, ...) {
  .arg <- class(x)[1]
  .env <- x[[1]]
  if (exists(.arg, envir = .env)) {
    .ret <- get(.arg, envir = .env)
    if (inherits(.ret, "raw")) .ret <- qs::qdeserialize(.ret)
    return(.ret)
  }
  # Now get the ui, install the control object temporarily and use `rxUiGet`
  .ui <- get("ui", envir=.env)
  .ui <- rxode2::rxUiDecompress(.ui)
  on.exit({
    assign("ui", rxode2::rxUiCompress(.ui), envir=.env)
  })
  .ctl <- nmObjGetControl(.createEstObject(x[[1]]), ...)
  if (!is.null(.ctl)) {
    assign("control", .ctl, envir=.ui)
  }
  on.exit(suppressWarnings(try(rm(list="control", envir=.ui), silent=TRUE)))
  .lst <- list(.ui, x[[2]])
  class(.lst) <- c(.arg, "rxUiGet")
  .ret <- rxUiGet(.lst)
  .ret
}

#'@rdname nmObjGet
#' @export
nmObjGet.modelName <- function(x, ...) {
  .obj <- x[[1]]
  .ui <- get("ui", .obj)
  .ui$modelName
}
attr(nmObjGet.modelName, "desc") <- "name of the model used for nlmixr2 model fit"

#' @rdname nmObjGet
#' @export
nmObjGet.cor <- function(x, ...) {
  .obj <- x[[1]]
  .cov <- .obj$cov
  .sd2 <- sqrt(diag(.cov))
  .cor <- stats::cov2cor(.cov)
  dimnames(.cor) <- dimnames(.cov)
  diag(.cor) <- .sd2
  .cor
}
attr(nmObjGet.cor, "desc") <- "correlation matrix of theta, calculated from covariance of theta"

#' @rdname nmObjGet
#' @export
nmObjGet.omegaR <- function(x, ...) {
  .obj <- x[[1]]
  .cov <- .obj$omega
  .sd2 <- sqrt(diag(.cov))
  .cor <- stats::cov2cor(.cov)
  dimnames(.cor) <- dimnames(.cov)
  diag(.cor) <- .sd2
  .cor
}
attr(nmObjGet.omegaR, "desc") <- "correlation matrix of omega"

#' @rdname nmObjGet
#' @export
nmObjGet.phiR <- function(x, ...) {
  .obj <- x[[1]]
  .phi <- .obj$phiC
  if (is.null(.phi)) {
    if (any(names(x[[1]]) !="CWRES")) warning("this requires 'CWRES' in fit (use `addCwres()`)", call.=FALSE)
    return(NULL)
  }
  .ret <- lapply(seq_along(.phi), function(i) {
    .cov <- .phi[[i]]
    .d <- diag(.cov)
    if (any(.d == 0) || any(!is.finite(.d))) {
      .d1 <- length(.d)
      return(matrix(rep(NA, .d1*.d1), .d1, .d1))
    }
    .sd2 <- sqrt(diag(.cov))
    .cor <- stats::cov2cor(.cov)
    dimnames(.cor) <- dimnames(.cov)
    diag(.cor) <- .sd2
    .cor
  })
  names(.ret) <- names(.phi)
  .ret
}
attr(nmObjGet.phiR, "desc") <- "correlation matrix of each individual's eta (if present)"

#' @rdname nmObjGet
#' @export
nmObjGet.phiSE <- function(x, ...) {
  .obj <- x[[1]]
  .phi <- .obj$phiC
  if (is.null(.phi)) {
    if (any(names(x[[1]]) !="CWRES")) warning("this requires 'CWRES' in fit (use `addCwres()`)", call.=FALSE)
    return(NULL)
  }
  .d1 <- dim(.phi[[1]])[1]
  .ret <- vapply(seq_along(.phi), function(i) {
    .cov <- .phi[[i]]
    suppressWarnings(sqrt(diag(.cov)))
  }, double(.d1), USE.NAMES=FALSE)
  dim(.ret) <- c(.d1, length(.phi))
  dimnames(.ret) <- list(colnames(.phi[[1]]), names(.phi))
  names(.ret) <- paste0("se(", names(.ret), ")")
  .ret <- as.data.frame(t(.ret))
  .id <- seq_along(.phi)
  if (!is.null(names(.phi))) {
    .id <- names(.phi)
  }
  cbind(data.frame(ID=.id), .ret)
}
attr(nmObjGet.phiSE, "desc") <- "standard error of each individual's eta (if present)"

#' @rdname nmObjGet
#' @export
nmObjGet.phiRSE <- function(x, ...) {
  .obj <- x[[1]]
  .phi <- .obj$phiC
  .eta <- .obj$eta[,-1, drop=FALSE]
  if (is.null(.phi)){
    if (any(names(x[[1]]) !="CWRES")) warning("this requires 'CWRES' in fit (use `addCwres()`)", call.=FALSE)
    return(NULL)
  }
  .d1 <-dim(.phi[[1]])[1]
  .ret <-vapply(seq_along(.phi), function(i) {
    .cov <- .phi[[i]]
    suppressWarnings(sqrt(diag(.cov))/unlist(.eta[i,, drop=FALSE])*100)
  }, double(.d1), USE.NAMES=FALSE)
  dim(.ret) <- c(.d1, length(.phi))
  dimnames(.ret) <- list(colnames(.phi[[1]]), names(.phi))
  .ret <- as.data.frame(t(.ret))
  names(.ret) <- paste0("rse(", names(.ret), ")%")
  .id <- seq_along(.phi)
  if (!is.null(names(.phi))) {
    .id <- names(.phi)
  }
  cbind(data.frame(ID=.id), .ret)
}
attr(nmObjGet.phiRSE, "desc") <- "relative standard error of each individual's eta (if present)"



#' @rdname nmObjGet
#' @export
nmObjGet.dataSav <- function(x, ...) {
  .obj <- x[[1]]
  .objEnv <- .obj$env
  if (exists("dataSav", .objEnv)) return(get("dataSav", envir=.objEnv))
  .data <- .obj$origData
  .env <- new.env(emptyenv())
  .foceiPreProcessData(.data, .env, .obj$ui)
  .env$dataSav
}
#attr(nmObjGet.dataSav, "desc") <- "data that focei sees for optimization"

#' @export
nmObjGet.foceiControl <- function(x, ...) {
  nmObjGetFoceiControl(.createEstObject(x[[1]]), ...)
}
attr(nmObjGet.foceiControl, "desc") <- "Get the focei control required for creating the nlmixr object"

#' @rdname nmObjGet
#' @export
nmObjGet.idLvl <- function(x, ...){
  .obj <- x[[1]]
  .objEnv <- .obj$env
  if (exists("idLvl", .objEnv)) return(get("idLvl", envir=.objEnv))
  .data <- .obj$origData
  .env <- new.env(emptyenv())
  .foceiPreProcessData(.data, .env, .obj$ui)
  .env$idLvl
}

#' @rdname nmObjGet
#' @export
nmObjGet.covLvl <- function(x, ...) {
  .obj <- x[[1]]
  .objEnv <- .obj$env
  if (exists("covLvl", .objEnv)) return(get("covLvl", envir=.objEnv))
  .data <- .obj$origData
  .env <- new.env(emptyenv())
  .foceiPreProcessData(.data, .env, .obj$ui)
  .env$covLvl
}
#attr(nmObjGet.dataSav, "desc") <- "data that focei sees for optimization"

.dataMergeStub <- function(obj) {
  .env      <- obj$env
  .origData <- obj$origData
  .origData$nlmixrRowNums <- seq_len(nrow(.origData))
  # add llikObs
  .llikObs <- FALSE
  if (exists("llikObs", obj$env)) {
    if (length(obj$env$llikObs) == length(.origData$nlmixrRowNums)) {
      .origData$nlmixrLlikObs <- obj$env$llikObs
      .llikObs <- TRUE
    } else {
      .dataSav <- obj$dataSav
      if (length(.dataSav$nlmixrRowNums) == length(obj$env$llikObs)) {
        .llik0 <- data.frame(nlmixrRowNums=.dataSav$nlmixrRowNums, nlmixrLlikObs=obj$env$llikObs)
        .llik0 <- .llik0[.llik0$nlmixrRowNums != 0,]
        .origData <- merge(.origData, .llik0, by="nlmixrRowNums", all.x=TRUE)
        .origData <- .origData[order(.origData$nlmixrRowNums),]
      } else {
        .nlmixrRowNums <- .dataSav[.dataSav$EVID == 0 | .dataSav$EVID == 2 |
                                     (.dataSav$EVID >= 9 & .dataSav$EVID <= 99),
                                   "nlmixrRowNums"]
        .llikObs <- obj$env$llikObs[!is.na(obj$env$llikObs)]
        if (length(.nlmixrRowNums) == length(.llikObs)) {
          .llik0 <- data.frame(nlmixrRowNums=.nlmixrRowNums, nlmixrLlikObs=.llikObs)
          .llik0 <- .llik0[.llik0$nlmixrRowNums != 0,]
          .origData <- merge(.origData, .llik0, by="nlmixrRowNums", all.x=TRUE)
          .origData <- .origData[order(.origData$nlmixrRowNums),]
        } else {
          warning("'nlmixrLlikObs' not added to dataset", call.=FALSE)
        }
      }
    }
  }
  .fitData <- as.data.frame(obj)
  if (is.null(.fitData$EVID)) .fitData$EVID <- 0
  if (is.null(.fitData$AMT))  .fitData$AMT  <- 0
  .names <- tolower(names(.origData))
  .wid <- which(.names == "id")
  names(.origData)[.wid] <- "ID"
  .fitData$nlmixrRowNums <- .env$.rownum
  .share <- setdiff(intersect(names(.origData), names(.fitData)), c("ID", "nlmixrRowNums"))
  .fitData <- .fitData[, !(names(.fitData) %in% .share)]
  if (inherits(.fitData$ID, "factor")) {
    .origData$ID <- factor(paste(.origData$ID), levels = levels(.fitData$ID))
  }
  list(.origData, .fitData)
}

#' @rdname nmObjGetData
#' @export
nmObjGetData.dataMergeLeft <- function(x, ...) {
  .obj <- x[[1]]
  .lst <- .dataMergeStub(.obj)
  .ret <- merge(.lst[[1]], .lst[[2]], by=c("ID", "nlmixrRowNums"), all.x=TRUE)
  .ret <- .ret[, names(.ret) != "nlmixrRowNums"]
  .ret
}
attr(nmObjGetData.dataMergeLeft, "desc") <- "left join between original and fit dataset"

#' @rdname nmObjGetData
#' @export
nmObjGetData.dataMergeRight <- function(x, ...) {
  .obj <- x[[1]]
  .lst <- .dataMergeStub(.obj)
  .ret <- merge(.lst[[1]], .lst[[2]], by=c("ID", "nlmixrRowNums"), all.y=TRUE)
  .ret <- .ret[, names(.ret) != "nlmixrRowNums"]
  .ret
}
attr(nmObjGetData.dataMergeRight, "desc") <- "right join between original and fit dataset"

#' @rdname nmObjGetData
#' @export
nmObjGetData.dataMergeInner <- function(x, ...) {
  .obj <- x[[1]]
  .lst <- .dataMergeStub(.obj)
  .ret <- merge(.lst[[1]], .lst[[2]], by=c("ID", "nlmixrRowNums"))
  .ret <- .ret[, names(.ret) != "nlmixrRowNums"]
  .ret
}
attr(nmObjGetData.dataMergeRight, "desc") <- "inner join between original and fit dataset"

#' @rdname nmObjGetData
#' @export
nmObjGetData.dataMergeFull <- function(x, ...) {
  .obj <- x[[1]]
  .lst <- .dataMergeStub(.obj)
  .ret <- merge(.lst[[1]], .lst[[2]], by=c("ID", "nlmixrRowNums"), all.x=TRUE, all.y=TRUE)
  .ret <- .ret[, names(.ret) != "nlmixrRowNums"]
  .ret
}
attr(nmObjGetData.dataMergeFull, "desc") <- "full join between original and fit dataset"


#' @rdname nmObjGet
#' @export
nmObjGet.saemTransformedData <- function(x, ...) {
  .dataSav <- nmObjGet.dataSav(x, ...)
  .obj <- x[[1]]
  .ui <- .obj$ui
  .saemGetDataForFit(.dataSav, .ui)
}
#attr(nmObjGet.saemTransformedData, "desc") <- "data that saem sees for optimization"


#' @rdname nmObjGet
#' @export
nmObjGet.parHistStacked <- function(x, ...) {
  .obj <- x[[1]]
  .env <- .obj$env
  if (exists("parHist", envir=.env)) {
    .parHist <- .env$parHist
    .iter <- .parHist$iter
    .ret <- data.frame(iter=.iter, stack(.parHist[, -1]))
    names(.ret) <- sub("values", "val",
                       sub("ind", "par", names(.ret)))
    .ret
  } else {
    NULL
  }
}
attr(nmObjGet.parHistStacked, "desc") <- "stacked parameter history"



#' @rdname nmObjGet
#' @export
nmObjGet.md5 <- function(x, ...) {
  .nlmixr2Md5(x[[1]])
}

#' @rdname nmObjGet
#' @export
nmObjGet.notes <- function(x, ...) {
  .notesFit(x[[1]])
}

#' @rdname nmObjGet
#' @export
nmObjGet.sigma <- function(x, ...) {
  .sigma(x[[1]])
}

#' @rdname nmObjGet
#' @export
nmObjGet.coefficients <- function(x, ...) {
  list(fixed = fixef(x[[1]]),
       random = ranef(x[[1]]))
}

#' @rdname nmObjGet
#' @export
nmObjGet.env <- function(x, ...) {
  x[[1]]
}

#' @rdname nmObjGet
#' @export
nmObjGet.condition <- function(x, ...) {
  .env <- x[[1]]
  .objDf <- .env$objDf
  if (any(names(.objDf) == "Condition#(Cor)")) {
    .cn <- .objDf[, "Condition#(Cor)"]
    .cn <- .cn[!is.na(.cn)]
    return(.cn)
  }
  return(NULL)
}

#' @rdname nmObjGet
#' @export
nmObjGet.simInfo <- function(x, ...) {
  .simInfo(x[[1]])
}

#' @rdname nmObjGet
#' @export
nmObjGet.seed <- function(x, ...) {
  .env <- x[[1]]
  if (exists("saem", .env)) {
    attr(get("saem", .env), "saem.cfg")$seed
  }
  NULL
}
#' @rdname nmObjGet
#' @export
nmObjGet.saemCfg <- function(x, ...) {
  .env <- x[[1]]
  if (exists("saem", .env)) {
    return(attr(get("saem", .env), "saem.cfg"))
  }
}

#' @export
nmObjGet.saemNmc <- function(x, ...) {
  .obj <- x[[1]]
  .env <- .obj$env
  if (exists("saemControl", envir=.env)) {
    .saemControl <- get("saemControl", envir=.env)
    return(.saemControl$mcmc$nmc)
  } else if (exists("control", envir=.env)) {
    .saemControl <- get("control", envir=.env)
    if (any(names(.saemControl) == "mcmc")) return(.saemControl$mcmc$nmc)
  }
  NA_integer_
}

#' @export
nmObjGet.saemEvtDf <- function(x, ...) {
  .obj <- x[[1]]
  .evt <- nmObjGet.dataSav(x, ...)
  .evt$ID <- .evt$ID - 1
  .evt
}
#attr(nmObjGet.saemEvtDf, "desc") <- "event data frame as seen by saem"

#' @export
nmObjGet.saemEvt <- function(x, ...) {
  as.matrix(nmObjGet.saemEvtDf(x, ...))
}
#attr(nmObjGet.saemEvtDf, "desc") <- "event matrix as seen by saem; stored in saem.cfg"

#' @export
nmObjGet.saemEvtMDf <- function(x, ...) {
  .nmc <- nmObjGet.saemNmc(x, ...)
  if (is.na(.nmc)) stop("cannot figure out the number of mcmc simulations", call.=FALSE)
  .evt <- nmObjGet.saemEvtDf(x, ...)
  .evtM <- .evt[rep(seq_len(dim(.evt)[1]), .nmc), ]
  .evtM$ID <- cumsum(c(FALSE, diff(.evtM$ID) != 0))
  .evtM
}
#attr(nmObjGet.saemEvtDf, "desc") <- "saem event data frame evtM for mcmc"

#' @export
nmObjGet.saemEvtM <- function(x, ...) {
  as.matrix(nmObjGet.saemEvtMDf(x, ...))
}
#attr(nmObjGet.saemEvtM, "desc") <- "saem event matrix evtM for mcmc; stored in saem.cfg"

#' @export
nmObjGet.saem <- function(x, ...) {
  .obj <- x[[1]]
  if (!exists("saem0", .obj)) return(NULL)
  .saem <- .obj$saem0
  .saemCfg <- attr(.saem, "saem.cfg")
  .saemCfg$evtM <- .obj$saemEvtM
  .saemCfg$evt <- .obj$saemEvt
  attr(.saem, "saem.cfg") <- .saemCfg
  .saem
}

#' @export
nmObjGet.innerModel <- function(x, ...) {
  .obj <- x[[1]]
  .env <- .obj$env
  if (exists("foceiModel", envir=.env)) {
    .model <- get("foceiModel", envir=.env)
  } else if (exists("model", envir=.env)) {
    .model <- get("model", envir=.env)
  } else {
    return(NULL)
  }
  .model$inner
}

#' @export
nmObjGet.innerModelForce <- function(x, ...) {
  .inner <- nmObjGet.innerModel(x, ...)
  if (is.null(.inner)) {
    .obj <- x[[1]]
    .env <- .obj$env
    .env$model <- .obj$ui$focei
    .inner <- .env$model$inner
    nmObjHandleModelObject(.env$model, .env)
  }
  .inner
}


#' @export
nmObjGet.ipredModel <- function(x, ...) {
  .obj <- x[[1]]
  .env <- .obj$env
  .est <- .env$est
  .env <- list(.env)
  class(.env) <- .est
  nmObjGetIpredModel(.env)
}
attr(nmObjGet.ipredModel, "desc") <- "rxode2 ipred model for fit"

#' @export
nmObjGet.predOnlyModel <- function(x, ...) {
  .obj <- x[[1]]
  .env <- .obj$env
  .est <- .env$est
  .env <- list(.env)
  class(.env) <- .est
  nmObjGetPredOnly(.env)
}
attr(nmObjGet.ipredModel, "desc") <- "rxode2 pred only model for fit"


#' Get the pred-only model for a fit depending on the object type
#'
#' By default it gets the focei models if available
#'
#' @param x nlmixr fit object
#'
#' @return rxode2 pred-only model
#'
#' @export
nmObjGetPredOnly <- function(x) {
  UseMethod("nmObjGetPredOnly")
}

#' @rdname nmObjGetPredOnly
#' @export
nmObjGetPredOnly.saem <- function(x) {
  .env <- x[[1]]
  .model <- NULL
  if (exists("saemModel", envir=.env)) {
    .model <- get("saemModel", envir=.env)
  } else if (exists("model", envir=.env)) {
    .model <- get("model", envir=.env)
  } else {
    stop("cannot find saem model components", call.=FALSE)
  }
  .model$predOnly
}

#' @rdname nmObjGetPredOnly
#' @export
nmObjGetPredOnly.default <- function(x) {
  .env <- x[[1]]
  .model <- NULL
  if (exists("foceiModel", envir=.env)) {
    .model <- get("foceiModel", envir=.env)
  } else if (exists("model", envir=.env)) {
    .model <- get("model", envir=.env)
  }
  .model$predOnly
}


#' Get the ipred model for a fit object depending on the object type
#'
#' By default it gets the focei models if available.
#'
#' @param x nlmixr fit object
#'
#' @return ipred `rxode2` model
#'
#' @export
nmObjGetIpredModel <- function(x) {
  UseMethod("nmObjGetIpredModel")
}

#' @rdname nmObjGetIpredModel
#' @export
nmObjGetIpredModel.saem <- function(x) {
  .env <- x[[1]]
  .model <- NULL
  if (exists("saemModel", envir=.env)) {
    .model <- get("saemModel", envir=.env)
  } else if (exists("model", envir=.env)) {
    .model <- get("model", envir=.env)
  } else {
    stop("cannot find saem model components", call.=FALSE)
  }
  .model$predOnly
}

#' @rdname nmObjGetIpredModel
#' @export
nmObjGetIpredModel.default <- function(x) {
  .env <- x[[1]]
  .model <- NULL
  if (exists("foceiModel", envir=.env)) {
    .model <- get("foceiModel", envir=.env)
  } else if (exists("model", envir=.env)) {
    .model <- get("model", envir=.env)
  }
  .inner <- .model$inner
  if (is.null(.inner)) return(.model$predOnly)
  .inner
}


#' @rdname nmObjGet
#' @export
nmObjGet.estimationModel <- function(x, ...) {
  .obj <- x[[1]]
  .env <- .obj$env
  .est <- .env$est
  .env <- list(.env)
  class(.env) <- .est
  nmObjGetEstimationModel(.env)
}
attr(nmObjGet.ipredModel, "desc") <- "rxode2 estimation model for fit"

#' Get the estimation model for a fit object depending on the object type
#'
#' By default it gets the focei models if available.
#'
#' @param x nlmixr fit object
#'
#' @return returns the estimation `$model` for the estimation type
#'
#' @export
nmObjGetEstimationModel <- function(x) {
  UseMethod("nmObjGetEstimationModel")
}

#' @rdname nmObjGetIpredModel
#' @export
nmObjGetEstimationModel.saem <- function(x) {
  .env <- x[[1]]
  attr(.env$saem, "saem.cfg")$opt$.rx
}

#' @rdname nmObjGetIpredModel
#' @export
nmObjGetEstimationModel.default <- function(x) {
  .env <- x[[1]]
  .model <- NULL
  if (exists("foceiModel", envir=.env)) {
    .model <- get("foceiModel", envir=.env)
  } else if (exists("model", envir=.env)) {
    .model <- get("model", envir=.env)
  }
  .inner <- .model$inner
  if (is.null(.inner)) return(.model$predOnly)
  .inner
}
#' Create an estimation object
#'
#' @param x nlmixr2 object
#' @return  list(nlmixr2 environment) with class of the estimation procedure ran.
#' @author Matthew L. Fidler
#' @noRd
.createEstObject <- function(x) {
  if (inherits(x, "nlmixr2FitData")) {
    .env <- x$env
  } else {
    .env <- x
  }
  if (exists("est", envir=.env)) {
    .est <- get("est", envir=.env)
    .ret <- list(.env)
    class(.ret) <- .est
    return(.ret)
  } else {
    stop("Cannot figure out the estimation method", call.=FALSE)
  }
}

#' @rdname nmObjGet
#' @export
nmObjGet.atol <- function(x, ...) {
  nmObjGetRxSolve(.createEstObject(x[[1]]), "atol")
}

#' @rdname nmObjGet
#' @export
nmObjGet.rtol <- function(x, ...) {
  nmObjGetRxSolve(.createEstObject(x[[1]]), "rtol")
}

#' @rdname nmObjGet
#' @export
nmObjGet.maxstepsOde <- function(x, ...) {
  nmObjGetRxSolve(.createEstObject(x[[1]]), "maxsteps")
}

#' @rdname nmObjGet
#' @export
nmObjGet.hmin <- function(x, ...) {
  nmObjGetRxSolve(.createEstObject(x[[1]]), "hmin")
}

#' @rdname nmObjGet
#' @export
nmObjGet.hmax <- function(x, ...) {
  nmObjGetRxSolve(.createEstObject(x[[1]]), "hmax")
}

#' @rdname nmObjGet
#' @export
nmObjGet.hini <- function(x, ...) {
  nmObjGetRxSolve(.createEstObject(x[[1]]), "hini")
}

#' @rdname nmObjGet
#' @export
nmObjGet.maxordn <- function(x, ...) {
  nmObjGetRxSolve(.createEstObject(x[[1]]), "maxordn")
}

#' @rdname nmObjGet
#' @export
nmObjGet.maxords <- function(x, ...) {
  nmObjGetRxSolve(.createEstObject(x[[1]]), "maxords")
}

#' @rdname nmObjGet
#' @export
nmObjGet.methodOde <- function(x, ...) {
  nmObjGetRxSolve(.createEstObject(x[[1]]), "method")
}

#' @rdname nmObjGet
#' @export
nmObjGet.covsInterpolation <- function(x, ...) {
  nmObjGetRxSolve(.createEstObject(x[[1]]), "covsInterpolation")
}

#' @rdname nmObjGet
#' @export
nmObjGet.control <- function(x, ...) {
  nmObjGetControl(.createEstObject(x[[1]]), ...)
}

#' Get an option for the estimation method
#'
#' By default it gets the focei models if available.
#'
#' @param x nlmixr fit object in a list.  The class is the
#'   estimation method used.
#' @param what What part of the rx solve are you attempting to get?
#' @return The estimation option based on `what`, for example
#'   `nlmixrObjGetRxSolve(x, "atol")` will get a double vector of
#'   absolute tolerances
#' @keywords internal
#' @export
nmObjGetRxSolve <- function(x, what) {
  UseMethod("nmObjGetRxSolve")
}

#' @rdname nmObjGetRxSolve
#' @export
nmObjGetRxSolve.default <- function(x, what) {
  .env <- x[[1]]
  .control <- nmObjGetControl(x)
  if (!any(names(.control) == "rxControl")) return(NULL)
  .lst <- .control$rxControl
  if (is.null(what)) return(.lst)
  .lst[[what]]
}

#' @rdname nmObjGet
#' @export
nmObjGet.simulationModel <- function(x, ...) {
  eval(rxode2::getBaseSimModel(x[[1]]))
}

#' @rdname nmObjGet
#' @export
nmObjGet.rxControl <- function(x, ...) {
  nmObjGetRxSolve(.createEstObject(x[[1]]), NULL)
}
attr(nmObjGet.rxControl, "desc") <- "rxode2 solving options"
