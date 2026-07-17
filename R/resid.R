#' Get theta/eta parameters needed for residuals/shrinkage calculations
#'
#' @param fit focei style fit
#' @return list with rxode2 `params` datasets for `pred`/`ipred`, plus
#'   `eta.lst` (mean/sd/variance/kurtosis/skewness per ETA; more filled in
#'   later for shrinkage)
#' @author Matthew Fidler
#' @noRd


# Since it can be accessed by the object, simply export it
#' @rdname nmObjGet
#' @export
nmObjGet.foceiThetaEtaParameters <- function(x, ...) {
  .fit <- x[[1]]
  .etas <- .fit$ranef
  .thetas <- .fit$fixef
  .w <- which(names(.etas) %in% c("mixnum", "MIXEST"))
  if (length(.w) > 0L) {
    .etas <- .etas[, -.w, drop=FALSE]
  }
  .Call(`_nlmixr2est_nlmixr2Parameters`, .thetas, .etas)
}
#attr(nmObjGet.foceiThetaEtaParameters, "desc") <- "nmObjGet.foceiThetaEtaParameters"

#' This adjusts the names in the IPRED data frame to calculate censoring output correctly
#'
#' @param df ipred data frame
#'
#' @return ipred data frame with lower names dv, evid, cens, and limit
#'   in lower case (regardless of input)
#'
#' @author Matthew L. Fidler
#' @noRd
.residAdjustIpredNames <- function(df) {
  for (.v in c("dv", "evid", "cens", "limit")) {
    .w <- which(tolower(names(df)) == .v)
    if (length(.w) == 1L) {
      names(df)[.w] <- .v
    }
  }
  df
}


#' Solve making sure that ID is dropped
#'
#' This also suppresses the warning for id sorting
#' @param ... All parameters set to `rxode2::rxSolve()`
#' @param fit focei style fit
#' @return solved tataset
#' @author Matthew Fidler
#' @noRd
.foceiSolveWithId <- function(...) {
  .ret <- rxode2::rxSolve(..., warnIdSort = FALSE)
  if (names(.ret)[1] == "time") {
    ## For single subject ID is dropped.
    .ret <- data.frame(ID = 1L, .ret)
  }
  .w <- which(tolower(names(.ret)) == "dv")
  names(.ret)[.w] <- "dv"
  return(.ret)
}


#' Solve for pred/ipred types of calculations (including residuals)
#'
#' @param fit focei style fit
#' @param model rxode2 model
#' @param pars parameters to solve
#' @param keep vector of columns to keep
#' @param what string of what type of calculation is being performed
#' @inheritParams rxode2::rxSolve
#' @return Solved rxode2 data
#' @author Matthew Fidler
#' @noRd
.foceiSolvePars <- function(fit, model, pars=NULL, returnType="data.frame", keep=NULL, what="pred",
                            addDosing=FALSE, subsetNonmem=TRUE, addCov=FALSE) {
  if (is.null(model)) {
    stop("cannot solve with `model` NULL", call.=FALSE)
  }
  keep <- unique(c(keep, "nlmixrRowNums"))
  # Use character method names, not numeric codes, to avoid staying in sync
  # with rxode2 internals.
  currentOdeMethod <- fit$methodOde
  if (!inherits(currentOdeMethod, "character")) {
    cur <- as.integer(currentOdeMethod)+1L
    attr(cur, "levels") <- c("dop853", "lsoda", "liblsoda", "indLin")
    attr(cur, "class") <- "factor"
    currentOdeMethod <- as.character(cur)
  }
  allOdeMethods <-
    setdiff(
      eval(formals(rxode2::rxSolve)$method),
      # ignore indLin for now
      "indLin"
    )
  # Fallback ODE methods, see nlmixr2/nlmixr2est#254
  if (currentOdeMethod %in% "dop853") {
    allOdeMethods <- "liblsoda"
  } else if (currentOdeMethod %in% c("liblsoda", "lsoda")) {
    allOdeMethods <- "dop853"
  } # otherwise, use all the methods
  odeMethods <-
    append(
      list(currentOdeMethod),
      as.list(setdiff(allOdeMethods, currentOdeMethod))
    )
  failedMethods <- character()
  isFirstFit <- TRUE
  recalc <- TRUE
  maxAtolRtol <- fit$foceiControl$rxControl$maxAtolRtolFactor
  recalcFactor <- fit$foceiControl$odeRecalcFactor
  .tolFactor <- fit$env$tolFactor
  # For mixture models, pass per-subject mixture assignments via iCov so
  # rxode2 sets ind->mixest correctly during the table solve. `fit` is the
  # nlmixr2FitCore environment, accessed directly.
  .iCov <- NULL
  .env <- fit
  if (!is.environment(.env) && is.environment(fit$env)) {
    .env <- fit$env
  }
  if (is.environment(.env) && exists("mixIcov", envir=.env, inherits = FALSE)) {
    .iCov <- get("mixIcov", envir=.env, inherits = FALSE)
  }
  # Fallback flag: if rxode2 rejects iCov (older versions), retry without
  # it; .mixFixTable() post-corrects me/mn/mu from mixNum.
  .iCovOK <- !is.null(.iCov)
  while (recalc & length(odeMethods) > 0) {
    recalcN <- 0
    currentOdeMethod <- odeMethods[[1]]
    odeMethods <- odeMethods[-1]
    .atol <- fit$atol[1]
    .rtol <- fit$rtol[1]
    ## message(currentOdeMethod)
    while (recalc & recalcN < fit$foceiControl$stickyRecalcN) {
      # Iterate up atol/rtol
      ## message("\t", .atol, " ", .rtol)
      .res <- if (.iCovOK) {
        tryCatch(
          .foceiSolveWithId(model, pars, fit$dataSav,
                            returnType = returnType,
                            atol = .atol, rtol = .rtol,
                            maxsteps = fit$maxstepsOde,
                            hmin = fit$hmin, hmax = fit$hmax, hini = fit$hini,
                            maxordn = fit$maxordn, maxords = fit$maxords,
                            method = rxode2::odeMethodToInt(currentOdeMethod),
                            tolFactor = .tolFactor,
                            iCov = .iCov,
                            keep=keep, addDosing=addDosing, subsetNonmem=subsetNonmem, addCov=addCov),
          error = function(e) {
            if (grepl("time.varying|mixest must be", conditionMessage(e), ignore.case=TRUE)) {
              .iCovOK <<- FALSE
              .foceiSolveWithId(model, pars, fit$dataSav,
                                returnType = returnType,
                                atol = .atol, rtol = .rtol,
                                maxsteps = fit$maxstepsOde,
                                hmin = fit$hmin, hmax = fit$hmax, hini = fit$hini,
                                maxordn = fit$maxordn, maxords = fit$maxords,
                                method = rxode2::odeMethodToInt(currentOdeMethod),
                                tolFactor = .tolFactor,
                                iCov = NULL,
                                keep=keep, addDosing=addDosing, subsetNonmem=subsetNonmem, addCov=addCov)
            } else stop(e)
          })
      } else {
        .foceiSolveWithId(model, pars, fit$dataSav,
                          returnType = returnType,
                          atol = .atol, rtol = .rtol,
                          maxsteps = fit$maxstepsOde,
                          hmin = fit$hmin, hmax = fit$hmax, hini = fit$hini,
                          maxordn = fit$maxordn, maxords = fit$maxords,
                          method = rxode2::odeMethodToInt(currentOdeMethod),
                          tolFactor = .tolFactor,
                          iCov = NULL,
                          keep=keep, addDosing=addDosing, subsetNonmem=subsetNonmem, addCov=addCov)
      }
      rxode2::rxSolveFree()
      recalc <- any(is.na(.res$rx_pred_))
      recalcN <- recalcN + 1
      if (recalc) {
        .atol <- min(.atol*recalcFactor, maxAtolRtol)
        .rtol <- min(.rtol*recalcFactor, maxAtolRtol)
        if (.atol == maxAtolRtol && .rtol == maxAtolRtol) {
          recalcN <- fit$foceiControl$stickyRecalcN + 1
        }
      }
    }
    if (recalc) {
      failedMethods <- c(failedMethods, currentOdeMethod)
    }
    if (isFirstFit) {
      isFirstFit <- FALSE
      .resFirst <- .res
    }
  }
  if (recalc) {
    .res <- .resFirst
    warning("Problems solving ", what, " with ", paste(failedMethods, collapse = ", "), ", returning results from the first method")
  } else if (length(failedMethods) > 0) {
    warning("Problems solving ", what, " with ", paste(failedMethods, collapse = ", "), ", returning results from ", currentOdeMethod)
  }
  .res
}

#' Create a ipred/pred list from the focei style model
#'
#' @param fit focei style fit
#' @param thetaEtaParameters Theta/eta parameter list generated from `nmObjGet.foceiThetaEtaParameters()`
#' @param predOnly Pred Only for .ipred model (useful for mean/population models)
#' @inheritParams rxode2::rxSolve
#' @return list with ipred and pred datasets
#' @author Matthew Fidler
#' @noRd
.foceiPredIpredList <- function(fit, data=fit$dataSav,
                                thetaEtaParameters=fit$foceiThetaEtaParameters,
                                keep=NULL,
                                predOnly=is.null(fit$innerModel),
                                addDosing=FALSE, subsetNonmem=TRUE) {
  keep <- unique(c(keep, "nlmixrRowNums"))
  if (!predOnly && is.null(fit$innerModel)) {
    # Add inner problem calculation for cwres calculation
    fit$innerModelForce
  }
  .keep <- keep
  .names <- names(data)
  .lowerNames <- tolower(.names)
  for (.n in c("dv", "cens", "limit")) {
    .w <- which(.lowerNames == .n)
    if (length(.w) == 1L) .keep <- c(.keep, .names[.w])
  }
  .keep <- unique(.keep)
  .ipredModel <- fit$innerModel
  if (is.null(.ipredModel)) {
    predOnly <- TRUE
  }
  if (predOnly) {
    # Use predOnlyModel (which outputs the user-defined lhs like ka/cl/v and
    # tad/dosenum) instead of ipredModel (only rx_pred_/rx_r_ + sensitivities).
    # This keeps the output columns consistent with the cwres=TRUE path, which
    # solves predOnlyModel for these display columns; see nlmixr2/nlmixr2est#497.
    .ipredModel <- fit$predOnlyModel
  }
  .ret <- list(ipred = .residAdjustIpredNames(
    .foceiSolvePars(fit, .ipredModel, thetaEtaParameters$ipred,
                    returnType="data.frame.TBS", keep=.keep, what="ipred",
                    addDosing=addDosing, subsetNonmem=subsetNonmem, addCov=predOnly)),
               pred = .foceiSolvePars(fit, .ipredModel, thetaEtaParameters$pred,returnType="data.frame", what="pred",
                                      addDosing=addDosing, subsetNonmem=subsetNonmem),
               etaLst=thetaEtaParameters$eta.lst)
  if (!predOnly) {
    .ret <- c(.ret, list(predOnly=.foceiSolvePars(fit, fit$predOnlyModel, thetaEtaParameters$ipred,
                                                   returnType="data.frame", keep=.keep, what="ebe",
                                                   addDosing=addDosing, subsetNonmem=subsetNonmem, addCov=TRUE)))
  }
  .ret
}

.getRelevantLhs <- function(fit, keep=NULL, ipred=NULL) {
  .ret <- setdiff(fit$predOnlyModel$lhs,fit$ui$ini$name)
  .w <- which(regexpr("^rx", .ret) == -1)
  .ret <- unique(c(.ret[.w], keep))
  if (any(.ret == "tad")) {
    if (all(is.na(ipred$tad))) {
      .ret <- setdiff(.ret, c("tad", "dosenum"))
    }
  }
  .ret
}

.calcCwres0 <- function(fit, data=fit$dataSav, thetaEtaParameters=fit$foceiThetaEtaParameters,
                        table=tableControl(), dv=NULL, predOnly=FALSE,
                        addDosing=FALSE, subsetNonmem=TRUE, keep=NULL, npde=FALSE,
                        .prdLst) {
  assertNlmixrFit(fit)
  checkmate::assertDataFrame(data)
  checkmate::assertLogical(predOnly, len=1, any.missing=FALSE)
  checkmate::assertLogical(addDosing, len=1, any.missing=FALSE)
  checkmate::assertLogical(subsetNonmem, len=1, any.missing=FALSE)
  checkmate::assertLogical(npde, len=1, any.missing=FALSE)
  keep <- unique(c(keep, "nlmixrRowNums"))
  if (!inherits(dv, "numeric")) {
    dv <- .prdLst$ipred$dv
    table$doSim <- TRUE
  } else {
    table$doSim <- FALSE
  }

  if (npde) {
    .ni <- fit$dataNormInfo
    .sim <- vpcSim(fit, n = table$nsim, seed = table$seed,
                   addDosing=addDosing, subsetNonmem=subsetNonmem)
    .w <- which(names(.sim) == "ipred")
    if (length(.w) == 1) .sim <- .sim[, -.w]
    .w <- which(names(.sim) == "sim")
    .n0 <- c(names(.sim)[seq(1, .w)], "rxLambda", "rxYj", "rxLow", "rxHi")
    .sim <- .sim[, .n0]
    .ipred <- .prdLst$ipred
    .ipred <- .ipred[.ipred$nlmixrRowNums %in% .ni$nlmixrRowNums, ]
    .ipred <- .ipred[order(.ipred$nlmixrRowNums), ]
    .ret <- .Call(`_nlmixr2est_npdeCalc`, .sim, .ipred$dv, .ipred$evid,
                  .prdLst$ipred$cens, .prdLst$ipred$limit, table)
    .df <- data.frame(nlmixrRowNums=.prdLst$ipred$nlmixrRowNums)
    .ret1 <- .ret[[1]]
    .ret2 <- .ret[[2]]
    .ret2$nlmixrRowNums <- .ipred$nlmixrRowNums
    .ret2 <- merge(.df, .ret2, all.x=TRUE, by="nlmixrRowNums")
    .ret2 <- .ret2[,names(.ret2) != "nlmixrRowNums"]
    .ret1 <- as.data.frame(.ret1)
    .ret1$nlmixrRowNums <- .ipred$nlmixrRowNums
    .ret1 <- merge(.df, .ret1, all.x=TRUE, by="nlmixrRowNums")
    .ret1 <- .ret1[,names(.ret1) != "nlmixrRowNums"]
    .ret1 <- as.matrix(.ret1)
    list(.ret1, .ret2)
  } else {
    if (predOnly) {
      .state <- c(fit$predOnlyModel$state, fit$predOnlyModel$stateExtra)
      .lhs <- setdiff(unique(.getRelevantLhs(fit, keep, .prdLst$ipred)), .state)
      .params <- setdiff(intersect(names(fit$dataSav),fit$predOnlyModel$params),
                         c("CMT","cmt","Cmt", .state, .lhs))
      .Call(`_nlmixr2est_resCalc`, .prdLst, fit$omega,
            fit$eta, .prdLst$ipred$dv, .prdLst$ipred$evid, .prdLst$ipred$cens,
            .prdLst$ipred$limit, .lhs, .state, .params, fit$IDlabel, table)
    } else {
      .state <- c(fit$predOnlyModel$state, fit$predOnlyModel$stateExtra)
      .lhs <- setdiff(unique(.getRelevantLhs(fit, keep, .prdLst$predOnly)), .state)
      .params <- setdiff(intersect(names(fit$dataSav),fit$predOnlyModel$params),c("CMT","cmt","Cmt", .state, .lhs))
      .Call(`_nlmixr2est_cwresCalc`, .prdLst, fit$omega,
            fit$eta, .prdLst$ipred$dv, .prdLst$ipred$evid, .prdLst$ipred$cens,
            .prdLst$ipred$limit, .lhs, .state, .params, fit$IDlabel, table)
    }
  }
}

.calcCwres <- function(fit, data=fit$dataSav, thetaEtaParameters=fit$foceiThetaEtaParameters,
                       table=tableControl(), dv=NULL, predOnly=TRUE,
                       addDosing=FALSE, subsetNonmem=TRUE, keep=NULL, npde=FALSE,
                       .prdLst=NULL) {
  if (!inherits(table, "tableControl")) table <- do.call(tableControl, table)
  keep <- unique(c(keep, "nlmixrRowNums"))
  if (is.null(.prdLst)) {
    .prdLst <- .foceiPredIpredList(fit, data=data, keep=keep, thetaEtaParameters=thetaEtaParameters, predOnly=predOnly,
                                   addDosing=addDosing, subsetNonmem=subsetNonmem)
  }
  ## Split out so that .prdLst can be shared between npde/cwres npde/res
  .ret <- .calcCwres0(fit, data, thetaEtaParameters, table, dv=dv, predOnly,
                      addDosing, subsetNonmem, keep, npde, .prdLst=.prdLst)
  .dups <- which(duplicated(names(.ret)))
  if (length(.dups) > 0) {
    warning("some duplicate columns were dropped", call.=FALSE)
    .ret <- .ret[, -.dups]
  }
  .ret
}

.calcRes <- function(..., predOnly=TRUE) {
  .calcCwres(..., predOnly=predOnly)
}

.calcNpde <- function(..., npde=TRUE, predOnly=TRUE) {
  .calcCwres(..., npde=npde, predOnly=predOnly)
}

.calcIres <- function(fit, data=fit$dataSav, table=tableControl(), dv=NULL,
                      addDosing=FALSE, subsetNonmem=TRUE, keep=NULL) {
  keep <- unique(c(keep, "nlmixrRowNums"))
  if (!inherits(table, "tableControl")) table <- do.call(tableControl, table)
  .keep <- keep
  .names <- names(data)
  .lowerNames <- tolower(.names)
  for (.n in c("dv", "cens", "limit")) {
    .w <- which(.lowerNames == .n)
    if (length(.w) == 1L) .keep <- c(.keep, .names[.w])
  }
  .thetas <- fit$fixef
  names(.thetas) <- paste0("THETA[", seq_along(.thetas), "]")
  .eta <- fit$eta
  if (inherits(.eta, "data.frame")) {
    .n <- length(.eta) - 1
    .thetas <- c(.thetas, setNames(rep(0, .n), paste0("ETA[", seq_len(.n), "]")))
  }
  .pars <- fit$ipredModel$params
  .cmt <- which(tolower(.pars) == "cmt")
  if (length(.cmt) == 1) {
    .cmt <-.pars[.cmt]
    .keep <- c(.cmt, .keep)
  }
  .ipred <- .residAdjustIpredNames(.foceiSolvePars(fit, fit$ipredModel, .thetas,
                                                   returnType="data.frame.TBS", keep=.keep, what="ipred",
                                                   addDosing=addDosing, subsetNonmem=subsetNonmem))
  if (!inherits(dv, "numeric")) {
    dv <- .ipred$dv
    table$doSim <- TRUE
  } else {
    table$doSim <- FALSE
  }
  .state <- c(fit$ipredModel$state, fit$ipredModel$stateExtra)
  .lhs <- setdiff(unique(.getRelevantLhs(fit, keep, .ipred)), .state)
  .params <- setdiff(intersect(names(fit$dataSav),fit$ipredModel$params),c("CMT","cmt","Cmt", .state, .lhs))
  .ret <- .Call(`_nlmixr2est_iresCalc`, .ipred, dv, .ipred$evid, .ipred$cens, .ipred$limit,
                .lhs, .state, .params, fit$IDlabel, table)
  .dups <- which(duplicated(names(.ret)))
  if (length(.dups) > 0) {
    warning("some duplicate columns were dropped", call.=FALSE)
    .ret <- .ret[, -.dups]
  }
  .addLevels(fit, .ret)
}

.calcShrinkOnly <- function(fit, thetaEtaParameters=fit$foceiThetaEtaParameters) {
  .omega <- fit$omega
  if (exists("etaExpected", envir=fit$env)) {
    .etas <- fit$env$etaExpected
    .w <- which(names(.etas) %in% c("mixnum", "MIXEST"))
    if (length(.w) > 0L) {
      .etas <- .etas[, -.w, drop=FALSE]
    }
    .pars <- .Call(`_nlmixr2est_nlmixr2Parameters`, fit$fixef, .etas)
    .ret <- .Call(`_nlmixr2est_calcShrinkOnly`, .omega, .pars$eta.lst, length(.etas[,1]))
  } else {
    .ret <- .Call(`_nlmixr2est_calcShrinkOnly`, .omega, thetaEtaParameters$eta.lst, length(fit$eta[,1]))
  }
  .ret[, -dim(.omega)[1] - 1]
}

#' Add factor levels to data based on fit object
#'
#' @param fit A list with a `ui` element containing `levels` to apply.
#' @param data A data frame to modify.
#' @return Data frame with the specified variables converted to factors;
#'   out-of-range values become `NA`.
#' @author Matthew L. Fidler
#' @noRd
.addLevels <- function(fit, data) {
  .levels <-  fit$ui$levels
  if (!is.null(.levels)) {
    for (i in seq_along(.levels)) {
      .cur <- .levels[[i]] # levels() expression
      .var <- deparse1(.cur[[2]][[2]]) # levels variable name
      .w <- which(names(data) == .var) # does an output var match?
      if (length(.w) == 1) {
        # convert to factor
        data[[.var]] <- as.integer(data[[.var]])
        .w <- which(data[[.var]] < 1L)
        data[[.var]][.w] <- NA_integer_
        .lvls <- eval(.cur[[3]])
        .w <- which(data[[.var]] > length(.lvls))
        data[[.var]][.w] <- NA_integer_
        attr(data[[.var]], "levels") <- .lvls
        attr(data[[.var]], "class") <- "factor"
      }
    }
  }
  data
}

.calcTables <- function(fit, data=fit$dataSav, thetaEtaParameters=fit$foceiThetaEtaParameters,
                        table=tableControl(), keep=NULL) {
  keep <- unique(c(keep, "nlmixrRowNums"))

  if (!inherits(table, "tableControl")) table <- do.call(tableControl, table)
  if (is.null(table$cwres)) {
    table$cwres <- !is.null(fit$innerModel)
  }
  if (table$cwres) {
    fit$innerModelForce
  }
  if (is.null(table$npde)) {
    table$npde <- FALSE
  }
  .predOnly <- !table$cwres
  .censMethod <- table$censMethod
  .ret <- vector("list",2)
  .thetaEtaParameters <- fit$foceiThetaEtaParameters
  .prdLst <- .foceiPredIpredList(fit, data=fit$dataSav, keep=keep, thetaEtaParameters=.thetaEtaParameters, predOnly=.predOnly,
                                 addDosing=table$addDosing, subsetNonmem=table$subsetNonmem)
  if (.censMethod %in% c(2L, 6L)) {
    if (!table$npde) {
      warning("censoring method requires npde, adding npde", call.=FALSE)
      table$npde <- TRUE
    }
    .npde1 <- TRUE
    .npde2 <- FALSE
  } else {
    .npde1 <- FALSE
    .npde2 <- TRUE
  }
  if ((.npde1 & table$npde) | !.npde1)
    .ret[[1]] <- .calcCwres(fit, data=fit$dataSav, thetaEtaParameters=.thetaEtaParameters, table=table,
                            predOnly=.predOnly, addDosing=table$addDosing, subsetNonmem=table$subsetNonmem,
                            keep=keep, .prdLst=.prdLst, npde=.npde1)
  if ((.npde2 & table$npde) | !.npde2)
    .ret[[2]] <- .calcCwres(fit, data=fit$dataSav, thetaEtaParameters=.thetaEtaParameters, table=table, dv=.ret[[1]][[1]],
                            predOnly=.predOnly, addDosing=table$addDosing, subsetNonmem=table$subsetNonmem,
                            keep=keep, .prdLst=.prdLst, npde=.npde2)
  .ret <- .Call(`_nlmixr2est_popResFinal`, .ret)
  .dups <- which(duplicated(names(.ret)))
  if (length(.dups) > 0) {
    warning("some duplicate columns were dropped", call.=FALSE)
    .ret <- .ret[, -.dups]
  }
  .ret[[1]] <- .addLevels(fit, .ret[[1]])
  .ret
}

#' Re-insert subjects dropped during preprocessing into an output table
#'
#' Adds subjects dropped by `.foceiPreProcessData()` back into the table with a population `PRED` and NA individual columns.
#'
#' @param df assembled output data.frame (one row per output record for the
#'   subjects that were estimated)
#' @param object the nlmixr2 fit object
#' @param table the `tableControl()` list
#' @return `df` with the dropped subjects re-inserted (population `PRED`, NA
#'   individual columns), `ID` factor extended to the full `origData` levels, and
#'   rows re-sorted by ID/TIME.  Unchanged when nothing was dropped.
#' @author Matthew L. Fidler
#' @noRd
.reinsertNoObsSubjects <- function(df, object, table) {
  .orig <- object$origData
  if (is.null(.orig) || is.null(.orig$ID) || is.null(df$ID)) {
    return(df)
  }
  # normalize origData names to match the output table (uppercase except covariates)
  .cov <- tryCatch(object$ui$covariates, error = function(e) character(0))
  names(.orig) <- .nmUpcaseNonCov(names(.orig), .cov)
  .have <- levels(df$ID)
  if (is.null(.have)) .have <- unique(as.character(df$ID))
  .full <- unique(as.character(.orig$ID))
  .miss <- setdiff(.full, .have)
  if (length(.miss) == 0L) {
    return(df)
  }
  .rows <- .orig[as.character(.orig$ID) %in% .miss, , drop = FALSE]
  # only observation/other records (EVID 0 or 2), plus dosing if addDosing is requested
  if (!is.null(.rows$EVID) && !isTRUE(table$addDosing)) {
    .rows <- .rows[.rows$EVID %in% c(0, 2), , drop = FALSE]
  }
  if (nrow(.rows) == 0L) {
    return(df)
  }
  # NA rows with output column types, filled from origData below
  .add <- df[rep(1L, nrow(.rows)), , drop = FALSE]
  for (.cn in names(.add)) {
    .add[[.cn]][] <- NA
  }
  for (.cn in intersect(names(.add), names(.rows))) {
    if (.cn == "ID") next
    if (is.factor(df[[.cn]])) {
      .add[[.cn]] <- factor(as.character(.rows[[.cn]]), levels = levels(df[[.cn]]))
    } else {
      .add[[.cn]] <- as.vector(.rows[[.cn]])
    }
  }
  .add <- .fillNoObsPred(.add, .rows, object, table)
  # rebuild ID factor on both halves with the full origData level order for rbind
  .add$ID <- factor(as.character(.rows$ID), levels = .full)
  df$ID <- factor(as.character(df$ID), levels = .full)
  .out <- rbind(df, .add)
  .out <- .out[order(as.integer(.out$ID), .out$TIME), , drop = FALSE]
  rownames(.out) <- NULL
  .out
}
#' Fill the population PRED column for re-inserted no-observation subjects
#'
#' Solves at the population estimates (eta = 0) and writes `PRED` into `add`; best-effort, leaves `PRED` NA on failure.
#'
#' @param add the NA-filled re-insertion rows (output-table columns)
#' @param rows the matching origData records for the dropped subjects
#' @param object the nlmixr2 fit object
#' @param table the `tableControl()` list
#' @return `add` with `PRED` filled where it could be solved
#' @author Matthew L. Fidler
#' @noRd
.fillNoObsPred <- function(add, rows, object, table) {
  if (!("PRED" %in% names(add)) || is.null(rows$ID) || is.null(rows$TIME)) {
    return(add)
  }
  .orig <- object$origData
  if (is.null(.orig) || is.null(.orig$ID)) return(add)
  .cov <- tryCatch(object$ui$covariates, error = function(e) character(0))
  names(.orig) <- .nmUpcaseNonCov(names(.orig), .cov)
  .dvCol <- which(names(.orig) == "DV")
  # solve one subject at a time to align rows by order and avoid ID renumbering;
  # dummy DV avoids the predict run dropping these subjects again
  for (.id in unique(as.character(rows$ID))) {
    .subj <- .orig[as.character(.orig$ID) == .id, , drop = FALSE]
    if (length(.dvCol) == 1L) {
      .subj[[.dvCol]][is.na(.subj[[.dvCol]])] <- 0
    }
    .p <- tryCatch(suppressMessages(suppressWarnings(
      as.data.frame(stats::predict(object, newdata = .subj, level = "population")))),
      error = function(e) NULL)
    if (is.null(.p)) next
    .pc <- which(tolower(names(.p)) == "pred")
    if (length(.pc) != 1L) next
    .w <- which(as.character(rows$ID) == .id)
    if (length(.w) == nrow(.p)) {
      add$PRED[.w] <- .p[[.pc]]
    } else {
      .tc <- which(tolower(names(.p)) == "time")
      if (length(.tc) == 1L) {
        add$PRED[.w] <- .p[[.pc]][match(rows$TIME[.w], .p[[.tc]])]
      }
    }
  }
  add
}
#' Add table information to nlmixr2 fit object without tables
#'
#' @param object nlmixr2 family of objects
#' @param updateObject Update the object (default FALSE)
#' @param data Saved data from
#' @param thetaEtaParameters Internal theta/eta parameters
#' @param table a `tableControl()` list of options
#' @param keep Character Vector of items to keep
#' @param drop Character Vector of items to drop or NULL
#' @param envir Environment to search for updating
#' @return Fit with table information attached
#' @author Matthew Fidler
#' @export
#' @examples
#'
#' \donttest{
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
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' # run without tables step
#' f <- nlmixr2(one.cmt, theo_sd, "saem", control=list(calcTables=FALSE))
#'
#' print(f)
#'
#' # Now add the tables
#'
#' f <- addTable(f)
#'
#' print(f)
#'
#' }
addTable <- function(object, updateObject = FALSE,
                     data=object$dataSav,
                     thetaEtaParameters=object$foceiThetaEtaParameters,
                     table=tableControl(),
                     keep=NULL, drop=NULL,
                     envir = parent.frame(1)) {
  nlmixr2global$finalUiCompressed <- FALSE
  on.exit(nlmixr2global$finalUiCompressed <- TRUE)
  nlmixrWithTiming("table", {
    keep <- unique(c(keep, "nlmixrRowNums"))
    .malert("Calculating residuals/tables")
    .objName <- substitute(object)
    if (!inherits(object, "nlmixr2FitCore")) {
      stop("requires a nlmixr2 fit object",
           call.=FALSE)
    }
    .fit <- object$env
    if (exists("origControl", .fit)) {
      .control <- .fit$origControl
    } else if (exists("control", .fit)) {
      .control <- .fit$control
    } else {
      .control <- foceiControl()
    }
    if (is.null(.fit$omega)) {
      .df <- .calcIres(.fit, data=data, table=table, dv=NULL,
                       addDosing=table$addDosing, subsetNonmem=table$subsetNonmem, keep=keep)
    } else {
      .tabs <- .calcTables(.fit, data=data, table=table, keep=keep)
      assign("shrink", .tabs$shrink, .fit)
      .df <- .tabs$resid
    }
    .rownum <- as.integer(.df$nlmixrRowNums)
    assign(".rownum", .rownum, envir=.fit)
    drop <- c(drop, "rxLambda", "rxYj", "nlmixrRowNums",
              "rx__sens_central_BY_p1",
                "rx__sens_central_BY_v1",
                "rx__sens_central_BY_p2",
                "rx__sens_central_BY_p3",
                "rx__sens_central_BY_p4",
                "rx__sens_central_BY_ka",
                "rx__sens_peripheral1_BY_p1",
                "rx__sens_peripheral1_BY_v1",
                "rx__sens_peripheral1_BY_p2",
                "rx__sens_peripheral1_BY_p3",
                "rx__sens_peripheral1_BY_p4",
                "rx__sens_peripheral1_BY_ka",
                "rx__sens_peripheral2_BY_p1",
                "rx__sens_peripheral2_BY_v1",
                "rx__sens_peripheral2_BY_p2",
                "rx__sens_peripheral2_BY_p3",
                "rx__sens_peripheral2_BY_p4",
                "rx__sens_peripheral2_BY_ka",
                "rx__sens_depot_BY_ka")
    .w <- -which(names(.df) %in% drop)
    if (length(.w) > 0) .df <- .df[, .w, drop=FALSE]
    class(.df) <- "data.frame"
    .id <- .df$ID
    if (is.null(.id)) {
      .df$ID <- 1L
    } else {
      attr(.id, "levels") <- object$idLvl
      class(.id) <- "factor"
      .df$ID <- .id
    }
    .covLvl <- object$covLvl
    for (.v in names(.covLvl)) {
      .l <- as.integer(.df[[.v]])
      attr(.l, "levels") <- .covLvl[[.v]]
      class(.l) <- "factor"
      .df[[.v]] <- .l
    }
    # re-insert subjects dropped for having no usable observation
    .df <- .reinsertNoObsSubjects(.df, object, table)
    .isDplyr <- requireNamespace("tibble", quietly = TRUE)
    if (!.isDplyr) {
      .isDataTable <- requireNamespace("data.table", quietly = TRUE)
      if (.isDataTable) {
        .df <- data.table::data.table(.df)
      }
    } else {
      .df <- tibble::as_tibble(.df)
    }
    .cls <- class(.df)
    if (!any(names(.control) == "interaction")) {
      .control$interaction <- FALSE
    }
    if (.fit$method == "population only") {
      .cls <- c("nlmixr2FitData", "nlmixr2FitCore", "pop", paste0("nlmixr2.", .fit$env$est),  .cls)
    } else {
      .cls <- c("nlmixr2FitData", "nlmixr2FitCore", paste0("nlmixr2.", .fit$env$est), .cls)
    }
    if (inherits(updateObject, "logical")) {
      if (!updateObject) {
        .fit <- .cloneEnv(.fit)
      }
    }
    class(.fit) <- "nlmixr2FitCoreSilent"
    attr(.cls, ".foceiEnv") <- .fit
    class(.df) <- .cls
    if (inherits(updateObject, "logical")) {
      if (updateObject) {
        .parent <- envir
        .bound <- do.call("c", lapply(ls(.parent, all.names = TRUE), function(.cur) {
          if (.cur == .objName && identical(.parent[[.cur]]$env, .fit$env)) {
            return(.cur)
          }
          return(NULL)
        }))
        if (length(.bound) == 1) {
          if (exists(.bound, envir = .parent)) {
            assign(.bound, .df, envir = .parent)
          }
        }
      }
    }
    .msuccess("done")
    .df
  }, envir=object)
}

#' Output table/data.frame options
#'
#' @param npde When TRUE, request npde regardless of the algorithm used.
#'
#' @param cwres When TRUE, request CWRES and FOCEi likelihood
#'     regardless of the algorithm used.
#'
#' @param censMethod Handle censoring method:
#'
#'  - `"truncated-normal"` Simulates from a truncated normal distribution under the assumption of the model and censoring.
#'
#'  - `"cdf"` Use the cdf-method for censoring with npde and use this for any other residuals (`cwres` etc)
#'
#'  - `"omit"` omit the residuals for censoring
#'
#' @param ties When `TRUE` jitter prediction-discrepancy points to discourage ties in cdf.
#'
#' @param cholSEtol The tolerance for the `rxode2::choleSE` function
#'
#' @param eta is a Boolean indicating if `eta` values will be included (default `TRUE`)
#'
#' @param state is a Boolean indicating if `state` values will be included (default `TRUE`)
#'
#' @param lhs is a Boolean indicating if remaining `lhs` values will be included (default `TRUE`)
#'
#' @param covariates is a Boolean indicating if covariates will be included (default `TRUE`)
#'
#' @param keep is the keep sent to the table
#'
#' @param drop is the dropped variables sent to the table
#'
#' @inheritParams addNpde
#' @inheritParams rxode2::rxSolve
#'
#' @details Use \code{\link{addCwres}} to add CWRES/FOCEi objective
#'   function, or \code{\link{addNpde}} to add NPDE/EPRED columns.
#'
#' @return A list of table options for nlmixr2
#' @author Matthew L. Fidler
#' @export
tableControl <- function(npde = NULL,
                         cwres = NULL,
                         nsim = 300, ties = TRUE,
                         censMethod=c("truncated-normal", "cdf", "ipred", "pred", "epred", "omit"),
                         seed = 1009,
                         cholSEtol=(.Machine$double.eps)^(1/3),
                         state=TRUE,
                         lhs=TRUE,
                         eta=TRUE,
                         covariates=TRUE,
                         addDosing=FALSE, subsetNonmem = TRUE,
                         cores=NULL,
                         keep=NULL,
                         drop=NULL) {
  checkmate::assertLogical(npde, any.missing=FALSE, len=1, null.ok=TRUE)
  checkmate::assertLogical(cwres, any.missing=FALSE, len=1, null.ok=TRUE)
  checkmate::assertLogical(ties, any.missing=FALSE, len=1, null.ok=FALSE)
  checkmate::assertIntegerish(nsim, lower=0, len=1)
  checkmate::assertIntegerish(seed, lower=0, len=1)
  checkmate::assertNumeric(cholSEtol, lower=0, len=1)
  checkmate::assertLogical(state, len=1, any.missing=FALSE)
  checkmate::assertLogical(lhs, len=1, any.missing=FALSE)
  checkmate::assertLogical(eta, len=1, any.missing=FALSE)
  checkmate::assertLogical(covariates, len=1, any.missing=FALSE)
  checkmate::assertLogical(addDosing, len=1, any.missing=FALSE)
  checkmate::assertLogical(subsetNonmem, len=1, any.missing=FALSE)
  checkmate::assertCharacter(keep, null.ok=TRUE, pattern = "^[.]*[a-zA-Z]+[a-zA-Z0-9._]*$",
                             any.missing = FALSE,min.chars=1)
  .invalidKeep <- c("id", "sim.id", "resetno", "time", "nlmixrRowNums")
  .invalidKeep <- intersect(tolower(keep), tolower(.invalidKeep))
  if (length(.invalidKeep) > 0) {
    .w <- which(tolower(keep) %in% .invalidKeep)
    keep <- keep[-.w]
    warning("'keep' contains ", paste(.invalidKeep, collapse=", "), "\nwhich are output when needed, ignoring these items", call.=FALSE)
  }
  .invalidKeep <- c("evid",  "ss", "amt", "rate", "dur", "ii")
  .invalidKeep <- intersect(tolower(keep), tolower(.invalidKeep))
  if (length(.invalidKeep) > 0) {
    stop("'keep' cannot contain ", paste(.invalidKeep, collapse=", "), "\nconsider using addDosing=TRUE or merging to original dataset\nfor a fit the merge can be called by fit$dataMergeLeft fit$dataMergeRight or fit$dataMergeInner", call.=FALSE)
  }
  .invalidKeep <- c ("rxLambda", "rxYj", "rxLow", "rxHi")
  .invalidKeep <- intersect(tolower(keep), tolower(.invalidKeep))
  if (length(.invalidKeep) > 0) {
    stop("'keep' cannot contain ", paste(.invalidKeep, collapse=", "), call.=FALSE)
  }

  checkmate::assertCharacter(drop, null.ok=TRUE)
  if (inherits(censMethod, "character")) {
    .censMethod <- setNames(c("truncated-normal"=3L, "cdf"=2L, "omit"=1L, "pred"=5L, "ipred"=4L, "epred"=6L)[match.arg(censMethod)], NULL)
  } else {
    checkmate::assertIntegerish(censMethod)
    .censMethod <- as.integer(censMethod)
  }
  if (is.null(cores)) {
    cores <- rxode2::rxCores()
  } else {
    checkmate::assertIntegerish(cores, len=1, lower=1)
  }
  .ret <- list(
    npde = npde, cwres = cwres, nsim = nsim, ties = ties, seed = seed,
    censMethod=.censMethod,
    cholSEtol=cholSEtol, state=state, lhs=lhs, eta=eta, covariates=covariates, addDosing=addDosing, subsetNonmem=subsetNonmem, cores=cores, keep=keep, drop=drop)
  class(.ret) <- "tableControl"
  .ret
}

#' @export
rxUiDeparse.tableControl <- function(object, var) {
  .default <- tableControl()
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}
