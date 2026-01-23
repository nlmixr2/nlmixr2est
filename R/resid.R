#'  theta/eta parameters needed for residuals/shrinkage calculations
#'
#' @param fit focei style fit
#'
#' @return list with:
#'
#'  - A rxode2 `params` dataset `pred` predictions
#'
#'  - A rxode2 `params` dataset for `ipred` predictions
#'
#'  - `eta.lst` is a numerical vector for each of the ETAs listed. The first 5 components are the mean, sd, variance, kurtosis and
#'    skewness statistics.  The rest of the components will be filled in later when calculating the shrinkage dataframe
#'
#' @author Matthew Fidler
#' @noRd
#'


# Since it can be accessed by the object, simply export it
#' @rdname nmObjGet
#' @export
nmObjGet.foceiThetaEtaParameters <- function(x, ...) {
  .fit <- x[[1]]
  .etas <- .fit$ranef
  .thetas <- .fit$fixef
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
  # The numeric versions are at
  # https://github.com/nlmixr2/rxode2/blob/7e27a7842ca0b5dd849ea75833bc7c34be729e31/R/rxsolve.R#L804,
  # but keeping them in sync will be fragile.  Only using the character
  # versions.
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
  # Fallback methods based on discussion in
  # https://github.com/nlmixr2/nlmixr2est/issues/254
  # Optimized: Use == for single element comparison instead of %in%
  if (currentOdeMethod == "dop853") {
    allOdeMethods <- "liblsoda"
  } else if (currentOdeMethod == "liblsoda" || currentOdeMethod == "lsoda") {
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
  while (recalc & length(odeMethods) > 0) {
    # Iterate through ODE methods
    recalcN <- 0
    currentOdeMethod <- odeMethods[[1]]
    odeMethods <- odeMethods[-1]
    .atol <- fit$atol[1]
    .rtol <- fit$rtol[1]
    ## message(currentOdeMethod)
    while (recalc & recalcN < fit$foceiControl$stickyRecalcN) {
      # Iterate up atol/rtol
      ## message("\t", .atol, " ", .rtol)
      .res <- .foceiSolveWithId(model, pars, fit$dataSav,
                                returnType = returnType,
                                atol = .atol, rtol = .rtol,
                                maxsteps = fit$maxstepsOde,
                                hmin = fit$hmin, hmax = fit$hmax, hini = fit$hini,
                                maxordn = fit$maxordn, maxords = fit$maxords,
                                method = rxode2::odeMethodToInt(currentOdeMethod),
                                keep=keep, addDosing=addDosing, subsetNonmem=subsetNonmem, addCov=addCov)
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
    .ipredModel <- fit$ipredModel
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
  .ret <- .Call(`_nlmixr2est_calcShrinkOnly`, .omega, thetaEtaParameters$eta.lst, length(fit$eta[,1]))
  .ret[, -dim(.omega)[1] - 1]
}

#' Add Levels to Data Based on Fit Object
#'
#' This function modifies a data frame by adding levels to a specified
#' variable based on the levels defined in a fit$ui$levels
#'
#' @param fit A list object that contains a `ui` element with `levels`
#'   to be added to the data.
#' @param data A data frame that will be modified by adding levels to
#'   one or more of its variables.
#' @return A modified data frame with levels added to the specified
#'   variables. If the variable's values are out of the defined range,
#'   they are set to `NA`.
#' @details The function checks if the `fit` object contains
#'   levels.
#'
#' If levels are present, it iterates through them and modifies the
#' corresponding variable in the data frame:
#'
#' - Converts the variable to integer type.
#'
#' - Sets values less than 1 to `NA_integer_`.
#'
#' - Sets values greater than the number of levels to `NA_integer_`.
#'
#' - Assigns the levels and sets the class of the variable to
#'   "factor".
#'
#'
#' @author Matthew L. Fidler
#'
#' @noRd
.addLevels <- function(fit, data) {
  .levels <-  fit$ui$levels
  if (!is.null(.levels)) {
    for (i in seq_along(.levels)) {
      .cur <- .levels[[i]] # language expression of levels() declaration
      .var <- deparse1(.cur[[2]][[2]]) # levels variable
      .w <- which(names(data) == .var) # does one of the output
                                       # variables match?
      if (length(.w) == 1) {
        # now change to a factor
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
#' @details
#'
#' If you ever want to add CWRES/FOCEi objective function you can use the \code{\link{addCwres}}
#'
#' If you ever want to add NPDE/EPRED columns you can use the \code{\link{addNpde}}
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
