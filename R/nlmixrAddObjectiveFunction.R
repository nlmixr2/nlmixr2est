#'  Add objective function data frame to the current objective function
#'
#' @param fit nlmixr fit object
#' @param objDf nlmixr objective function data frame which has column
#'   names "OBJF", "AIC", "BIC", "Log-likelihood" and
#'   "Condition#(Cov)" "Condition#(Cor)"
#' @param type Objective Function Type
#' @param etaObf Eta objective function table to add (with focei) to
#'   give focei objective function
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @export
nlmixrAddObjectiveFunctionDataFrame <- function(fit, objDf, type, etaObf=NULL) {
  assertNlmixrFit(fit)
  checkmate::assertCharacter(type, len=1, any.missing=FALSE)
  .inRow <- assertNlmixrObjDataFrameRow(objDf, allowNa=FALSE)
  .cur <- fit$objDf
  .rownames <- row.names(.cur)
  if (!is.null(etaObf)) {
    assign("etaObf", etaObf, envir=fit$env)
  }
  if (length(.cur$OBJF) == 1) {
    .inRow2 <- assertNlmixrObjDataFrameRow(.cur, allowNa=TRUE)
    .cn <- NA_real_
    if (!is.na(.inRow2[[2]])) {
      .cn <- .inRow2[[2]]
    } else if (!is.na(.inRow[[2]])) {
      .cn <- .inRow[[2]]
    }
    .cnr <- NA_real_
    if (!is.na(.inRow2[[3]])) {
      .cnr <- .inRow2[[3]]
    } else if (!is.na(.inRow[[3]])) {
      .cnr <- .inRow[[3]]
    }
    if (is.na(.inRow2[[1]][[1]])) {
      # Here the original data frame is NA, that is the objective function has not been calculated
      .tmp <- cbind(.inRow[[1]], data.frame("Condition#(Cov)"=.cn, "Condition#(Cor)"=.cnr, check.names=FALSE))
      row.names(.tmp) <- type
      assign("objDf", .tmp, envir=fit$env)
      setOfv(fit, type)
    } else {
      if (any(.rownames == type)) {
        stop("objective function '", type, "' already present", call.=FALSE)
      }
      # Now the original data frame is not NA.
      .tmp <- rbind(.inRow[[1]], .inRow2[[1]])
      .tmp[["Condition#(Cov)"]] <- .cn
      .tmp[["Condition#(Cor)"]] <- .cnr
      row.names(.tmp) <- c(type, .rownames)
      assign("objDf", .tmp, envir=fit$env)
      setOfv(fit, type)
    }
  } else {
    if (any(.rownames == type)) stop("objective function '", type, "' already present", call.=FALSE)
    ## Now there is at least one interesting objective function
    .cn <- .cur[["Condition#(Cov)"]][1]
    if (is.null(.cn)) {
      .cn <- NA_real_
    }
    if (is.na(.cn) & !is.na(.inRow[[2]])) {
      .cn <- .inRow[[2]]
    }
    .cnr <- .cur[["Condition#(Cor)"]][1]
    if (is.null(.cnr)) {
      .cnr <- NA_real_
    }
    if (is.na(.cnr) & !is.na(.inRow[[3]])) {
      .cnr <- .inRow[[3]]
    }
    .cur <- rbind(.cur[, !(names(.cur) %in% c("Condition#(Cor)", "Condition#(Cov)"))], .inRow[[1]])
    .cur[["Condition#(Cov)"]] <- .cn
    .cur[["Condition#(Cor)"]] <- .cnr
    row.names(.cur) <- c(.rownames, type)
    assign("objDf", .cur, envir=fit$env)
    setOfv(fit, type)
  }
}
