#' Setup the par history information
#'
#' @param .ret Return data
#' @return Nothing called for side effects
#' @author Matthew L. Fidler
#' @noRd
.nlmixr2setupParHistData <- function(.ret) {
  if (exists("parHistData", .ret)) {
    .tmp <- .ret$parHistData
    .tmp <- .tmp[.tmp$type == "Unscaled", names(.tmp) != "type"]
    .iter <- .tmp$iter
    .tmp <- .tmp[, names(.tmp) != "iter"]
    .ret$parHistStacked <- data.frame(stack(.tmp), iter = .iter)
    names(.ret$parHistStacked) <- c("val", "par", "iter")
    .ret$parHist <- data.frame(iter = .iter, .tmp)
  }
}

.getBackTransformationFunction <- function(par, ui) {
  # This has a specified back-transformation
  .w <- which(ui$iniDf$name == par)
  if (length(.w) == 1L) {
    .b <- ui$iniDf$backTransform
    if (!is.na(.b)) {
      return(.b)
    }
  }
  # If this has extra information in the mu-ref
  # (ie exp(tv + eta.v + 2)), then this shouldn't
  # have a default back-transformation
  .w <- which(ui$muRefExtra$parameter == par)
  if (length(.w) == 1L) {
    return("")
  }
  # Covariates should be reported without back-transformation
  .w <- which(ui$muRefCovariateDataFrame$covariateParameter == par)
  if (length(.w) == 1L) {
    return("")
  }
}
#' Get the parameter label and apply to parameter dataset
#'
#' @param ret Final fit environment
#' @param ui UI fit information
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.updateParFixedAddParameterLabel <- function(ret, ui) {
  .lab <- paste(ui$iniDf$label[!is.na(ui$iniDf$ntheta)])
  .lab[.lab == "NA"] <- ""
  .lab <- gsub(" *$", "", gsub("^ *", "", .lab))
  if (!all(.lab == "")) {
    ret$popDf <- data.frame(Parameter = .lab, ret$popDf, check.names = FALSE)
    ret$popDfSig <- data.frame(Parameter = .lab, ret$popDfSig, check.names = FALSE)
  }
}
#'  This applies the manually specified back-transformations
#'
#' @param .ret focei environment
#' @return Nothing, called for side effecs
#' @author Matthew L. Fidler
#' @examples
.updateParFixedApplyManualBacktransformations <- function(.ret) {
  .qn <- qnorm(1.0-(1-0.95)/2)
  .btName <- names(.ret$popDfSig)[4]
  .sigdig <- rxode2::rxGetControl(.ui, "sigdig", 3L)
  .fmt <- paste0("%", .sigdig, "g (%", .sigdig, "g, %", .sigdig, "g)")
  .fmt2 <- paste0("%", .sigdig, "g")
  lapply(seq_along(.ret$popDf$Estimate), function(i) {
    theta <- row.names(.ret$popDf)[i]
    .w <- which(ui$iniDf$name == theta)
    if (length(.w) == 1L) {
      .b <- ui$iniDf$backTransform[.w]
      if (!is.na(.b)) {
        .est <- .ret$popDf$Estimate[i]
        .se <- .ret$popDf$SE[i]
        .bt <- .ret$popDf[["Back-transformed"]]
        .bt[i] <- get(.b, envir=globalenv())(.est)
        .ret$popDf[["Back-transformed"]] <- .bt
        if (is.na(.se)) {
          .i1 <- .ret$popDf[["CI Lower"]]
          .low <- get(.b, envir=globalenv())(.est - .se * .qn)
          .i1[i] <- .low
          .ret$popDf[["CI Lower"]] <- .i1
          .i1 <- .ret$popDf[["CI Upper"]]
          .hi <- get(.b, envir=globalenv())(.est + .se * .qn)
          .i1[i] <- .hi
          .ret$popDf[["CI Upper"]] <- .i1
          .bt2 <- .ret$popDfSig[[.btName]]
          .bt2[i] <- sprintf(.fmt, .est, .low, .hi)
          .ret$popDfSig[[.btName]] <- .bt2
        } else {
          .bt2 <- .ret$popDfSig[[.btName]]
          .bt2[i] <- sprintf(.fmt, .est, .low, .hi)
          .ret$popDfSig[[.btName]] <- .bt2
          .bt2 <- .ret$popDfSig[[.btName]]
          .bt2[i] <- sprintf(.fmt2, .est)
          .ret$popDfSig[[.btName]] <- .bt2
        }
      }
    }
  })
}

#' This gest the CV/SD for a single ETA
#'
#'
#' @param .eta Eta Name
#' @param .env Environment where the indicators of `.sdOnly`, `.cvOnly` are stored so the column name can be changed to match the data
#' @param .ome Omega fixed vector
#' @param .muRefCurEval The current mu ref evaluation.  This determines if the ETA is logit normal and %CV should be calculated.
#' @param .sigdig is the number of significant digits used in the evaulation
#' @return Data frame row with ch= the charaacter representation and v is the vector representation of the CV or sd
#' @author Matthew L. Fidler
#' @noRd
.updateParFixedGetEtaRow <- function(.eta, .env, .ome, .omegaFix, .muRefCurEval, .sigdig) {
  .v <- .ome[.eta, .eta]
  .w <- which(.muRefCurEval$parameter == .eta)
  if (.muRefCurEval$curEval[.w] == "exp") {
    assign(".sdOnly", FALSE, envir=.env)
    return(data.frame(
      ch = paste0(
        ifelse(.omegaFix[.y], "fix(", ""),
        formatC(signif(sqrt(exp(.v) - 1) * 100, digits = .sigdig),
                digits = .digs, format = "fg", flag = "#"),
        ifelse(.omegaFix[.y], ")", "")
      ),
      v = sqrt(exp(.v) - 1) * 100))
  } else {
    assign(".cvOnly", FALSE, envir=.env)
    return(data.frame(
      ch = paste0(
        ifelse(.omegaFix[.y], "fix(", ""),
        formatC(signif(sqrt(.v), digits = .sigdig),
                digits = .digs, format = "fg", flag = "#"),
        ifelse(.omegaFix[.y], ")", "")),
      v = .v))
  }
}

#'  This will add the between subject varaibility to the mu-referenced theta.  It also expands the table to include non-mu referenced ETAs
#'
#'
#' @param .ret The focei return environment
#' @param .ui The rxode2 ui environment
#' @return Nothing called for side effects on popDf and popDfSig in the .ret environment
#' @author Matthew L. Fidler
#' @noRd
.updateParFixedAddBsv <- function(.ret, .ui) {
  .ome <- .ret$omega
  .omegaFix <- as.data.frame(.ui$ini)
  .omegaFix <- .omegaFix[is.na(.omegaFix$ntheta), ]
  .omegaFix <- setNames(.omegaFix$fix, paste(.omegaFix$name))
  .sigdig <- rxode2::rxGetControl(.ui, "sigdig", 3L)

  .muRefDataFrame <- .ui$muRefDataFrame
  .muRefCurEval   <- .ui$muRefCurEval
  .env <- new.env(parent=emptyenv())
  .env$.cvOnly <- TRUE
  .env$.sdOnly <- TRUE
  .env$.muRefVars <- NULL
  .cvp <- lapply(row.names(.ret$popDfSig), function(x) {
    .w <- which(.muRefDataFrame$theta == x)
    if (length(.w) != 1) {
      return(data.frame(ch = " ", v = NA_real_))
    }
    .eta <- .muRefDataFrame$eta[.w]
    assign(".muRefVars", c(.env$.muRefVars, .eta), envir=.env)
    .updateParFixedGetEtaRow(.eta, .env, .ome, .omegaFix, .muRefCurEval, .sigdig)
  })
  .cvp <- do.call("rbind", .cvp)
  .nonMuRef <- setdiff(dimnames(.ome)[[1]], .env$.muRefVars)
  if (length(.nonMuRef) > 0) {
    .ret$popDfSig2 <- as.data.frame(lapply(names(.ret$popDfSig), function(x) { rep("", length(.noMuRef))}))
    .ret$popDf2 <- as.data.frame(lapply(names(.ret$popDf), function(x) { rep(NA_real_, length(.noMuRef))}))
    row.names(.ret$popDfSig2) <- .nonMuRef
    row.names(.ret$popDf2) <- .nonMuRef
  }
  .ret$popDfSig <- data.frame(.ret$popDfSig, "BSD" = .cvp$ch, check.names = FALSE)
  .ret$popDf <- data.frame(.ret$popDf, "BSD" = .cvp$v, check.names = FALSE)
  if (length(.nonMuRef) > 0) {
    .cvp <- lapply(row.names(.ret$popDfSig2), function(x) {
      .updateParFixedGetEtaRow(x, .env, .ome, .omegaFix, .muRefCurEval, .sigdig)
    })
    .cvp <- do.call("rbind", .cvp)
    .ret$popDfSig2 <- data.frame(.ret$popDfSig2, "BSD" = .cvp$ch, check.names = FALSE)
    .ret$popDf2 <- data.frame(.ret$popDf2, "BSD" = .cvp$v, check.names = FALSE)
    .ret$popDfSig <- rbind(.ret$popDfSig, .ret$popDfSig2)
    .ret$popDf <- rbind(.ret$popDf, .ret$popDf2)
    rm(c("popDf2", "popDfSig2"), envir=.ret)
  }
  .w <- which(names(.ret$popDfSig) == "BSD")
  if (length(.w) == 1) {
    names(.ret$popDfSig)[.w] <- ifelse(.env$.sdOnly, "BSV(SD)", ifelse(.env$.cvOnly, "BSV(CV%)", "BSV(CV% or SD)"))
  }
  .w <- which(names(.ret$popDf) == "BSD")
  if (length(.w) == 1) {
    names(.ret$popDf)[.w] <- ifelse(.env$.sdOnly, "BSV(SD)", ifelse(.env$.cvOnly, "BSV(CV%)", "BSV(CV% or SD)"))
  }
}

.updateParFixedAddShrinkage <- function(.ret, .ui) {
  .shrink <- .ret$shrink
  .errs <- as.data.frame(.ui$ini)
  .errs <- paste(.errs[which(!is.na(.errs$err)), "name"])
  .muRefDataFrame <- .ui$muRefDataFrame
  .sh <- lapply(row.names(.ret$popDfSig), function(x) {
    .w <- which(.muRefDataFrame$theta == x)
    if (length(.w) != 1) {
      .w <- which(names(.shrink) == x)
      if (length(.w) != 1) {
        if (any(x == .errs)) {
          x <- "IWRES"
        } else {
          return(data.frame(ch = " ", v = NA_real_))
        }
      }
      .eta <- x
    } else {
      .eta <- .muRefDataFrame$eta[.w]
    }
    .v <- .shrink[7, .y]
    if (length(.v) != 1) {
      return(data.frame(ch = " ", v = NA_real_))
    }
    if (is.na(.v)) {
      return(data.frame(ch = " ", v = NA_real_))
    }
    .t <- ">"
    if (.v < 0) {
    } else if (.v < 20) {
      .t <- "<"
    } else if (.v < 30) {
      .t <- "="
    }
    return(data.frame(
      ch = sprintf("%s%%%s", formatC(signif(.v, digits = .digs),
                                     digits = .digs, format = "fg", flag = "#"
                                     ), .t),
      v = .v
    ))
  })
  .sh <- do.call("rbind", .sh)
  .ret$popDfSig <- data.frame(.ret$popDfSig, "Shrink(SD)%" = .sh$ch, check.names = FALSE)
  .ret$popDf <- data.frame(.ret$popDf, "Shrink(SD)%" = .sh$v, check.names = FALSE)
}
#' Create the parFixed dataset
#'
#' @param .ret focei style environment
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.updateParFixed <- function(.ret) {
  .ui <- .ret$ui
  .updateParFixedApplyManualBacktransformations(.ret)
  .updateParFixedAddParameterLabel(.ret, .ui)
  .updateParFixedAddBsv(.ret, .ui)
  .updateParFixedAddShrinkage(.ret, .ui)
  .ret$parFixed <- .ret$popDfSig
  .ret$parFixedDf <- .ret$popDf
  class(.ret$parFixed) <- c("nlmixr2ParFixed", "data.frame")
}
