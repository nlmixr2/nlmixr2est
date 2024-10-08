#' Get the zero etas from the model
#'
#' @param ui rxode2 ui
#' @return Names of zero estimated etasi
#' @author Matthew L. Fidler
#' @noRd
.getZeroEtasFromModel <- function(ui) {
  .iniDf <- ui$iniDf[is.na(ui$iniDf$ntheta),, drop = FALSE]
  if (length(.iniDf$neta1) == 0) return(character(0))
  .r <- range(.iniDf$neta1)
  .r <- seq(.r[1], .r[2])
  .etaNames <- dimnames(ui$omega)[[1]]
  .zeroEta <- vapply(.r, function(i) {
    all(.iniDf[(.iniDf$neta1 == i) | (.iniDf$neta2 == i), "est"] == 0)
  }, logical(1), USE.NAMES=FALSE)
  .etaNames[.zeroEta]
}
#' Add back interesting mu etas, replace remaining zero etas with 0
#'
#' @param x expression
#' @param muRefDataFrame rxode2 muRefDataFRame
#' @param zeroEtas rxode2 zero etas that will be dropped
#' @return expression with interesting mus re-inserted
#' @author Matthew L. Fidler
#' @noRd
.addBackInterestingMuEtas <- function(x, muRefDataFrame, zeroEtas) {
  if (is.call(x)) {
    return(as.call(lapply(x, .addBackInterestingMuEtas, muRefDataFrame=muRefDataFrame,
                          zeroEtas=zeroEtas)))
  } else if (is.name(x)) {
    .n <- as.character(x)
    if (.n %in% zeroEtas) {
      return(0)
    }
    .w <- which(muRefDataFrame$theta == .n)
    if (length(.w) == 1L) {
      .mu <- muRefDataFrame[.w, ]
      .eta <- .mu$eta
      if (.eta %in% zeroEtas) {
        return(x)
      } else {
        return(str2lang(paste0(.mu$theta, "+", .mu$eta)))
      }
    }
  }
  x
}


.getUiFunFromIniAndModel <- function(ui, ini, model) {
  .ls <- ls(ui$meta, all.names=TRUE)
  .ret <- vector("list", length(.ls) + 3)
  .ret[[1]] <- quote(`{`)
  for (.i in seq_along(.ls)) {
    .ret[[.i + 1]] <- eval(parse(text=paste("quote(", .ls[.i], "<-", deparse1(ui$meta[[.ls[.i]]]), ")")))
  }
  .len <- length(.ls)
  .ret[[.len + 2]] <- ini
  .ret[[.len + 3]] <- model
  .retf <- function(){}
  body(.retf) <- as.call(.ret)
  .retf
}

#' This downgrades the UI for any of the zero etas in the model
#'
#' @param ui  rxode2 User interface function
#' @param zeroEtas The names of the zero etas in the model
#' @return New rxode2 ui with the zero etas removed
#' @author Matthew L. Fidler
#' @noRd
.downgradeEtas <- function(ui, zeroEtas=character(0)) {
  .lst <- .saemDropMuRefFromModel(ui, noCovs=TRUE)
  .model <- str2lang(
    paste0("model({",
           paste(vapply(lapply(.lst, .addBackInterestingMuEtas,
                               muRefDataFrame=ui$muRefDataFrame,
                               zeroEtas=zeroEtas),
                        function(x) {
                          deparse1(x)
                        }, character(1), USE.NAMES=FALSE),
                 collapse="\n"),
           "})"))
  .iniDf <- ui$iniDf
  .etas <- .iniDf[.iniDf$name %in% zeroEtas, "neta1"]
  .w <- which(.iniDf$neta1 %in% .etas | .iniDf$neta2 %in% .etas)
  if (length(.w) > 0) {
    .iniDf <- .iniDf[-.w, ]
    .thetas <- .iniDf[!is.na(.iniDf$ntheta), ]
    .etas <- .iniDf[is.na(.iniDf$ntheta),, drop = FALSE]
    if (length(.etas$neta1) > 0) {
      .fct <- factor(c(.etas$neta1, .etas$neta2))
      .etas$neta1 <- as.integer(.fct[seq_along(.etas$neta1)])
      .fct <- .fct[-seq_along(.etas$neta1)]
      .etas$neta2 <- as.integer(.fct)
      .iniDf <- rbind(.thetas, .etas)
    } else {
      .iniDf <- .thetas
    }
  }
  .ini <- as.expression(lotri::as.lotri(.iniDf))
  .ini[[1]] <- quote(`ini`)
  .mod <- .getUiFunFromIniAndModel(ui, .ini, .model)
  .mod()
}


#' This preprocesses the UI with any needed modifications
#'
#'
#' This is used to apply the mu referencing bug fix for `rxode2`
#'
#' @param ui rxode2 UI
#' @return correct ui for say a fit (possibly a simulation)
#' @author Matthew L. Fidler
#' @noRd
.nlmixrPreprocessUi <- function(ui, control) {
  ui <- rxode2::assertRxUi(ui)
  .ret <- rxode2::rxUiDecompress(ui)
  .checkLiteralFix <- TRUE
  if (is.null(control)) {
  } else if (checkmate::testLogical(control$literalFix, any.missing=FALSE, len=1, null.ok=FALSE)) {
    .checkLiteralFix <- control$literalFix
  }
  if (.checkLiteralFix) {
    .ui <- try(rxode2::rxFixPop(ui, returnNull=TRUE))
    if (inherits(.ui, "try-error")) .ui <- NULL
    if (!is.null(.ui)) {
      .ret <- rxode2::rxUiDecompress(.ui)
      nlmixr2global$nlmixr2EstEnv$uiUnfix <- ui
    }
  }
  .zeroEtas <- .getZeroEtasFromModel(.ret)
  if (length(.zeroEtas) > 0) {
    nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi <- rxode2::rxUiDecompress(.ret)
    .minfo(paste0("the following etas are removed from the model since their initial estimates are zero: ",
           paste(.zeroEtas, collapse=", ")))
    .ret <- .downgradeEtas(ui, zeroEtas=.zeroEtas)
  }
  .ret
}
