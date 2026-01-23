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
#' @param zeroEtasLookup optional pre-computed lookup table for O(1) access
#' @return expression with interesting mus re-inserted
#' @author Matthew L. Fidler
#' @noRd
.addBackInterestingMuEtas <- function(x, muRefDataFrame, zeroEtas, zeroEtasLookup = NULL) {
  if (is.call(x)) {
    return(as.call(lapply(x, .addBackInterestingMuEtas, muRefDataFrame=muRefDataFrame,
                          zeroEtas=zeroEtas, zeroEtasLookup=zeroEtasLookup)))
  } else if (is.name(x)) {
    .n <- as.character(x)
    # Optimized: Use pre-computed lookup table for O(1) access if available
    if (!is.null(zeroEtasLookup)) {
      if (isTRUE(zeroEtasLookup[[.n]])) {
        return(0)
      }
    } else if (.n %in% zeroEtas) {
      return(0)
    }
    .w <- which(muRefDataFrame$theta == .n)
    if (length(.w) == 1L) {
      .mu <- muRefDataFrame[.w, ]
      .eta <- .mu$eta
      # Optimized: Use pre-computed lookup table for O(1) access if available
      if (!is.null(zeroEtasLookup)) {
        if (isTRUE(zeroEtasLookup[[.eta]])) {
          return(x)
        }
      } else if (.eta %in% zeroEtas) {
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
#'
#' @param zeroEtas The names of the zero etas in the model
#'
#' @return New rxode2 ui with the zero etas removed
#'
#' @author Matthew L. Fidler
#'
#' @keywords internal
#'
#' @export
#'
.downgradeEtas <- function(ui, zeroEtas=character(0)) {
  .lst <- .saemDropMuRefFromModel(ui, noCovs=TRUE)
  # Optimized: Pre-compute lookup table for O(1) access in recursive function
  .zeroEtasLookup <- if (length(zeroEtas) > 0) setNames(rep(TRUE, length(zeroEtas)), zeroEtas) else NULL
  .model <- str2lang(
    paste0("model({",
           paste(vapply(lapply(.lst, .addBackInterestingMuEtas,
                               muRefDataFrame=ui$muRefDataFrame,
                               zeroEtas=zeroEtas,
                               zeroEtasLookup=.zeroEtasLookup),
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
#' Remove an eta from the model
#'
#'
#' @param ui rxode2 user interface
#' @param eta eta to remove
#' @return ui model with eta removed
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' mod <- function ()  {
#'  description <- "One compartment PK model with linear clearance"
#'  ini({
#'    lka <- 0.45
#'    lcl <- 1
#'    lvc <- 3.45
#'     propSd <- c(0, 0.5)
#'     etaKa ~ 0.1
#'   })
#'  model({
#'    ka <- exp(lka + etaKa)
#'    cl <- exp(lcl)
#'    vc <- exp(lvc)
#'    Cc <- linCmt()
#'    Cc ~ prop(propSd)
#'  })
#' }
#'
#' mod |> rmEta("etaKa")
#'
#' # This can also remove more than one eta
#'
#' mod <- function ()  {
#'  description <- "One compartment PK model with linear clearance"
#'  ini({
#'    lka <- 0.45
#'    lcl <- 1
#'    lvc <- 3.45
#'    propSd <- c(0, 0.5)
#'    etaKa ~ 0.1
#'    etaCl ~ 0.2
#'    etaVc ~ 0.3
#'   })
#'  model({
#'    ka <- exp(lka + etaKa)
#'    cl <- exp(lcl + etaCl)
#'    vc <- exp(lvc + etaVc)
#'    Cc <- linCmt()
#'    Cc ~ prop(propSd)
#'  })
#' }
#'
#' mod |> rmEta(c("etaKa", "etaCl"))
#'
rmEta <- function(ui, eta) {
  ui <- rxode2::assertRxUi(ui, " for the 'rmEta()' function")
  .eta0 <- as.character(substitute(eta))
  .eta <- try(eta, silent=TRUE)
  if (inherits(.eta, "try-error")) {
    eta <- .eta0
  } else if (is.character(.eta)) {
    eta <- .eta
  }
  checkmate::assertCharacter(eta, any.missing=FALSE, min.len=1)
  for (e in eta)
    rxode2::assertExists(ui, e)
  .downgradeEtas(ui, eta)
}


#' Preprocess the zero omegas
#'
#' @param ui rxode2 ui model
#' @inheritParams nlmixr2
#' @return list with the ui (possibly modified)
#' @export
#' @author Matthew L. Fidler
.preProcessZeroOmega <- function(ui, est, data, control) {
  .ui <- ui
  .zeroEtas <- .getZeroEtasFromModel(.ui)
  if (length(.zeroEtas) > 0) {
    nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi <- rxode2::rxUiDecompress(.ui)
    .minfo(paste0("the following etas are removed from the model since their initial estimates are zero: ",
                  paste(.zeroEtas, collapse=", ")))
    .ui <- .downgradeEtas(.ui, zeroEtas=.zeroEtas)
  }
  list(ui=.ui)
}
preProcessHooksAdd(".preProcessZeroOmega", .preProcessZeroOmega)
