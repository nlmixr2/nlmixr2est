.getSimModel <- function(obj, hideIpred=FALSE) {
  .lines <- rxode2::rxCombineErrorLines(obj$ui)
  .f <- function(x) {
    if (is.atomic(x) || is.name(x) || is.pairlist(x)) {
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`<-`)) ||
            identical(x[[1]], quote(`=`))) {
        if (identical(x[[2]], quote(`ipredSim`))) {
          x[[2]] <- quote(`ipred`)
          if (hideIpred) {
            x[[1]] <- quote(`~`)
          } else {
            x[[1]] <- quote(`<-`)
          }
        } else if (identical(x[[2]], quote(`sim`))) {
          x[[2]] <- quote(`sim`)
          x[[1]] <- quote(`<-`)
        } else {
          x[[1]] <- quote(`~`)
        }
      }
      return(as.call(lapply(x, .f)))
    }
  }
  .f(.lines)
}


## Add rxode2 THETA/ETA replacement mini DSL
.repSim <- function(x, theta = c(), eta = c(), lhs = c()) {
  ret <- eval(parse(text = sprintf("quote({%s})", x)))
  f <- function(x) {
    if (is.atomic(x)) {
      return(x)
    } else if (is.name(x)) {
      return(x)
    } else if (is.pairlist(x)) {
      return(x)
    } else if (is.call(x)) {
      if (identical(x[[1]], quote(`[`))) {
        type <- tolower(as.character(x[[2]]))
        if (type == "theta") {
          return(eval(parse(text = sprintf("quote(%s)", theta[as.numeric(x[[3]])]))))
        } else if (type == "eta") {
          return(eval(parse(text = sprintf("quote(%s)", eta[as.numeric(x[[3]])]))))
        }
        stop("Only theta/eta translation supported.")
      } else if (length(x[[2]]) == 1 &&
        ((identical(x[[1]], quote(`=`))) ||
          (identical(x[[1]], quote(`~`))))) {
        if (any(as.character(x[[2]]) == lhs)) {
          if (any(as.character(x[[2]]) == c("rx_pred_", "rx_r_"))) {
            x[[1]] <- quote(`~`)
            return(as.call(lapply(x, f)))
          } else {
            return(NULL)
          }
        }
        return(as.call(lapply(x, f)))
      } else {
        return(as.call(lapply(x, f)))
      }
    } else {
      stop("Don't know how to handle type ", typeof(x),
        call. = FALSE
      )
    }
  }
  ret <- deparse(f(ret))[-1]
  ret <- ret[regexpr("^ *NULL$", ret) == -1]
  ret <- paste(ret[-length(ret)], collapse = "\n")
  return(rxode2::rxNorm(rxode2::rxGetModel(ret)))
}

.simInfo <- function(object) {
  .mod <- .getSimModel(object, hideIpred=FALSE)
  .omega <- object$omega
  .etaN <- dimnames(.omega)[[1]]
  .params <- nlme::fixed.effects(object)
  .params <- .params
  .dfObs <- object$nobs
  .nlmixr2Data <- nlme::getData(object)
  .dfSub <- object$nsub
  .env <- object$env
  if (exists("cov", .env)) {
    .thetaMat <- nlme::getVarCov(object)
  } else {
    .thetaMat <- NULL
  }
  if (all(is.na(object$ui$ini$neta1))) {
    .omega <- NULL
    .dfSub <- 0
  }
  .sigma <- object$ui$simulationSigma
  return(list(
    rx = .mod, params = .params, events = .nlmixr2Data,
    thetaMat = .thetaMat, omega = .omega, sigma = .sigma, dfObs = .dfObs, dfSub = .dfSub
  ))
}


#' Simulate a nlmixr2 solved system
#'
#' This takes the uncertainty in the model parameter estimates and to
#' simulate a number of theoretical studies.  Each study simulates a
#' realization of the parameters from the uncertainty in the fixed
#' parameter estimates.  In addition the omega and sigma matrices are
#' simulated from the uncertainty in the Omega/Sigma matrices based
#' on the number of subjects and observations the model was based on.
#'
#' @param object nlmixr2 object
#' @param ... Other arguments sent to \code{rxSolve}
#' @return A rxode2 solved object
#' @inheritParams rxode2::rxSolve
#' @export
nlmixr2Sim <- function(object, ...) {
  rxode2::rxSolveFree()
  rxode2::.setWarnIdSort(FALSE)
  on.exit({
    rxode2::.setWarnIdSort(TRUE)
  })
  save <- getOption("nlmixr2.save", FALSE)
  .si <- .simInfo(object)
  .xtra <- list(...)
  if (any(names(.xtra) == "rx")) {
    .si$rx <- .xtra$rx
  }
  if (!is.null(.xtra$modelName)) {
    message(sprintf("Compiling %s model...", .xtra$modelName), appendLF = FALSE)
  } else {
    message("Compiling model...", appendLF = FALSE)
  }
  .newobj <- eval(.si$rx)
  on.exit({
    rxode2::rxUnload(.newobj)
  })
  message("done")
  if ((any(names(.xtra) == "nStud") && .xtra$nStud <= 1) || !any(names(.xtra) == "nStud")) {
    .si$thetaMat <- NULL
    .si$dfSub <- NULL
    .si$dfObs <- NULL
  } else {
    if (rxode2::rxIs(.xtra$thetaMat, "matrix")) {
      .si$thetaMat <- .xtra$thetaMat
    } else if (!is.null(.xtra$thetaMat)) {
      if (is.na(.xtra$thetaMat)) {
        .si$thetaMat <- NULL
      }
    }
    if (any(names(.xtra) == "dfSub")) {
      .si$dfSub <- .xtra$dfSub
    }
    if (any(names(.xtra) == "dfObs")) {
      .si$dfObs <- .xtra$dfObs
    }
  }


  if (any(names(.xtra) == "omega")) {
    .si$omega <- .xtra$omega
    if (any(is.na(.xtra$omega))) {
      .si$omega <- NULL
    }
  }
  if (any(names(.xtra) == "sigma")) {
    .si$sigma <- .xtra$sigma
    if (any(is.na(.xtra$sigma))) {
      .si$sigma <- NULL
    }
  }
  if (any(names(.xtra) == "events") &&
    rxode2::rxIs(.xtra$events, "rx.event")) {
    .si$events <- .xtra$events
  }
  if (any(names(.xtra) == "params")) {
    .si$params <- .xtra$params
  }
  .xtra$object <- .newobj
  .xtra$params <- .si$params
  .xtra$events <- .si$events
  if (rxode2::rxIs(.xtra$thetaMat, "matrix")) {
    .xtra$thetaMat <- NULL
  } else {
    .xtra$thetaMat <- .si$thetaMat
  }
  .xtra$dfObs <- .si$dfObs
  .xtra$omega <- .si$omega
  .xtra$dfSub <- .si$dfSub
  .xtra$sigma <- .si$sigma
  if (save) {
    .modName <- ifelse(is.null(object$uif$model.name), "", paste0(object$uif$model.name, "-"))
    if (.modName == ".-") .modName <- ""
    .dataName <- ifelse(is.null(object$uif$data.name), "", paste0(object$uif$data.name, "-"))
    if (.dataName == ".-") .dataName <- ""
    .digest <- digest::digest(list(
      gsub("<-", "=", gsub(" +", "", object$uif$fun.txt)),
      as.data.frame(object$uif$ini),
      .xtra,
      as.character(utils::packageVersion("nlmixr2")),
      as.character(utils::packageVersion("rxode2"))
    ))
    .saveFile <- file.path(
      getOption("nlmixr2.save.dir", getwd()),
      paste0("nlmixr2-nlmixr2Sim-", .modName, .dataName, "-", .digest, ".rds")
    )
    if (file.exists(.saveFile)) {
      message(sprintf("Loading nlmixr2Sim already run (%s)", .saveFile))
      .ret <- readRDS(.saveFile)
      return(.ret)
    }
  }
  .ret <- do.call(getFromNamespace("rxSolve", "rxode2"), .xtra, envir = parent.frame(2))
  if (inherits(.ret, "rxSolve")) {
    .rxEnv <- attr(class(.ret), ".rxode2.env")
    if (!is.null(.xtra$nsim)) {
      .rxEnv$nSub <- .xtra$nsim
    }
    if (!is.null(.xtra$nSub)) {
      .rxEnv$nSub <- .xtra$nSub
    }
    if (is.null(.xtra$nStud)) {
      .rxEnv$nStud <- 1
    } else {
      .rxEnv$nStud <- .xtra$nStud
    }
    .cls <- c("nlmixr2Sim", class(.ret))
    attr(.cls, ".rxode2.env") <- .rxEnv
    if (any(names(.ret) == "CMT") && any(names(object) == "CMT")) {
      if (is(object$CMT, "factor")) {
        .ret$CMT <- as.integer(.ret$CMT)
        levels(.ret$CMT) <- levels(object$CMT)
        class(.ret$CMT) <- "factor"
      }
    }
    class(.ret) <- .cls
  }
  if (save) {
    saveRDS(.ret, file = .saveFile)
  }
  return(.ret)
}

#' @export
plot.nlmixr2Sim <- function(x, y, ...) {
  p1 <- eff <- Percentile <- sim.id <- id <- p2 <- p50 <- p05 <- p95 <- . <- NULL
  .args <- list(...)
  save <- getOption("nlmixr2.save", FALSE)
  rxode2::rxReq("dplyr")
  rxode2::rxReq("tidyr")
  if (is.null(.args$p)) {
    .p <- c(0.05, 0.5, 0.95)
  } else {
    .p <- .args$p
  }
  if (save) {
    .digest <- digest::digest(list(
      .args,
      as.character(utils::packageVersion("nlmixr2")),
      as.character(utils::packageVersion("rxode2"))
    ))
    .saveFile <- file.path(
      getOption("nlmixr2.save.dir", getwd()),
      paste0("nlmixr2SimPlot-", .digest, ".rds")
    )
    if (file.exists(.saveFile)) {
      message(sprintf("Loading nlmixr2SimPlot already summarized (%s)", .saveFile))
      .ret <- readRDS(.saveFile)
      return(.ret)
    }
  }
  if (x$env$nStud <= 1) {
    if (x$env$nSub < 2500) {
      warning("In order to put confidence bands around the intervals, you need at least 2500 simulations.")
      message("Summarizing data for plot")
      .ret <- x %>%
        dplyr::group_by(time) %>%
        dplyr::do(data.frame(p1 = .p, eff = quantile(.$sim, probs = .p))) %>%
        dplyr::mutate(Percentile = factor(sprintf("%02d%%", round(p1 * 100))))
      message("done.")
      .ret <- ggplot2::ggplot(.ret, aes(time, eff, col = Percentile, fill = Percentile)) +
        ggplot2::geom_line(size = 1.2)
      return(.ret)
    } else {
      .n <- round(sqrt(x$env$nSub))
    }
  } else {
    .n <- x$env$nStud
  }
  message("Summarizing data for plot")
  .ret <- x %>%
    dplyr::mutate(id = sim.id %% .n) %>%
    dplyr::group_by(id, time) %>%
    dplyr::do(data.frame(p1 = .p, eff = quantile(.$sim, probs = .p))) %>%
    dplyr::group_by(p1, time) %>%
    dplyr::do(data.frame(p2 = .p, eff = quantile(.$eff, probs = .p))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(p2 = sprintf("p%02d", (p2 * 100))) %>%
    tidyr::spread(p2, eff) %>%
    dplyr::mutate(Percentile = factor(sprintf("%02d%%", round(p1 * 100))))
  message("done.")
  .ret <- ggplot2::ggplot(.ret, aes(time, p50, col = Percentile, fill = Percentile)) +
    ggplot2::geom_ribbon(aes(ymin = p05, ymax = p95), alpha = 0.5) +
    ggplot2::geom_line(size = 1.2)
  if (save) {
    saveRDS(.ret, file = .saveFile)
  }
  return(.ret)
}

.predOnlyRxDsl <- function(x) {
  if (is.atomic(x)) {
    return(x)
  } else if (is.name(x)) {
    return(x)
  } else if (is.pairlist(x)) {
    return(x)
  } else if (is.call(x)) {
    if (length(x) >= 2) {
      if (length(x[[2]]) == 1 &&
            ((identical(x[[1]], quote(`=`))) ||
               (identical(x[[1]], quote(`~`))))) {
        x[[1]] <- quote(`~`)
        return(as.call(lapply(x, .predOnlyRxDsl)))
      } else {
        if (length(x) >= 2) {
          if (((identical(x[[1]], quote(`=`))) ||
                 (identical(x[[1]], quote(`~`)))) &&
                length(x[[2]] == 3)) {
            if (identical(x[[2]][[1]], quote(`/`))) {
              x[[1]] <- quote(`~`)
              return(as.call(lapply(x, .predOnlyRxDsl)))
            }
          }
        }
      }
    }
    return(as.call(lapply(x, .predOnlyRxDsl)))
  } else {
    stop("Don't know how to handle type ", typeof(x),
         call. = FALSE
         )
  }
}

## Mini DSL to fix pred-only models
.predOnlyRx <- function(object) {
  .ret <- eval(parse(text = sprintf("quote({%s})", rxode2::rxNorm(object$model$pred.only))))
  .ret <- deparse(.predOnlyRxDsl(.ret))[-1]
  .ret <- .ret[regexpr("^ *NULL$", .ret) == -1]
  .ret <- .ret[-length(.ret)]
  .w <- rev(which(regexpr("^ *cmt[(].*[)] *$", .ret) != -1))
  .cur <- 1
  ## Remove cmt(#) at the beginning
  while (any(.w == .cur)) {
    .w <- .w[.w != .cur]
    .cur <- .cur + 1
  }
  ## Put rx_pred_ at the end, but before cmt(name) statements
  if (length(.w) > 2) {
    while (.w[1] == .w[2] + 1) {
      .w <- .w[-1]
      if (length(.w) <= 2) break
    }
  }
  if (length(.w) > 0) {
    .w <- .w[1]
    .ret[.w] <- paste0("pred=rx_pred_\n", .ret[.w])
  } else {
    .ret <- c(.ret, "pred=rx_pred_")
  }
  .ret <- paste(.ret, collapse = "\n")
  return(rxode2::rxode2(.ret))
}

#' Predict a nlmixr2 solved system
#'
#' @param ipred Flag to calculate individual predictions. When
#'     \code{ipred} is \code{TRUE}, calculate individual predictions.
#'     When \code{ipred} is \code{FALSE}, set calculate typical population predictions.
#'     When \code{ipred} is \code{NA}, calculate both individual and
#'     population predictions.
#'
#' @inheritParams rxode2::rxSolve
#'
#' @return an rxode2 solved data frame with the predictions
#'
#' @export
#'
#' @export
nlmixr2Pred <- function(object, ..., ipred = FALSE) {
  rxode2::.setWarnIdSort(FALSE)
  on.exit(rxode2::.setWarnIdSort(TRUE))
  lst <- as.list(match.call()[-1])
  if (rxode2::rxIs(lst$params, "rx.event")) {
    if (!is.null(lst$events)) {
      tmp <- lst$events
      lst$events <- lst$params
      lst$params <- tmp
    } else {
      lst$events <- lst$params
      lst$params <- NULL
    }
  }
  if (!rxode2::rxIs(lst$events, "rx.event")) {
    lst$events <- nlmixr2Data(getData(object))
  }
  message("Compiling model...", appendLF = FALSE)
  lst$object <- .predOnlyRx(object)
  message("done")
  params <- fixed.effects(object)
  names(params) <- sprintf("THETA[%d]", seq_along(params))
  do.ipred <- FALSE
  do.pred <- FALSE
  if (is.na(ipred)) {
    do.ipred <- TRUE
    do.pred <- TRUE
  } else if (ipred) {
    do.ipred <- TRUE
  } else {
    do.pred <- TRUE
  }
  if (do.ipred) {
    re <- random.effects(object)
    if (is.null(re)) {
      .tmp <- lst$events
      .w <- which(tolower(names(.tmp)) == "id")
      .nid <- length(unique(.tmp[[.w]]))
      ipred.par <- data.frame(t(params),
        rx_err_ = rep(0, .nid),
        check.names = FALSE
      )
    } else {
      re <- re[, -1]
      names(re) <- sprintf("ETA[%d]", seq_along(names(re)))
      ipred.par <- data.frame(re, t(params),
        rx_err_ = 0, check.names = FALSE
      )
    }
  }
  if (do.pred) {
    neta <- dim(object$omega)[1]
    if (neta == 0) {
      pred.par <- c(params, rx_err_ = 0)
    } else {
      pred.par <- c(params, setNames(rep(0, neta + 1), c(sprintf("ETA[%d]", seq(1, neta)), "rx_err_")))
    }
  }
  on.exit(
    {
      rxode2::rxUnload(lst$object)
    },
    add = TRUE
  )
  if (!is.na(ipred)) {
    if (do.pred) {
      lst$params <- pred.par
    } else {
      lst$params <- ipred.par
    }
    ret <- suppressWarnings(do.call(getFromNamespace("rxSolve", "rxode2"), lst, envir = parent.frame(2)))
    if (do.ipred) {
      names(ret) <- sub("pred", "ipred", names(ret))
    }
    return(ret)
  } else {
    lst$params <- pred.par
    ret.pred <- suppressWarnings(do.call(getFromNamespace("rxSolve", "rxode2"), lst, envir = parent.frame(2)))
    lst$params <- ipred.par
    ret.pred$ipred <- suppressWarnings(do.call(getFromNamespace("rxSolve", "rxode2"), lst, envir = parent.frame(2)))$pred
    return(ret.pred)
  }
}
#' @rdname nlmixr2Pred
#' @export
predict.nlmixr2FitData <- function(object, ...) {
  nlmixr2Pred(object, ...)
}

#' @rdname nlmixr2Sim
#' @export
rxSolve.nlmixr2FitData <- function(object, params = NULL, events = NULL, inits = NULL,
                                  scale = NULL, method = c("liblsoda", "lsoda", "dop853", "indLin"),
                                  transitAbs = NULL, atol = 1.0e-8, rtol = 1.0e-6,
                                  maxsteps = 70000L, hmin = 0, hmax = NA_real_,
                                  hmaxSd = 0, hini = 0, maxordn = 12L, maxords = 5L, ...,
                                  cores,
                                  covsInterpolation = c("locf", "linear", "nocb", "midpoint"),
                                  addCov = FALSE, matrix = FALSE, sigma = NULL, sigmaDf = NULL,
                                  sigmaLower = -Inf, sigmaUpper = Inf,
                                  nCoresRV = 1L, sigmaIsChol = FALSE,
                                  sigmaSeparation = c("auto", "lkj", "separation"),
                                  sigmaXform = c("identity", "variance", "log", "nlmixr2Sqrt", "nlmixr2Log", "nlmixr2Identity"),
                                  nDisplayProgress = 10000L,
                                  amountUnits = NA_character_, timeUnits = "hours", stiff,
                                  theta = NULL,
                                  thetaLower = -Inf, thetaUpper = Inf,
                                  eta = NULL, addDosing = FALSE,
                                  stateTrim = Inf, updateObject = FALSE,
                                  omega = NULL, omegaDf = NULL, omegaIsChol = FALSE,
                                  omegaSeparation = c("auto", "lkj", "separation"),
                                  omegaXform = c("variance", "identity", "log", "nlmixr2Sqrt", "nlmixr2Log", "nlmixr2Identity"),
                                  omegaLower = -Inf, omegaUpper = Inf,
                                  nSub = 1L, thetaMat = NULL, thetaDf = NULL, thetaIsChol = FALSE,
                                  nStud = 1L, dfSub = 0.0, dfObs = 0.0, returnType = c("rxSolve", "matrix", "data.frame", "data.frame.TBS", "data.table", "tbl", "tibble"),
                                  seed = NULL, nsim = NULL,
                                  minSS = 10L, maxSS = 1000L,
                                  infSSstep = 12,
                                  strictSS = TRUE,
                                  istateReset = TRUE,
                                  subsetNonmem = TRUE,
                                  maxAtolRtolFactor = 0.1,
                                  from = NULL,
                                  to = NULL,
                                  by = NULL,
                                  length.out = NULL,
                                  iCov = NULL,
                                  keep = NULL,
                                  indLinPhiTol = 1e-7,
                                  indLinPhiM = 0L,
                                  indLinMatExpType = c("expokit", "Al-Mohy", "arma"),
                                  indLinMatExpOrder = 6L,
                                  drop = NULL,
                                  idFactor = TRUE,
                                  mxhnil = 0,
                                  hmxi = 0.0,
                                  warnIdSort = TRUE,
                                  warnDrop = TRUE,
                                  ssAtol = 1.0e-8,
                                  ssRtol = 1.0e-6,
                                  safeZero = TRUE,
                                  sumType = c("pairwise", "fsum", "kahan", "neumaier", "c"),
                                  prodType = c("long double", "double", "logify"),
                                  sensType = c("advan", "autodiff", "forward", "central"),
                                  linDiff=c(tlag=1.5e-5, f=1.5e-5, rate=1.5e-5, dur=1.5e-5, tlag2=1.5e-5, f2=1.5e-5, rate2=1.5e-5, dur2=1.5e-5),
                                  linDiffCentral=c(tlag=TRUE, f=TRUE, rate=TRUE, dur=TRUE, tlag2=TRUE, f2=TRUE, rate2=TRUE, dur2=TRUE),
                                  resample=NULL,
                                  resampleID=TRUE) {
  do.call("nlmixr2Sim", as.list(match.call()[-1]), envir = parent.frame(2))
}

#' @rdname nlmixr2Sim
#' @export
simulate.nlmixr2FitData <- function(object, nsim = 1, seed = NULL, ...) {
  nlmixr2::nlmixr2Sim(object, ..., nsim = nsim, seed = seed)
}

#' @rdname nlmixr2Sim
#' @export
solve.nlmixr2FitData <- function(a, b, ...) {
  lst <- as.list(match.call()[-1])
  n <- names(lst)
  if (!missing(a)) {
    n[n == "a"] <- ""
  }
  if (!missing(b)) {
    n[n == "b"] <- ""
  }
  names(lst) <- n
  do.call("nlmixr2Sim", lst, envir = parent.frame(2))
}

#' @importFrom rxode2 rxSolve
#' @export
rxode2::rxSolve
