# Choosing a covariance from what the covariance step computed.  The C++ FOCEi step
# computes the R and S matrices a request needs and calls .covSelectFocei() to choose.

# The FOCEi covariance methods that can fall back, and what they can fall back to
.covFallbackMethods <- c("r,s", "r", "s", "analytic")
.covFallbackTargets <- c("r,s", "r", "s")

#' A control's covariance fallbacks
#'
#' A control saved before `covFallback=` existed has none recorded; it gets the default
#' of the control function, the fallbacks it was computed with.
#' @param control a control list
#' @param type control function name, e.g. "foceiControl" or "saemControl"
#' @return named list
#' @noRd
.covFallbackOf <- function(control, type = "foceiControl") {
  .fb <- control[["covFallback"]]
  if (is.list(.fb)) {
    return(.fb)
  }
  eval(formals(get(type, envir = asNamespace("nlmixr2est")))$covFallback)
}

#' Check a control's `covFallback=`
#'
#' @param covFallback named list, one element per covariance method, each the
#'   ordered methods it falls back to
#' @param methods the methods that can have fallbacks
#' @param targets the methods that can be fallbacks
#' @return the list, each element a character vector
#' @noRd
.covFallbackCheck <- function(covFallback, methods = .covFallbackMethods, targets = .covFallbackTargets) {
  if (is.null(covFallback)) {
    return(list())
  }
  if (!is.list(covFallback) || (length(covFallback) > 0L && is.null(names(covFallback)))) {
    stop("'covFallback' must be a named list, e.g. list(\"r,s\" = c(\"r\", \"s\"))", call. = FALSE)
  }
  .n <- names(covFallback)
  if (any(!nzchar(.n)) || anyDuplicated(.n)) {
    stop("'covFallback' needs one uniquely named element per covariance method", call. = FALSE)
  }
  .bad <- setdiff(.n, methods)
  if (length(.bad) > 0L) {
    stop(
      sprintf(
        "'covFallback' names a method without fallbacks: %s (allowed: %s)",
        paste(dQuote(.bad, FALSE), collapse = ", "),
        paste(dQuote(methods, FALSE), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  for (.m in .n) {
    .v <- covFallback[[.m]]
    if (is.null(.v)) {
      .v <- character(0)
    }
    if (!is.character(.v) || anyNA(.v) || anyDuplicated(.v)) {
      stop(sprintf("'covFallback$%s' must be distinct method names", .m), call. = FALSE)
    }
    .bad <- setdiff(.v, setdiff(targets, .m))
    if (length(.bad) > 0L) {
      stop(
        sprintf(
          "'covFallback$%s' cannot fall back to %s (allowed: %s)",
          .m,
          paste(dQuote(.bad, FALSE), collapse = ", "),
          paste(dQuote(setdiff(targets, .m), FALSE), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    covFallback[[.m]] <- .v
  }
  covFallback
}

#' Which of the sandwich, R and S covariances the "r,s" request installs
#'
#' The sandwich, unless an R or S repair, a pseudo-inverse (`checkSandwich`) or a
#' sandwich variance below `covSmall` makes it doubtful; then the one not repaired, or
#' the R, S and sandwich diagonals are compared on the scale the rule was calibrated
#' to (`covR * 2`, `covS * 4`, before issue 666 rescaled them).
#' @param covRS,covR,covS the three covariances
#' @param rstr,sstr labels of R and S ("r", "r+", "|r|", ...)
#' @param checkSandwich whether a repair or pseudo-inverse was needed
#' @param covSmall `foceiControl(covSmall=)`
#' @return "covRS", "covR" or "covS"
#' @noRd
.covSandwichChoice <- function(covRS, covR, covS, rstr, sstr, checkSandwich, covSmall) {
  .rsSmall <- any(abs(diag(covRS)) < covSmall)
  .check2 <- !checkSandwich && .rsSmall
  if (!(checkSandwich || .check2)) {
    return("covRS")
  }
  if (!.check2 && rstr == "r") {
    return("covR")
  }
  if (!.check2 && sstr == "s") {
    return("covS")
  }
  .rsd <- sum(diag(covRS))
  .rSmall <- any(abs(2 * diag(covR)) < covSmall)
  .rd <- sum(2 * diag(covR))
  .sSmall <- any(abs(4 * diag(covS)) < covSmall)
  .sd <- sum(4 * diag(covS))
  if (.rsSmall && .sSmall && .rSmall) {
    "covRS"
  } else if (.rsSmall && .sSmall && !.rSmall) {
    "covR"
  } else if (.rsSmall && !.sSmall && .rSmall) {
    "covS"
  } else if (.rsd > .rd) {
    if (.rd > .sd) "covS" else "covR"
  } else if (.rsd > .sd) {
    "covS"
  } else {
    "covRS"
  }
}

#' What the state of an R or S matrix says about it (see `.covSelectFocei()`)
#' @param st 0 not computed, 1 usable, 2 not positive definite, 3 not computable
#' @return description
#' @noRd
.covMatState <- function(st) {
  c("not computed", "usable", "not positive definite", "could not be computed")[st + 1L]
}

#' Add a method and its outcome to the record of a covariance choice
#' @param acc environment holding `tried`, the list of rows
#' @param method covariance method
#' @param outcome "used" or why it was not
#' @return invisibly `acc`
#' @noRd
.covTriedAdd <- function(acc, method, outcome) {
  acc$tried[[length(acc$tried) + 1L]] <- data.frame(method = method, outcome = outcome)
  invisible(acc)
}

#' Choose the FOCEi covariance from the R and S matrices the C++ step computed
#'
#' The request first; when its matrices are not usable, the methods `fallback` lists,
#' in order (a request that needs R falls to "s" when R is not usable, "r,s" to "r" when
#' S is not); a doubtful sandwich is checked by `.covSandwichChoice()`, which may only
#' pick a listed method; a covariance whose variances are all below 1e-7 is none.
#' Installs `env$cov`, `env$covMethod` (the label, e.g. "r,s", "r+", "|s|") and
#' `env$covTried` (each method tried and why it was not used) and gives the warnings.
#' @param env fit environment (`covR`, `covS`, `covRS` and the S inverse `.covSinv`)
#' @param req requested slot: 1 "r,s", 2 "r", 3 "s"
#' @param rState,sState each matrix: 0 not computed, 1 usable, 2 not positive definite,
#'   3 its computation failed
#' @param rstr,sstr labels of R and S
#' @param checkSandwich whether a repair or pseudo-inverse was needed
#' @param sHasZero whether a subject's S score had to be substituted
#' @param covSmall `foceiControl(covSmall=)`
#' @param fallback methods the request may fall back to (`foceiControl(covFallback=)`)
#' @return list(slot = 0 (none), 1, 2 or 3, label)
#' @noRd
.covSelectFocei <- function(env, req, rState, sState, rstr, sstr, checkSandwich, sHasZero, covSmall,
                            fallback = c("r", "s")) {
  .names <- c("r,s", "r", "s")
  .acc <- new.env(parent = emptyenv())
  .acc$tried <- list()
  .cur <- req
  .which <- ""
  if (req %in% c(1L, 2L) && rState != 1L) {
    .covTriedAdd(.acc, .names[req], paste("R", .covMatState(rState)))
    .cur <- if ("s" %in% fallback) 3L else 0L
  }
  .orig <- .cur
  if (.cur == 2L) {
    .which <- "covR"
  }
  if (.cur %in% c(1L, 3L)) {
    if (sState == 2L || sState == 3L) {
      .covTriedAdd(.acc, .names[.cur], paste("S", .covMatState(sState)))
      if (sState == 3L && .cur == 1L && "r" %in% fallback) {
        cat("\rS matrix calculation failed; Switch to R-matrix covariance.\n")
      }
      if (.cur == 1L && "r" %in% fallback) {
        .which <- "covR"
        .cur <- 2L
      } else {
        if (sState == 3L) {
          cat("\rCould not calculate covariance matrix.\n")
        }
        warning("cannot calculate covariance", call. = FALSE)
        .cur <- 0L
      }
    } else if (sState == 1L) {
      if (.cur == 1L) {
        .which <- .covSandwichChoice(env$covRS, env$covR, env$covS, rstr, sstr, checkSandwich, covSmall)
        if (.which != "covRS" && !(c(covR = "r", covS = "s")[[.which]] %in% fallback)) {
          .which <- "covRS"
        }
        if (.which != "covRS") {
          .covTriedAdd(.acc, "r,s", "sandwich not used (covSmall check)")
        }
        .cur <- c(covRS = 1L, covR = 2L, covS = 3L)[[.which]]
      } else {
        .which <- ".covSinv"
      }
    }
  }
  if (.cur != 0L && nzchar(.which)) {
    .cov <- get(.which, envir = env, inherits = FALSE)
    if (all(diag(.cov) < 1e-7)) {
      warning("The variance of all elements are unreasonably small, <1e-7", call. = FALSE)
      .covTriedAdd(.acc, .names[.cur], "all variances below 1e-7")
      .cur <- 0L
      if (exists("cov", envir = env, inherits = FALSE)) rm(list = "cov", envir = env)
    } else {
      assign("cov", .cov, envir = env)
      .covTriedAdd(.acc, .names[.cur], "used")
    }
  }
  assign("covTried", do.call(rbind, .acc$tried), envir = env)
  if (.cur == 0L) {
    warning("covariance step failed", call. = FALSE)
    return(list(slot = 0L, label = "failed"))
  }
  if (sHasZero) {
    warning("S matrix had problems solving for some subject and parameters", call. = FALSE)
  }
  .doWarn <- FALSE
  if (.cur != 3L) {
    if (rstr == "|r|") {
      warning("R matrix non-positive definite but corrected by R = sqrtm(R%*%R)", call. = FALSE)
      .doWarn <- TRUE
    } else if (rstr == "r+") {
      warning("R matrix non-positive definite but corrected (because of cholAccept)", call. = FALSE)
      .doWarn <- TRUE
    }
  }
  if (.cur == 1L) {
    if (sstr == "|s|") {
      warning("S matrix non-positive definite but corrected by S = sqrtm(S%*%S)", call. = FALSE)
      .doWarn <- TRUE
    } else if (sstr == "s+") {
      warning("S matrix non-positive definite but corrected (because of cholAccept)", call. = FALSE)
      .doWarn <- TRUE
    }
    if (.doWarn) {
      warning("since sandwich matrix is corrected, you may compare to $covR or $covS if you wish", call. = FALSE)
    }
    .label <- paste0(rstr, ",", sstr)
  } else if (.cur == 2L) {
    .label <- rstr
    if (.orig != 2L) {
      warning(
        if (checkSandwich) {
          "using R matrix to calculate covariance, can check sandwich or S matrix with $covRS and $covS"
        } else {
          "using R matrix to calculate covariance"
        },
        call. = FALSE
      )
    }
  } else {
    .label <- sstr
    if (.orig != 2L) {
      warning(
        if (checkSandwich) {
          "using S matrix to calculate covariance, can check sandwich or R matrix with $covRS and $covR"
        } else {
          "using S matrix to calculate covariance"
        },
        call. = FALSE
      )
    }
  }
  assign("covMethod", .label, envir = env)
  list(slot = .cur, label = .label)
}
