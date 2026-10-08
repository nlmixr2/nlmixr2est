# Choosing a covariance from what the covariance step computed.  The C++ FOCEi step
# computes the R and S matrices a request needs and calls .covSelectFocei() to choose.

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

#' Choose the FOCEi covariance from the R and S matrices the C++ step computed
#'
#' A request that needs R falls to S when R is not usable; "r,s" falls to R when S is
#' not usable and is checked by `.covSandwichChoice()`; a covariance whose variances
#' are all below 1e-7 is none.  Installs `env$cov` and `env$covMethod` (the label, e.g.
#' "r,s", "r+", "|s|") and gives the warnings.
#' @param env fit environment (`covR`, `covS`, `covRS` and the S inverse `.covSinv`)
#' @param req requested slot: 1 "r,s", 2 "r", 3 "s"
#' @param rState,sState each matrix: 0 not computed, 1 usable, 2 not positive definite,
#'   3 its computation failed
#' @param rstr,sstr labels of R and S
#' @param checkSandwich whether a repair or pseudo-inverse was needed
#' @param sHasZero whether a subject's S score had to be substituted
#' @param covSmall `foceiControl(covSmall=)`
#' @return list(slot = 0 (none), 1, 2 or 3, label)
#' @noRd
.covSelectFocei <- function(env, req, rState, sState, rstr, sstr, checkSandwich, sHasZero, covSmall) {
  .cur <- req
  .which <- ""
  if (req %in% c(1L, 2L) && rState != 1L) {
    .cur <- 3L
  }
  .orig <- .cur
  if (.cur == 2L) {
    .which <- "covR"
  }
  if (.cur %in% c(1L, 3L)) {
    if (sState == 2L) {
      if (.cur == 1L) {
        .which <- "covR"
        .cur <- 2L
      } else {
        warning("cannot calculate covariance", call. = FALSE)
        .cur <- 0L
      }
    } else if (sState == 1L) {
      if (.cur == 1L) {
        .which <- .covSandwichChoice(env$covRS, env$covR, env$covS, rstr, sstr, checkSandwich, covSmall)
        .cur <- c(covRS = 1L, covR = 2L, covS = 3L)[[.which]]
      } else {
        .which <- ".covSinv"
      }
    } else if (sState == 3L) {
      if (.cur == 1L) {
        cat("\rS matrix calculation failed; Switch to R-matrix covariance.\n")
        .which <- "covR"
        .cur <- 2L
      } else {
        cat("\rCould not calculate covariance matrix.\n")
        warning("cannot calculate covariance", call. = FALSE)
        .cur <- 0L
      }
    }
  }
  if (.cur != 0L && nzchar(.which)) {
    .cov <- get(.which, envir = env, inherits = FALSE)
    if (all(diag(.cov) < 1e-7)) {
      warning("The variance of all elements are unreasonably small, <1e-7", call. = FALSE)
      .cur <- 0L
      if (exists("cov", envir = env, inherits = FALSE)) rm(list = "cov", envir = env)
    } else {
      assign("cov", .cov, envir = env)
    }
  }
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
