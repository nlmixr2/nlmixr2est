#' @export
rxUiGet.saemErrDf <- function(x, ...) {
  .x <- x[[1]]
  .w <- !is.na(.x$iniDf$ntheta) & !is.na(.x$iniDf$err)
  if (length(.w) == 0L) return(NULL)
  .x$iniDf[.w, ]
}

#' @export
rxUiGet.saemErrMuNames <- function(x, ...) {
  .iniErr <- rxUiGet.saemErrDf(x, ...)
  if (is.null(.iniErr)) return(NULL)
  paste0("rx_err_", seq_along(.iniErr$name), "_")
}

#' @export
rxUiGet.saemErrMuEst <- function(x, ...) {
  .iniErr <- rxUiGet.saemErrDf(x, ...)
  if (is.null(.iniErr)) return(NULL)
  setNames(
    vapply(seq_along(.iniErr$est),
           function(i) {
             .low <- .iniErr$lower[i]
             .hi <- .iniErr$upper[i]
             .est <- .iniErr$est[i]
             if (!is.finite(.low) && !is.finite(.hi)) {
               return(logit(.est, .low, .hi))
             } else if (is.finite(.low) && !is.finite(.hi)) {
               return(log(.est - .low))
             } else if (!is.finite(.low) && is.finite(.hi)) {
               return(log(.low - .est))
             } else if (!is.finite(.low) && !is.finite(.hi)) {
               return(.est)
             }
           }, numeric(1), USE.NAMES=FALSE),
    rxUiGet.saemErrMuNames(x, ...)
    )
}

#' @export
rxUiGet.saemErrMuEstLog <- function(x, ...) {
  .iniErr <- rxUiGet.saemErrDf(x, ...)
  if (is.null(.iniErr)) return(NULL)
  setNames(
    vapply(seq_along(.iniErr$est),
           function(i) {
             .low <- .iniErr$lower[i]
             .hi <- .iniErr$upper[i]
             if (is.finite(.low) && !is.finite(.hi) && .low == 0) {
               return(TRUE)
             }
             FALSE
           }, logical(1), USE.NAMES=FALSE),
  rxUiGet.saemErrMuNames(x, ...))
}

#' @export
rxUiGet.saemErrMuEst <- function(x, ...) {
  .iniErr <- rxUiGet.saemErrDf(x, ...)
  if (is.null(.iniErr)) return(NULL)
  setNames(
    vapply(seq_along(.iniErr$est),
           function(i) {
             .low <- .iniErr$lower[i]
             .hi <- .iniErr$upper[i]
             .est <- .iniErr$est[i]
             if (!is.finite(.low) && !is.finite(.hi)) {
               return(logit(.est, .low, .hi))
             } else if (is.finite(.low) && !is.finite(.hi)) {
               return(log(.est - .low))
             } else if (!is.finite(.low) && is.finite(.hi)) {
               return(log(.low - .est))
             } else if (!is.finite(.low) && !is.finite(.hi)) {
               return(.est)
             }
           }, numeric(1), USE.NAMES=FALSE),
    rxUiGet.saemErrMuNames(x, ...)
    )
}

#' @export
rxUiGet.saemErrMuFix <- function(x, ...) {
  .iniErr <- rxUiGet.saemErrDf(x, ...)
  if (is.null(.iniErr)) return(NULL)
  setNames(
    vapply(seq_along(.iniErr$est),
           function(i) {
             .iniErr$fix[i]
           }, logical(1), USE.NAMES=FALSE),
    rxUiGet.saemErrMuNames(x, ...)
    )
}


#' @export
rxUiGet.saemErrMu <- function(x, ...) {
  .iniErr <- rxUiGet.saemErrDf(x, ...)
  if (is.null(.iniErr)) return(NULL)
  lapply(seq_along(.iniErr$est),
         function(i) {
           .name <- .iniErr$name[i]
           .low <- .iniErr$lower[i]
           .hi <- .iniErr$upper[i]
           if (!is.finite(.low) && !is.finite(.hi)) {
             return(str2lang(paste0(.name, "<- expit(rx_err_", i,
                                    "_, ", .low, ", ", .hi, ")")))
           } else if (is.finite(.low) && !is.finite(.hi)) {
             return(str2lang(paste0(.name, "<- ", .low, " + exp(rx_err_", i, "_)")))
           } else if (!is.finite(.low) && is.finite(.hi)) {
             return(str2lang(paste0(.name, "<- ", .hi, " - exp(rx_err_", i, "_)")))
           } else if (!is.finite(.low) && !is.finite(.hi)) {
             return(str2lang(paste0(.name, "<- rx_err_", i, "_)")))
           }
         })
}

#' @export
rxUiGet.saemModelNeedsLlik <- function(x, ...) {
  .ui <- x[[1]]
  .ret <- any(.ui$predDf$distribution != "norm")
  if (.ret) return(.ret)
  # FIXME: check for complicated error expressions
  FALSE
}

.saemIntegrateErrMu <- function(ui, lines) {
  if (ui$saemModelNeedsLlik) {
    .ref <- ui$saemErrMu
    .ret <- lines
    .ret[[2]] <- as.call(c(list(quote(`{`)),
                           .ref,
                           lapply(seq_along(lines[[2]])[-1], function(i) {
                             lines[[2]][[i]]
                           })))
    return(.ret)
  }
  lines
}

#' @export
rxUiGet.saemDistribution <- function(x, ...) {
  if (rxUiGet.saemModelNeedsLlik(x, ...)) return(2L)
  1L
}

.saemParamsToEstimateLlik <- function(pars, x, ...) {
  return(pars)
  if (!rxUiGet.saemModelNeedsLlik(x, ...)) return(pars)
  c(pars, rxUiGet.saemErrMuNames(x, ...))
}

.saemFixedLlik <- function(dft, x, ...) {
  return(dft)
  if (!rxUiGet.saemModelNeedsLlik(x, ...)) return(dft)
  .iniErr <- rxUiGet.saemErrDf(x, ...)
  c(dft, setNames(.iniErr$fix, rxUiGet.saemErrMuNames(x, ...)))
}

.saemInitThetaLlik <- function(est, x, ...) {
  return(est)
  if (!rxUiGet.saemModelNeedsLlik(x, ...)) return(est)
  c(est, rxUiGet.saemErrMuEst(x, ...))
}

#' @export
rxUiGet.saemIniDf <- function(x, ...) {
  .iniDf <- get("iniDf", x[[1]])
  if (!rxUiGet.saemModelNeedsLlik(x, ...)) return(.iniDf)
  .theta <- .iniDf[!is.na(.iniDf$ntheta), ]
  .eta <- .iniDf[is.na(.iniDf$ntheta), ]
  .maxTheta <- max(.theta$ntheta, na.rm=TRUE)
  .names <- rxUiGet.saemErrMuNames(x, ...)
  .est <- rxUiGet.saemErrMuEst(x, ...)
  .fix <- rxUiGet.saemErrMuFix(x, ...)
  do.call("rbind",
          c(list(.theta),
            lapply(seq_along(.names), function(i) {
              data.frame(ntheta=.maxTheta+i, neta1=NA_real_, neta2=NA_real_,
                         name=.names[i], lower=-Inf, est=.est[i], upper=Inf,
                         fix=.fix[i], label=NA_character_, backTransform=NA_character_,
                         condition=NA_character_, err=NA_character_)
            }),
            list(.eta)))
}

#' @export
rxUiGet.saemMuRefCurEval <- function(x, ...) {
  .muRefCurEval <- get("muRefCurEval", x[[1]])
  if (!rxUiGet.saemModelNeedsLlik(x, ...)) return(.muRefCurEval)
  .iniErr <- rxUiGet.saemErrDf(x, ...)
  do.call("rbind",
         c(list(.muRefCurEval),
           lapply(seq_along(.iniErr$est),
                  function(i) {
                    .nameTheta <- paste0("rx_err_", i,"_")
                    .low <- .iniErr$lower[i]
                    .hi <- .iniErr$upper[i]
                    if (!is.finite(.low) && !is.finite(.hi)) {
                      return(data.frame(parameter=.nameTheta, curEval="expit",
                                        low=.low, hi=.hi))
                    } else if (is.finite(.low) && !is.finite(.hi)) {
                      return(data.frame(parameter=.nameTheta, curEval="exp",
                                        low=NA_real_, hi=NA_real_))
                    } else if (!is.finite(.low) && is.finite(.hi)) {
                      return(data.frame(parameter=.nameTheta, curEval="exp",
                                        low=NA_real_, hi=NA_real_))
                    } else if (!is.finite(.low) && !is.finite(.hi)) {
                      return(NULL)
                    }
                    NULL
                  })))
}


