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
rxUiGet.saemResLlMod <- function(x, ...) {
  .iniErr <- rxUiGet.saemErrDf(x, ...)
  if (is.null(.iniErr)) return(numeric(0))
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
    rxUiGet.saemErrMuNames(x, ...))
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
