# This stores information about the IOV model that can be used
# in nlmixr2 fits
.uiIovEnv <- new.env(parent = emptyenv())
#' This applies the IOV method to the model based on the data used
#'
#' @param env environment to apply the IOV model transformation.  This should contain:
#'
#' - `ui`: the model to apply the IOV transformation to
#' - `data`: the data to use for the IOV transformation
#' @return nothing, called for side effects
#' @noRd
#' @author Matthew L. Fidler
.uiApplyIov <- function(env) {
  .ui <- env$ui
  .iniDf <- .ui$iniDf
  .lvls <- .iniDf$condition[which(!is.na(.iniDf$condition) &
                                    .iniDf$condition != "id" &
                                     is.na(.iniDf$err))]
  if (length(.lvls) > 0) {
    .n <- .iniDf[which(.iniDf$condition %in% .lvls), "name"]
    .ui <- suppressWarnings(eval(str2lang(paste0("rxode2::rxRename(.ui, ",
                                paste(paste0("rx.", .n, "=", .n),
                                      collapse=", "), ")"))))

    # For the new iniDf, we will take out all the level variables and
    # then renumber the etas
    .thetas <- .iniDf[is.na(.iniDf$neta1),, drop=FALSE]
    .etas <- .iniDf[is.na(.iniDf$ntheta),, drop=FALSE]
    if (length(.thetas$name) > 0) {
      .maxtheta <- max(.thetas$ntheta, na.rm = TRUE)
      .theta1 <- .thetas[1,]
      .theta1$ntheta <- .maxtheta
    } else {
      .maxtheta <- 0L
      .theta1 <- .etas[1,]
      .theta1$ntheta <- 0L
      .theta1$neta1 <- NA_real_
      .theta1$neta2 <- NA_real_
    }
    .theta1$label <- NA_character_
    .eta1 <- .etas[1, ]
    .eta1$fix <- TRUE
    .eta1$neta1 <- .eta1$neta2 <- 0
    .eta1$est <- 1

    .etas <- .etas[which(!(.etas$condition %in% .lvls)), , drop=FALSE]
    if (length(.etas$name) > 0) {
      .etas$neta1 <- factor(.etas$neta1, levels = sort(unique(.etas$neta1)))
      .etas$neta2 <- factor(.etas$neta2, levels = sort(unique(.etas$neta2)))
      .etas$neta1 <- as.integer(.etas$neta1)
      .etas$neta2 <- as.integer(.etas$neta2)
      .maxeta <- max(.etas$neta1, na.rm = TRUE)
    } else {
      .maxeta <- 0L
      .theta1 <- .etas[1,]
    }

    .data <- env$data
    .lvls <- setNames(lapply(.lvls, function(l) {
      .v <- sort(unique(.data[[l]]))
      if (is.null(.v)) {
        stop(paste0("IOV variable '", l, "' is not present in the data "),
             call. = FALSE)
      }
      if (!is.numeric(.v)) {
        stop(paste0("IOV variable '", l, "' must be numeric"),
             call. = FALSE)
      }
      .v
    }), .lvls)
    .env <- new.env(parent = emptyenv())
    .env$thetas <- .thetas
    .env$etas <- .etas
    .env$maxtheta <- .maxtheta
    .env$maxeta <- .maxeta
    .env$drop <- NULL
    # Now we have enough information to create the IOV variables
    # changed to etas on id
    .lines <- lapply(names(.lvls),
                     function(l1) {
                       .w <-which(.iniDf$condition == l1)
                       .var <- .iniDf$name[.w]
                       .fixed <- .iniDf$fix[.w]
                       .lst <- lapply(.var, function(v) {
                         # Add theta to dataset; represents variance of iov
                         .curTheta <- .theta1
                         .curTheta$est <- .iniDf[which(.iniDf$name == v &
                                                         is.na(.iniDf$ntheta)), "est"]
                         .curTheta$name <- v
                         .curTheta$fix <- .fixed
                         .env$maxtheta <- .curTheta$ntheta <- .env$maxtheta + 1L
                         .env$thetas <- rbind(.env$thetas, .curTheta)
                         for (n in .lvls[[l1]]) {
                           .curEta <- .eta1
                           .curEta$name <- paste0("rx.", v, ".", n)
                           .env$drop <- c(.env$drop, .curEta$name)
                           .env$maxeta <- .curEta$neta1 <-
                             .curEta$neta2 <- .env$maxeta + 1L
                           .env$etas <- rbind(.env$etas, .curEta)
                         }
                         str2lang(paste0("rx.", v, " <- sqrt(abs(", v, "))*(",
                                         paste(paste0("rx.", v, ".", .lvls[[l1]],
                                                      "*(", l1,
                                                      " == ", .lvls[[l1]], ")"),
                                               collapse="+"),
                                         ")"))
                       })
                       .lst
                     })
    .lines <- do.call(`c`, c(.lines, list(.ui$lstExpr)))
    .ui <- rxode2::rxUiDecompress(.ui)
    # Now the lines can be added to the model
    assign("iniDf", rbind(.env$thetas,.env$etas), envir = .ui)
    assign("lstExpr", .lines, envir = .ui)
    .uiIovEnv$iov <- env$ui
    .uiIovEnv$iovDrop <- .env$drop # extra variables to drop
    env$ui <- rxode2::rxUiDecompress(suppressWarnings(suppressMessages(.ui$fun())))
  } else {
    .uiIovEnv$iov <- NULL
    .uiIovEnv$iovDrop <- NULL
  }
}
#' Finalizes IOV model
#'
#'
#' @param ret data frame with some iov information dropped
#' @return fit with iov information dropped
#' @noRd
#' @author Matthew L. Fidler
.uiFinalizeIov <- function(ret) {
  if (!is.null(.uiIovEnv$iov)) {
    if (is.null(ret$ui)) return(ret)
    # In this approach the model is simply kept,
    # but the data drops the iovDrop
    if (inherits(ret, "data.frame")) {
      .w <- which(names(ret) %in% .uiIovEnv$iovDrop)
      if (length(.w) > 0L) {
        .cls <- class(ret)
        class(ret) <- "data.frame"
        ret <- ret[,-.w]
        class(ret) <- .cls
      }
    }
  }
  ret
}
