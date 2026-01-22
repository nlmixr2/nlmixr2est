.nlmixr2iov <- function(val, type, transform) {
  # First get the standard deviation
  if (transform == "logvar") {
    sd <- sqrt(exp(val))
  } else if (transform == "logsd") {
    sd <- exp(val)
  } else if (transform == "sd") {
    sd <- abs(val)
  } else if (transform == "var") {
    sd <- sqrt(abs(val))
  } else {
    stop("Unknown transform")
  }
  if (type == "exp") {
    100 * sqrt(exp(sd^2) - 1)
  } else {
    sd
  }
}

#' Transform the estimated value to \%CV for IOV
#'
#' @param val estimated value
#' @return IOV value
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
nlmixr2iovLogvarCv <- function(val) {
  .nlmixr2iov(val, "exp", "logvar")
}
#' @rdname nlmixr2iovLogvarCv
#' @export
nlmixr2iovLogvarSd <- function(val) {
  .nlmixr2iov(val, "", "logvar")
}
#' @rdname nlmixr2iovLogvarCv
#' @export
nlmixr2iovLogsdCv <- function(val) {
  .nlmixr2iov(val, "exp", "logsd")
}

#' @rdname nlmixr2iovLogvarCv
#' @export
nlmixr2iovLogsdSd <- function(val) {
  .nlmixr2iov(val, "", "logsd")
}

#' @rdname nlmixr2iovLogvarCv
#' @export
nlmixr2iovSdCv <- function(val) {
  .nlmixr2iov(val, "exp", "sd")
}

#' @rdname nlmixr2iovLogvarCv
#' @export
nlmixr2iovSdSd <- function(val) {
  .nlmixr2iov(val, "", "sd")
}

#' @rdname nlmixr2iovLogvarCv
#' @export
nlmixr2iovVarCv <- function(val) {
  .nlmixr2iov(val, "exp", "var")
}

#' @rdname nlmixr2iovLogvarCv
#' @export
nlmixr2iovVarSd <- function(val) {
  .nlmixr2iov(val, "", "var")
}
# This stores information about the IOV model that can be used
# in nlmixr2 fits
.uiIovEnv <- new.env(parent = emptyenv())
.uiIovEnv$iovVars <- NULL
#' This applies the IOV method to the model based on the data used
#'
#' @param env environment to apply the IOV model transformation.  This should contain:
#'
#' - `ui`: the model to apply the IOV transformation to
#' - `data`: the data to use for the IOV transformation
#' - `control`: the control object, which should contain `iovXform`
#' @return nothing, called for side effects
#' @noRd
#' @author Matthew L. Fidler
.uiApplyIov <- function(env) {
  .uiIovEnv$iovVars <- NULL
  .xform <- env$control$iovXform
  if (length(.xform)  != 1) {
    .xform <- "sd"
  }
  if (!(.xform %in% c("sd", "var", "logsd", "logvar"))) {
    .xform <- "sd"
  }
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
                       .lst <- c(lapply(.var, function(v) {
                         # Add theta to dataset; represents variance of iov
                         .curTheta <- .theta1
                         # This is in variance and needs to be converted
                         # based on the xform
                         .est <- .iniDf[which(.iniDf$name == v &
                                                is.na(.iniDf$ntheta)), "est"]
                         if (.xform == "var") {
                           .curTheta$est <- .est
                         } else if (.xform == "sd") {
                           .curTheta$est <- sqrt(.est)
                         } else if (.xform == "logvar") {
                           .curTheta$est <- log(.est)
                         } else if (.xform == "logsd") {
                           .curTheta$est <- log(sqrt(.est))
                         }
                         .curTheta$name <- v
                         .uiIovEnv$iovVars <- c(.uiIovEnv$iovVars, v)
                         .curTheta$fix <- .fixed

                         .w <- which(env$ui$muRefCurEval$parameter == v)
                         if (length(.w) == 1L) {
                           .curEval <- env$ui$muRefCurEval$curEval[.w]
                         } else {
                           .curEval <- ""
                         }
                         .curTheta$backTransform <-
                           paste0(switch(.xform,
                                  "sd" = "nlmixr2iovSd",
                                  "var" = "nlmixr2iovVar",
                                  "logsd" = "nlmixr2iovLogsd",
                                  "logvar" = "nlmixr2iovLogvar"),
                                  ifelse(.curEval=="exp", "Cv", "Sd"))
                         if (.xform %in% c("sd", "var")) {
                           .curTheta$lower <- 0 # doesn't work with saem
                         }
                         .env$maxtheta <- .curTheta$ntheta <- .env$maxtheta + 1L
                         .env$thetas <- rbind(.env$thetas, .curTheta)
                         for (n in .lvls[[l1]]) {
                           .curEta <- .eta1
                           .curEta$name <- paste0("rx.", v, ".", n)
                           .curEta$label <- paste0(v, "(", l1, "==", n, ")")
                           .env$drop <- c(.env$drop, .curEta$name)
                           .env$maxeta <- .curEta$neta1 <-
                             .curEta$neta2 <- .env$maxeta + 1L
                           .env$etas <- rbind(.env$etas, .curEta)
                         }
                         if (.xform == "logsd") {
                           str2lang(paste0("rx.", v, " <- exp(", v, ")*(",
                                           paste(paste0("rx.", v, ".", .lvls[[l1]],
                                                        "*(", l1,
                                                        " == ", .lvls[[l1]], ")"),
                                                 collapse="+"),
                                           ")"))
                         } else if (.xform == "logvar") {
                             str2lang(paste0("rx.", v, " <- sqrt(exp(", v, "))*(",
                                             paste(paste0("rx.", v, ".", .lvls[[l1]],
                                                          "*(", l1,
                                                          " == ", .lvls[[l1]], ")"),
                                                   collapse="+"),
                                             ")"))
                         } else if (.xform == "sd") {
                           str2lang(paste0("rx.", v, " <- abs(", v, ")*(",
                                           paste(paste0("rx.", v, ".", .lvls[[l1]],
                                                        "*(", l1,
                                                        " == ", .lvls[[l1]], ")"),
                                                 collapse="+"),
                                           ")"))
                         } else if (.xform == "var") {
                           str2lang(paste0("rx.", v, " <- sqrt(abs(", v, "))*(",
                                           paste(paste0("rx.", v, ".", .lvls[[l1]],
                                                        "*(", l1,
                                                        " == ", .lvls[[l1]], ")"),
                                                 collapse="+"),
                                           ")"))
                         }
                       }),
                       lapply(.var, function(v) {
                         str2lang(paste0(v, ".rx <- rx.", v))
                       }))
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
#' @param ret data frame with some iov information dropped
#' @return fit with iov information dropped
#' @noRd
#' @author Matthew L. Fidler
.uiFinalizeIov <- function(ret) {
  if (!is.null(.uiIovEnv$iov)) {
    if (is.null(ret$ui)) return(ret)

    if (is.environment(ret$env)) {
      .iniDf <- .uiIovEnv$iov$iniDf
      .finalDf <- ret$ui$iniDf
      .iovName <- new.env(parent=emptyenv())
      .iovName$var <- character(0)
      .est <- vapply(seq_along(.iniDf$name), function(i) {
        if (is.na(.iniDf$neta1[i])) {
          .w <- which(.finalDf$name == .iniDf$name[i])
          .finalDf[.w, "est"]
        } else if (.iniDf$condition[i] == "id") {
          .w <- which(.finalDf$name == .iniDf$name[i])
          .finalDf[.w, "est"]
        } else {
          .w <- which(.finalDf$name == .iniDf$name[i])
          .iovName$var <- c(.iovName$var, .iniDf$name[i])
          .fun <- sub("Cv$", "Sd", .finalDf[.w, "backTransform"])
          .fun <- get(.fun)
          .fun(.finalDf[.w, "est"])^2
        }
      }, double(1), USE.NAMES = FALSE)
      names(.est) <- .iniDf$name

      # Now we can update the finalDf
      assign("iniDf0", .iniDf, envir = ret$env)
      .finalDf <- .iniDf
      .finalDf$est <- .est
      .ui <- .uiIovEnv$iov
      suppressMessages(rxode2::ini(.ui) <- .finalDf)
      assign("ui", .ui, envir = ret$env)

      # Adjust Matrices to remove dummy IOV components
      .omega <- ret$env$omega
      .d1 <- dimnames(.omega)[[1]]
      .d1 <- .d1[!(.d1 %in% .uiIovEnv$iovDrop)]
      assign("omega", .ui$omega, envir = ret$env)

      .omega <- .ui$omega
      .n <- names(.omega)
      .n <- .n[.n != "id"]
      .omega <- lapply(.n, function(x) {
        .omega[[x]]
      })
      names(.omega) <- .n

      .nid <- length(ret$env$eta$ID)

      # Fix shrinkage, now split out by iov variable
      .shrink <- ret$env$shrink
      .w <- which(names(.shrink) %in% .uiIovEnv$iovDrop)
      .shrink0 <- .shrink[,-.w]
      .shrinkN <- cbind(data.frame(type=row.names(.shrink0)),
                        .shrink)
      .shrink1 <- lapply(.n, function(var) {
        .cur <- .omega[[var]]
        .dn <- dimnames(.cur)[[1]]
        do.call(`cbind`,
                lapply(.dn, function(d) {
                  .w <- c(which(grepl(d, names(.shrinkN), fixed=TRUE)))
                  .curd <- .shrinkN[,.w]

                  # n = the same for each eta (#id)
                  # mean = sum(mean)/n_means
                  # var = sum(var)/n_vars
                  # sd = sqrt(var)
                  # skewness = sum(skewness)/n_skewness
                  # kurtosis = sum(kurtosis+3)/n_kurtosis - 3
                  .nv <- length(.curd["var",])
                  .var <- sum(.curd["var", ])/.nv

                  .mean <- sum(.curd["mean",])/.nv
                  .sd <- sqrt(.var)

                  .tv <- .mean*sqrt(.nid*.nv)/.sd

                  .curi <- data.frame(v=c(.mean,
                                          .var,
                                          .sd,
                                          sum(.curd["skewness",])/.nv,
                                          sum(.curd["kurtosis",] + 3)/.nv - 3,
                                          (1-.var/.est[d])*100, # var shrinkage
                                          (1 - sqrt(.var)/sqrt(.est[d]))*100, # sd shrinkage
                                          .tv,
                                          2*pt(-abs(.tv), df = (.nid*.nv) -1)
                                          ),
                                      row.names=row.names(.curd))
                  names(.curi) <- d
                  .curd <- cbind(.curd, .curi)
                }))
      })
      names(.shrink1) <- .n
      .shrink <- c(list(id=.shrink0),
                   .shrink1)
      assign("shrink", .shrink, envir = ret$env)


      ## .phiC <- ret$env$phiC
      ## .phiC <- lapply(seq_along(.phiC),
      ##                 function(i) {
      ##                   .m <- .phiC[[i]]
      ##                   .m <- .m[.d1, .d1, drop=FALSE]
      ##                   .m
      ##                 })
      ## names(.phiC) <- names(ret$env$phiC)
      ## assign("phiC", .phiC, envir = ret$env)


      ## .phiH <- ret$env$phiH
      ## .phiH <- lapply(seq_along(.phiH),
      ##                 function(i) {
      ##                   .m <- .phiH[[i]]
      ##                   .m <- .m[.d1, .d1, drop=FALSE]
      ##                   .m
      ##                 })
      ## names(.phiH) <- names(ret$env$phiH)

      ## assign("phiH", .phiH, envir = ret$env)


      # Fix eta objective function; Maybe save the full one for
      # passing the etaMat information to the next estimation method

      # Now fix the random effect matrix
      .ranef <- ret$env$ranef

      .w <- which(names(.ranef) %in% .uiIovEnv$iovDrop)
      .iov <- .ranef
      .ranef <- .ranef[,-.w]
      assign("ranef", .ranef, envir = ret$env)

      .w <- which(names(.iov) %in% c(.uiIovEnv$iovDrop, "ID"))
      .iov <- .iov[,.w]

      .sdIov <- sqrt(.est)

      .dt <- NULL
      .iov <- lapply(.n, function(var) {
        .cur <- .omega[[var]]
        .dn <- dimnames(.cur)[[1]]
        for (d in .dn) {
          .w <- c(1L,which(grepl(d, names(.iov), fixed=TRUE)))
          .curd <- data.table::data.table(.iov[,.w])
          .curd <- data.table::melt(.curd,
                                    id.vars=names(.curd)[1],
                                    measure.vars=names(.curd)[-1],
                                    variable.name = var,
                                    value.name = d)
          # Since this is scaled by the standard deviation, we can
          # calculate it from the derived eta fixed to 1 by mutiplying
          # by the standard deviation of the IOV variable (calculated above)
          .curd[[d]] <- .curd[[d]] *.sdIov[d]
          if (is.null(.dt)) {
            .dt <- .curd
          } else {
            .dt[[d]] <- .curd[[d]]
          }
        }
        .dt[[var]] <- as.integer(sub(paste0("rx.", d, "."), "", as.character(.dt[[var]]), fixed=TRUE))
        .dt <- as.data.frame(.dt)
        .dt <- .dt[order(.dt[[1]], .dt[[2]]), , drop=FALSE]
        rownames(.dt) <- NULL
        .dt
      })
      names(.iov) <- .n
      assign("iov", .iov, envir = ret$env)

      # Now fixed effects
      .fixef <- ret$env$fixef
      .w <- which(names(.fixef) %in% .iovName$var)
      .fixef <- .fixef[-.w]
      assign("fixef",.fixef, envir = ret$env)

      .parFixedDf <- ret$env$parFixedDf
      .bck <- which(grepl("Back",names(.parFixedDf)))
      .bsv <- which(grepl("BSV", names(.parFixedDf)))
      .est <- which(grepl("Est", names(.parFixedDf)))

      .valCharPrep <-
        .parFixedDf[.uiIovEnv$iovVars,.bsv] <-
        .parFixedDf[.uiIovEnv$iovVars, .bck]
      .parFixedDf[.uiIovEnv$iovVars,.bsv] <- NA_real_
      .parFixedDf[.uiIovEnv$iovVars,.est] <- NA_real_

      .parFixedDf <- .parFixedDf[!grepl("^rx[.]", rownames(.parFixedDf)),]
      assign("parFixedDf", .parFixedDf, envir = ret$env)

      .parFixed <- ret$env$parFixed
      .bck2 <- which(grepl("Back", names(.parFixed)))
      .bsv2 <- which(grepl("BSV",  names(.parFixed)))
      .est2 <- which(grepl("Est", names(.parFixed)))

      .sigdig <- ret$control$sigdig
      .parFixed[.uiIovEnv$iovVars, .bck2] <- ""
      .parFixed[.uiIovEnv$iovVars, .est2] <- ""
      .parFixed[.uiIovEnv$iovVars, .bsv2] <- formatC(
        signif(.valCharPrep, digits = .sigdig),
        digits = .sigdig, format = "fg", flag = "#")
      .parFixed <- .parFixed[!grepl("^rx[.]", rownames(.parFixed)),]
      assign("parFixed", .parFixed, envir=ret$env)
    }
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
      .rename <- paste0(.uiIovEnv$iovVars, ".rx")
      names(ret) <- vapply(names(ret), function(n) {
        if (n %in% .rename) {
          sub("[.]rx$", "", n)
        } else {
          n
        }
      }, character(1), USE.NAMES = FALSE)
    }
  }
  ret
}
