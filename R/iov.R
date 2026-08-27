#' Check if an estimation method supports iov
#'
#' Uses the \code{"iov"} attribute on the \code{nlmixr2Est.<method>} S3
#' method (so external packages can register IOV support without editing
#' this file). The attribute may be \code{TRUE}/\code{FALSE}, or a
#' \code{function(control)} returning logical for conditional support.
#'
#' @param est estimation routine name
#' @param control control object
#' @return boolean
#' @noRd
#' @author Matthew Fidler
.isIovMethod <- function(est, control = NULL) {
  .v <- as.character(utils::methods("nlmixr2Est"))
  .method <- paste0("nlmixr2Est.", est)
  if (.method %in% .v) {
    .iov <- attr(utils::getS3method("nlmixr2Est", est), "iov")
    if (is.null(.iov)) return(FALSE)
    if (is.function(.iov)) return(isTRUE(.iov(control)))
    return(isTRUE(.iov))
  }
  FALSE
}

.nlmixr2iov <- function(val, type, transform) {
  # get the standard deviation
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
#' @return nothing, called for side effects
#' @noRd
#' @author Matthew L. Fidler
.uiApplyIov <- function(ui, est, data, control) {
  if (!.isIovMethod(est, control)) {
    .uiIovEnv$ui <- NULL
    .uiIovEnv$iovDrop <- NULL
    .uiIovEnv$iovVars <- NULL
    .uiIovEnv$iovRename <- NULL
    .uiIovEnv$lines <- NULL
    .uiIovEnv$muModel <- NULL
    return(NULL)
  }
  .uiIovEnv$iovVars <- NULL
  .uiIovEnv$muModel <- NULL
  .xform <- control$iovXform
  if (length(.xform)  != 1) {
    .xform <- "sd"
  }
  if (!(.xform %in% c("sd", "var", "logsd", "logvar"))) {
    .xform <- "sd"
  }
  .ui <- ui
  .iniDf <- .ui$iniDf
  .wOcc <- which(!is.na(.iniDf$condition) &
                   .iniDf$condition != "id" &
                   is.na(.iniDf$err))
  # the expansion gives each occasion parameter its OWN magnitude theta and
  # per-occasion unit-variance etas, which cannot represent a correlation
  # between two of them; an off-diagonal row here would otherwise be treated
  # as one more occasion parameter named "(iov.a,iov.b)"
  .wOff <- .wOcc[which(!is.na(.iniDf$neta1[.wOcc]) &
                         .iniDf$neta1[.wOcc] != .iniDf$neta2[.wOcc])]
  if (length(.wOff) > 0L) {
    stop("correlated inter-occasion random effects are not supported: ",
         paste0("'", .iniDf$name[.wOff], "'", collapse=", "),
         "; give each occasion parameter its own variance",
         call.=FALSE)
  }
  # one entry per occasion VARIABLE, however many parameters ride on it
  .lvls <- unique(.iniDf$condition[.wOcc])

  .uiIovEnv$iovRename <- NULL
  if (length(.lvls) > 0) {
    .n <- .iniDf[which(.iniDf$condition %in% .lvls), "name"]
    .ui <- suppressWarnings(eval(str2lang(paste0("rxode2::rxRename(.ui, ",
                                paste(paste0("rx.", .n, "=", .n),
                                      collapse=", "), ")"))))
    .uiIovEnv$iovRename <- str2lang(paste0("rxode2::rxRename(.ui, ",
                                          paste(paste0(.n, "=", "rx.", .n),
                                      collapse=", "), ")"))
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
    # Both rows are COPIES of an unrelated parameter's row used as a
    # template; a `prior` carried over from it belongs to that parameter,
    # not to the IOV magnitude / occasion eta this becomes.  The magnitude
    # gets the prior declared on the occasion eta itself below.
    if (any(names(.theta1) == "prior")) .theta1$prior <- NA_character_
    if (any(names(.eta1) == "prior")) .eta1$prior <- NA_character_

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

    .data <- data
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
    .env$extraThetas <- NULL
    .env$extraEtas <- NULL
    .lines <- lapply(names(.lvls),
                     function(l1) {
                       .w <-which(.iniDf$condition == l1)
                       .var <- .iniDf$name[.w]
                       .lst <- c(lapply(.var, function(v) {
                         # Add theta to dataset; represents variance of iov,
                         # converted below based on the xform
                         .curTheta <- .theta1
                         # this occasion parameter's OWN row -- `.var` can
                         # hold several parameters riding on one occasion
                         # variable, so every per-parameter field has to be
                         # read from here rather than from a vector over the
                         # whole condition
                         .wv <- which(.iniDf$name == v & is.na(.iniDf$ntheta))
                         if (length(.wv) != 1L) {
                           # rxode2 does build a ui with the same occasion
                           # parameter declared twice; every field read from
                           # .wv below would then be a vector, and the
                           # rewrite would die on "replacement has 2 rows,
                           # data has 1" without naming the parameter
                           stop("the inter-occasion parameter '", v,
                                "' has ", length(.wv), " variance ",
                                "declarations; it needs exactly one",
                                call.=FALSE)
                         }
                         .est <- .iniDf[.wv, "est"]
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
                         .curTheta$fix <- .iniDf$fix[.wv]
                         # the prior the user declared with `prior(iov.x)`
                         # describes this magnitude: carry it across the
                         # rewrite, which deletes the row it was written on.
                         # No `fix` guard is needed -- lotri refuses a prior
                         # on a fixed parameter while the ui is built, so a
                         # fixed magnitude reaching here always has none.
                         if (any(names(.curTheta) == "prior")) {
                           .curTheta$prior <- .iniDf$prior[.wv]
                         }

                         .w <- which(ui$muRefCurEval$parameter == v)
                         if (length(.w) == 1L) {
                           .curEval <- ui$muRefCurEval$curEval[.w]
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

                         .env$extraThetas <- c(.env$extraThetas, .curTheta)
                         for (n in .lvls[[l1]]) {
                           .curEta <- .eta1
                           .curEta$name <- paste0("rx.", v, ".", n)
                           .curEta$label <- paste0(v, "(", l1, "==", n, ")")
                           .env$drop <- c(.env$drop, .curEta$name)
                           .env$maxeta <- .curEta$neta1 <-
                             .curEta$neta2 <- .env$maxeta + 1L
                           .env$etas <- rbind(.env$etas, .curEta)
                           .env$extraEtas <- c(.env$extraEtas, .curEta)
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
                           # |theta| written as sqrt(theta^2): identical value,
                           # but symengine differentiates it EXACTLY
                           # (theta/sqrt(theta^2) = sign) -- abs() is rewritten
                           # to an indicator form whose derivative is a
                           # tanh-SMOOTHED step, making the theta-sensitivity
                           # column wrong by a theta-dependent factor (#952)
                           str2lang(paste0("rx.", v, " <- sqrt((", v, ")^2)*(",
                                           paste(paste0("rx.", v, ".", .lvls[[l1]],
                                                        "*(", l1,
                                                        " == ", .lvls[[l1]], ")"),
                                                 collapse="+"),
                                           ")"))
                         } else if (.xform == "var") {
                           # sqrt(|theta|) as (theta^2)^(1/4), same reasoning
                           # as the "sd" branch above (#952)
                           str2lang(paste0("rx.", v, " <- ((", v, ")^2)^0.25*(",
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
    .uiIovEnv$lines <- do.call(`c`, .lines)
    .lines <- c(.uiIovEnv$lines, .ui$lstExpr)
    .ui <- rxode2::rxUiDecompress(.ui)
    # Now the lines can be added to the model
    assign("iniDf", rbind(.env$thetas,.env$etas), envir = .ui)
    assign("lstExpr", .lines, envir = .ui)
    .uiIovEnv$ui <- ui
    .uiIovEnv$iovDrop <- .env$drop # extra variables to drop
    list(ui = rxode2::rxUiDecompress(suppressWarnings(suppressMessages(.ui$fun()))))
  } else {
    .uiIovEnv$ui <- NULL
    .uiIovEnv$iovDrop <- NULL
    NULL
  }
}
#' Finalizes IOV model
#'
#' @param ret data frame with some iov information dropped
#' @return fit with iov information dropped
#' @noRd
#' @author Matthew L. Fidler
.uiFinalizeIov <- function(ret) {
  if (!is.null(.uiIovEnv$ui)) {
    if (is.null(ret$ui)) return(ret)

    if (is.environment(ret$env)) {
      #.preFinalParTableHooksRun(.uiIovEnv)
      .ui <- ret$env$ui
      .lstExpr <- .ui$lstExpr
      .w <- which(vapply(seq_along(.lstExpr),
             function(i) {
               any(vapply(seq_along(.uiIovEnv$lines),
                      function(j) {
                        if (identical(.lstExpr[[i]], .uiIovEnv$lines[[j]])) {
                          return(TRUE)
                        }
                        FALSE
                      }, logical(1), USE.NAMES = FALSE))
             }, logical(1), USE.NAMES = FALSE))
      if (length(.w) > 0L) {
        .lstExpr <- lapply(seq_along(.lstExpr)[-.w],
                           function(i) {
                             .lstExpr[[i]]
                           })
      }
      # Get the IOV variables that are present as thetas in the model
      .iovName <- new.env(parent=emptyenv())
      .iovDf <- .uiIovEnv$ui$iniDf
      .iovDf <- .iovDf[!is.na(.iovDf$neta1) & .iovDf$condition != "id",, drop=FALSE]
      .iovName$var <- .iovDf$name

      getEstimateDf <- function(iniDf) {
        .iniDf <- iniDf
        # Final thetaDf & etaDf
        .thetaDf <- .iniDf[is.na(.iniDf$neta1),, drop=FALSE]
        .etaDf <- .iniDf[!is.na(.iniDf$neta1),, drop=FALSE]

        # Drop the dummy etas
        .etaDf <- .etaDf[!(.etaDf$name %in% .uiIovEnv$iovDrop),, drop=FALSE]

        # Renumber etas, just in case
        .etaDf$neta1 <- factor(.etaDf$neta1)
        .etaDf$neta2 <- factor(.etaDf$neta2, levels=levels(.etaDf$neta1))
        .etaDf$neta1 <- as.integer(.etaDf$neta1)
        .etaDf$neta2 <- as.integer(.etaDf$neta2)

        .maxEta <- if (nrow(.etaDf) == 0L) 0L else max(.etaDf$neta1)


        # Go through each IOV variable and calculate the variance from the back-transform
        # Add it to the .etaDf afterward, and remove from .thetaDf
        # Use a template row; fall back to .iniDf when .etaDf is empty (all ETAs
        # were IOV dummy ETAs that got dropped above).
        .etaTemplate <- if (nrow(.etaDf) > 0L) {
          .etaDf[1, , drop = FALSE]
        } else {
          .tmp <- .iniDf[1, , drop = FALSE]
          .tmp$ntheta <- NA_real_
          .tmp
        }
        for (i in seq_along(.iovDf$name)) {
          .w <- which(.thetaDf$name == .iovDf$name[i])
          .fun <- sub("Cv$", "Sd", .thetaDf[.w, "backTransform"])
          .fun <- get(.fun)
          .est <- .fun(.thetaDf[.w, "est"])^2
          .maxEta <- .maxEta + 1L
          .cur <- .etaTemplate
          .cur$neta1 <- .cur$neta2 <- .maxEta
          .cur$est <- .est
          .cur$fix <- .iovDf$fix[i]
          .cur$upper <- .iovDf$upper[i]
          .cur$lower <- .iovDf$lower[i]
          .cur$label <- .iovDf$label[i]
          .cur$backTransform <- .iovDf$backTransform[i]
          .cur$err <- .iovDf$err[i]
          .cur$name <- paste0("rx.", .iovDf$name[i]) # Matches replacement
          .cur$condition <- .iovDf$condition[i]
          # .etaTemplate is a COPY of the first remaining eta (or of the
          # first iniDf row): its prior belongs to that parameter.  This row
          # is the user's own `iov.x ~ v | occ` being restored, so it takes
          # the prior they wrote on it.
          if (any(names(.cur) == "prior")) {
            .cur$prior <- .iovDf$prior[i]
          }
          .etaDf <- rbind(.etaDf, .cur)
          .thetaDf <- .thetaDf[-.w, , drop=FALSE]
        }

        # Renumber
        .thetaDf$ntheta <- as.integer(factor(.thetaDf$ntheta))
        rbind(.thetaDf, .etaDf)
      }

      .finalDf <- getEstimateDf(.ui$iniDf)
      .iniDf0 <- getEstimateDf(ret$env$iniDf0)
      assign("iniDf0", .iniDf0, envir = ret$env)

      # Save with final estimates
      .ini <- as.expression(lotri::as.lotri(.finalDf))
      .ini[[1]] <- quote(`ini`)
      .model <- rxode2::as.model(.lstExpr)
      .ui <- .getUiFunFromIniAndModel(.ui, .ini, .model)
      # apply renames and evaluate new model
      .ui <- eval(.uiIovEnv$iovRename)
      assign("ui", .ui, envir = ret$env)

      .finalDf <- .ui$iniDf
      .est <- setNames(.finalDf$est, .finalDf$name)

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
                  ## nocov start
                  # simply testing edge cases with warnings
                  if(length(.w)==0L) {
                    warning("IOV variable '", d,
                            "' is not present in the shrinkage results, check your model and dataset",
                            call. = FALSE)
                    return(NULL)
                  }
                  if (length(.w)==1L) {
                    warning("IOV variable '", d,
                            "' has the same number of levels as the random effect, check your dataset",
                            call. = FALSE)
                  }
                  ## nocov end
                  .curd <- .shrinkN[,.w,  drop=FALSE]

                  # combine per-level mean/var/skewness/kurtosis by averaging
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
                                          2*stats::pt(-abs(.tv), df = (.nid*.nv) -1)
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
          # rescale the derived eta (fixed to 1) by the IOV variable's sd
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
      .cls <- class(ret)
      if (length(.w) > 0L) {
        class(ret) <- "data.frame"
        ret <- ret[,-.w]
      }
      .rename <- paste0(.uiIovEnv$iovVars, ".rx")
      names(ret) <- vapply(names(ret), function(n) {
        if (n %in% .rename) {
          sub("[.]rx$", "", n)
        } else {
          n
        }
      }, character(1), USE.NAMES = FALSE)
      class(ret) <- .cls
    }
    .stripFastmatchHash(ret$env)
  }
  ret
}

preProcessHooksAdd(".uiApplyIov", .uiApplyIov)
postFinalObjectHooksAdd(".uiFinalizeIov", .uiFinalizeIov)
