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

#' Reason a method's own IOV handling cannot cover this model
#'
#' Read from the optional \code{"iovNativeScope"} attribute on the
#' \code{nlmixr2Est.<method>} S3 method, a \code{function(ui, data, control)}
#' returning \code{NULL} when the method's own handling covers the model and a
#' short reason when it does not.  A reason makes the shared rewrite run after
#' all, so a model outside the newer handling still fits.
#'
#' @param ui rxode2 user interface model
#' @param est estimation routine name
#' @param data data set for the fit
#' @param control control object
#' @return \code{NULL}, or a short character reason
#' @noRd
#' @author Matthew L. Fidler
.iovNativeDecline <- function(ui, est, data, control) {
  .v <- as.character(utils::methods("nlmixr2Est"))
  if (!(paste0("nlmixr2Est.", est) %in% .v)) return(NULL)
  .f <- attr(utils::getS3method("nlmixr2Est", est), "iovNativeScope")
  if (!is.function(.f)) return(NULL)
  .f(ui, data, control)
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
#' Estimation methods that honour a repeated (`same()`) omega block
#'
#' `iovMethod="omega"` only means anything for a method that actually
#' shares the repeated block's parameters while estimating.  That
#' sharing lives in `.foceiOptEnvSetupBounds()` (`R/focei.R`), which
#' passes `same=` to `rxSymInvCholCreate()`, so the FOCEi family gets it
#' and nothing else does.
#'
#' SAEM in particular must NOT be given this expansion: its omega is a
#' moment M-step in C++ (`Gamma2_phi1`) with no notion of a shared
#' block, so each occasion's block would be estimated independently and
#' `.uiFinalizeIov()` would then report occasion one and silently
#' discard the rest.  The same goes for the variational (`vae`, `fbvi`,
#' `emvi`), nonparametric (`npag`, `npb`, ...) and importance-sampling
#' (`imp`, `impmap`, `qrpem`) methods, which build their own omega.
#'
#' A method left out here keeps the long standing behaviour: a
#' correlated occasion block is refused outright.
#'
#' @noRd
.iovSameMethods <- c("focei", "foce", "focep", "laplace", "agq",
                     "foceif", "focef", "focepf", "agqf",
                     "ifoce", "ifocep", "ifocei", "ilaplace", "iagq",
                     "ifoceif", "ifocef", "ifocepf", "iagqf",
                     "mfoce", "mfocep", "mfocei", "mlaplace", "magq",
                     "mfoceif", "mfocef", "mfocepf", "magqf")

#' Does this estimation method honour a repeated (`same()`) omega block?
#'
#' @param est estimation method name
#' @return TRUE when the shared-block expansion may be used
#' @noRd
#' @author Matthew L. Fidler
.isIovSameMethod <- function(est) {
  if (length(est) != 1L || is.na(est)) return(FALSE)
  ## an out-of-tree method can opt in the way it declares `iov` itself
  .v <- as.character(utils::methods("nlmixr2Est"))
  if (paste0("nlmixr2Est.", est) %in% .v) {
    .a <- attr(utils::getS3method("nlmixr2Est", est), "iovSame")
    if (isTRUE(.a)) return(TRUE)
  }
  est %in% .iovSameMethods
}

#' Build the per-occasion eta blocks for `iovMethod="omega"`
#'
#' Occasion one IS the estimated block: it takes the variances and
#' covariances the user wrote on the `| occ` rows.  Every later occasion
#' repeats it, recorded the way `lotri`'s `same()` records a repeated
#' block -- in the `condition` column, pointing at the element it
#' mirrors -- so only one set of parameters is estimated however many
#' occasions there are.
#'
#' The etas are laid out occasion-major, with the parameters inside each
#' occasion, so that each occasion's block is contiguous.  The "theta"
#' expansion lays them out parameter-major instead, which is fine there
#' because each is an independent 1x1.
#'
#' @param var occasion parameter names riding on this level
#' @param lvls the observed levels of the occasion variable
#' @param l1 the occasion variable's name
#' @param iniDf the ORIGINAL ini data frame, read for the user's block
#' @param eta1 template eta row
#' @param env accumulating environment (`etas`, `maxeta`, `drop`,
#'   `extraEtas`)
#' @return nothing, called for its effect on `env`
#' @noRd
#' @author Matthew L. Fidler
.uiIovOmegaEtas <- function(var, lvls, l1, iniDf, eta1, env) {
  .nm <- function(v, n) paste0("rx.", v, ".", n)
  ## the user's block for this level, read by name
  .diagRow <- function(v) which(iniDf$name == v & is.na(iniDf$ntheta))[1]
  .diagEst <- vapply(var, function(v) iniDf$est[.diagRow(v)],
                     double(1), USE.NAMES = TRUE)
  ## `iov.cl ~ fix(0.1) | occ` fixes the VARIANCE.  Under "theta" that
  ## rides on the magnitude theta; here the variance is the omega block,
  ## so the flag has to land on the block or the parameter is quietly
  ## estimated while `.uiFinalizeIov()` still reports it as fixed.
  .diagFix <- vapply(var, function(v) isTRUE(iniDf$fix[.diagRow(v)]),
                     logical(1), USE.NAMES = TRUE)
  ## `as.data.frame()` names an off diagonal "(smaller,larger)" by
  ## matrix position, so try both spellings rather than guessing
  .covRow <- function(a, b) {
    .w <- which(iniDf$name %in% c(paste0("(", a, ",", b, ")"),
                                  paste0("(", b, ",", a, ")")) &
                  is.na(iniDf$ntheta))
    if (length(.w) == 0L) return(NA_integer_)
    .w[1]
  }
  .covEst <- function(a, b) {
    .w <- .covRow(a, b)
    if (is.na(.w)) return(0)
    iniDf$est[.w]
  }
  .covFix <- function(a, b) {
    .w <- .covRow(a, b)
    if (is.na(.w)) return(FALSE)
    isTRUE(iniDf$fix[.w])
  }
  .masterOf <- character(0)
  for (.oi in seq_along(lvls)) {
    .n <- lvls[.oi]
    .first <- env$maxeta + 1L
    ## diagonal rows for this occasion
    for (.vi in seq_along(var)) {
      .v <- var[.vi]
      .cur <- eta1
      .cur$name <- .nm(.v, .n)
      ## NO per-occasion label: `same()` is only re-emitted when the
      ## copy block matches its master in values, `fix` AND labels, and a
      ## distinct "iov.cl(occ==2)" label would silently break the
      ## repetition on the next `.ui$fun()` round trip.  These etas are
      ## internal; `.uiFinalizeIov()` restores the user's own labeled rows.
      .cur$label <- NA_character_
      .cur$fix <- .diagFix[[.v]]
      .cur$est <- .diagEst[[.v]]
      env$drop <- c(env$drop, .cur$name)
      env$maxeta <- .cur$neta1 <- .cur$neta2 <- env$maxeta + 1L
      .cur$condition <- if (.oi == 1L) {
        "id"
      } else {
        paste0("id:same:", .masterOf[[.v]])
      }
      env$etas <- rbind(env$etas, .cur)
      env$extraEtas <- c(env$extraEtas, .cur)
      if (.oi == 1L) {
        .masterOf[[.v]] <- .cur$name
        ## the block `.uiFinalizeIov()` reads the estimate back off,
        ## since the magnitude theta is fixed at one in this mode
        .uiIovEnv$iovMaster[[.v]] <- .cur$name
      }
    }
    ## and the covariances within it
    if (length(var) > 1L) {
      for (.i in seq_along(var)[-1L]) {
        for (.j in seq_len(.i - 1L)) {
          .est <- .covEst(var[.i], var[.j])
          if (.est == 0) next
          .cur <- eta1
          .cur$name <- paste0("(", .nm(var[.j], .n), ",",
                              .nm(var[.i], .n), ")")
          .cur$label <- NA_character_
          .cur$fix <- .covFix(var[.i], var[.j])
          .cur$est <- .est
          .cur$neta1 <- .first + .i - 1L
          .cur$neta2 <- .first + .j - 1L
          .cur$condition <- if (.oi == 1L) {
            "id"
          } else {
            paste0("id:same:", .masterOf[[var[.j]]], ":",
                   .masterOf[[var[.i]]])
          }
          env$etas <- rbind(env$etas, .cur)
        }
      }
    }
  }
  invisible()
}

#' This applies the IOV method to the model based on the data used
#'
#' @return nothing, called for side effects
#' @noRd
#' @author Matthew L. Fidler
.uiApplyIov <- function(ui, est, data, control) {
  .fellBack <- FALSE
  if (!.isIovMethod(est, control)) {
    # the method's own IOV handling is in force -- unless this model is outside
    # its scope, in which case fall back to the rewrite below so the fit still
    # runs.  The reason is short because it is collected into the fit's $runInfo.
    .why <- .iovNativeDecline(ui, est, data, control)
    if (is.null(.why)) {
      .uiIovEnv$ui <- NULL
      .uiIovEnv$iovDrop <- NULL
      .uiIovEnv$iovVars <- NULL
      .uiIovEnv$iovRename <- NULL
      .uiIovEnv$lines <- NULL
      .uiIovEnv$muModel <- NULL
      .uiIovEnv$iovTwoLevel <- NULL
      .uiIovEnv$iovCollapsed <- NULL
      # this env outlives one fit, so a stale "omega" here would make the
      # NEXT model's finalize read estimates off a block that is not there
      .uiIovEnv$iovMethod <- NULL
      .uiIovEnv$iovMaster <- list()
      return(NULL)
    }
    warning(.why, "; used iovMethod='theta'", call.=FALSE)
    control$iovMethod <- "theta"
    .fellBack <- TRUE
  }
  .uiIovEnv$iovVars <- NULL
  .uiIovEnv$muModel <- NULL
  .uiIovEnv$iovTwoLevel <- NULL
  .uiIovEnv$iovCollapsed <- NULL
  .uiIovEnv$iovMaster <- list()
  .xform <- control$iovXform
  if (length(.xform)  != 1) {
    .xform <- "sd"
  }
  if (!(.xform %in% c("sd", "var", "logsd", "logvar"))) {
    .xform <- "sd"
  }
  .ui <- ui
  .iniDf <- .ui$iniDf
  ## the level is the BASE condition: a repeated (`same()`) block carries a
  ## `:same:<master>` suffix, which is not a level of variability
  .baseCnd <- lotri::lotriBaseCondition(.iniDf$condition)
  .wOcc <- which(!is.na(.iniDf$condition) &
                   .baseCnd != "id" &
                   is.na(.iniDf$err))
  .wOff <- .wOcc[which(!is.na(.iniDf$neta1[.wOcc]) &
                         .iniDf$neta1[.wOcc] != .iniDf$neta2[.wOcc])]
  ## How the occasion parameters are expanded before estimation.
  ##
  ## "theta" is the long standing shape: one magnitude theta per occasion
  ## parameter, with per-occasion unit-variance etas fixed to it.  It
  ## cannot represent a correlation between two occasion parameters -- an
  ## off diagonal row would be treated as one more occasion parameter
  ## named "(iov.a,iov.b)".
  ##
  ## "omega" fixes the magnitude at one and estimates the per-occasion eta
  ## blocks instead, occasion one being the block and the rest repeating
  ## it, so a correlation is carried by the block itself.
  ## `control$iovMethod` is overloaded: `saemControl()` uses it for SAEM's
  ## NATIVE handling ("twoLevel"/"collapsed"/"theta"), `foceiControl()` for
  ## which REWRITE to use ("auto"/"theta"/"omega").  Only read it as the
  ## latter for a method that can honour a shared block; for anything else
  ## the rewrite is the long standing one, whatever the native value said.
  ## (A native handler that declined has already set it to "theta" above.)
  .sameOk <- .isIovSameMethod(est)
  ## asked for outright on a method that ignores it: refuse by name rather
  ## than silently downgrade, since the fit would look fine and estimate
  ## each occasion's block on its own
  if (!.sameOk && identical(control$iovMethod, "omega")) {
    stop("'iovMethod=\"omega\"' repeats one estimated omega block per ",
         "occasion, which the estimation method '", est, "' does not ",
         "honour; use a FOCEi family method for correlated ",
         "inter-occasion random effects",
         call.=FALSE)
  }
  .iovMethod <- if (.sameOk) control$iovMethod else "theta"
  if (length(.iovMethod) != 1L ||
        !(.iovMethod %in% c("auto", "theta", "omega"))) {
    .iovMethod <- if (.sameOk) "auto" else "theta"
  }
  ## Resolved PER LEVEL of variability, not per model: a correlation on
  ## `occ` is no reason to change how an unrelated diagonal `occ2` is
  ## expanded, and "theta" is the better conditioned inner problem.
  .offLvls <- unique(.baseCnd[.wOff])
  .resolveLvl <- function(l1) {
    if (.iovMethod != "auto") return(.iovMethod)
    if (l1 %in% .offLvls && .sameOk) "omega" else "theta"
  }
  .anyOmega <- .iovMethod == "omega" ||
    (.iovMethod == "auto" && length(.offLvls) > 0L && .sameOk)
  ## a level whose expansion resolves to "theta" cannot carry a
  ## correlation, however it was chosen
  .badOff <- .wOff[vapply(.baseCnd[.wOff], function(l1) {
    .resolveLvl(l1) == "theta"
  }, logical(1), USE.NAMES=FALSE)]
  if (length(.badOff) > 0L) {
    stop("correlated inter-occasion random effects are not supported",
         if (.sameOk) "" else paste0(" by '", est, "'"), ": ",
         paste0("'", .iniDf$name[.badOff], "'", collapse=", "),
         "; give each occasion parameter its own variance",
         if (.sameOk) ", or use `iovMethod=\"omega\"`" else "",
         call.=FALSE)
  }
  .uiIovEnv$iovMethod <- if (.anyOmega) "omega" else "theta"
  if (.anyOmega) {
    ## the magnitude theta is fixed at one, so the way it is
    ## parameterized no longer means anything; "sd" is what the model
    ## line and the analytic covariance expect to see
    if (!identical(control$iovXform, "sd") &&
          !is.null(control$iovXform)) {
      .minfo(paste0("'iovXform' is ignored when 'iovMethod=\"omega\"'; ",
                    "the magnitude is fixed at one and the variability ",
                    "is estimated in the omega block"))
    }
  }
  # one entry per occasion VARIABLE, however many parameters ride on it
  .lvls <- unique(.baseCnd[.wOcc])

  .uiIovEnv$iovRename <- NULL
  if (length(.lvls) > 0) {
    ## DIAGONAL rows only: an occasion block may now carry covariances,
    ## and `(iov.a,iov.b)` is a cell of the block, not a parameter to
    ## rename
    .n <- .iniDf[which(.baseCnd %in% .lvls &
                         !is.na(.iniDf$neta1) &
                         .iniDf$neta1 == .iniDf$neta2), "name"]
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

    .etas <- .etas[which(!(lotri::lotriBaseCondition(.etas$condition) %in%
                             .lvls)), , drop=FALSE]
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
                       ## this LEVEL's expansion; `iovXform` parameterizes
                       ## the "theta" magnitude and means nothing once the
                       ## magnitude is fixed at one
                       .m1 <- .resolveLvl(l1)
                       .xf1 <- if (.m1 == "omega") "sd" else .xform
                       ## again diagonal rows only -- these are the
                       ## occasion PARAMETERS, not the cells of their block
                       .w <- which(.baseCnd == l1 &
                                     !is.na(.iniDf$neta1) &
                                     .iniDf$neta1 == .iniDf$neta2)
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
                         if (.xf1 == "var") {
                           .curTheta$est <- .est
                         } else if (.xf1 == "sd") {
                           .curTheta$est <- sqrt(.est)
                         } else if (.xf1 == "logvar") {
                           .curTheta$est <- log(.est)
                         } else if (.xf1 == "logsd") {
                           .curTheta$est <- log(sqrt(.est))
                         }
                         if (.m1 == "omega") {
                           ## the magnitude is carried by the omega block,
                           ## so the theta is a constant multiplier of one.
                           ## It is kept rather than dropped so the model
                           ## line, `rxUiGet.foceiSkipCov` and the theta
                           ## deletion in `.uiFinalizeIov()` all keep
                           ## working unchanged; `foceiSetupTheta_()`
                           ## removes a fixed theta from the optimizer.
                           .curTheta$est <- 1
                           .curTheta$fix <- TRUE
                           .curTheta$lower <- -Inf
                         }
                         .curTheta$name <- v
                         .uiIovEnv$iovVars <- c(.uiIovEnv$iovVars, v)
                         if (.m1 != "omega") {
                           ## in "omega" mode the magnitude is fixed at one
                           ## above and the user's `fix` belongs to the
                           ## omega block instead
                           .curTheta$fix <- .iniDf$fix[.wv]
                         }
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
                           paste0(switch(.xf1,
                                  "sd" = "nlmixr2iovSd",
                                  "var" = "nlmixr2iovVar",
                                  "logsd" = "nlmixr2iovLogsd",
                                  "logvar" = "nlmixr2iovLogvar"),
                                  ifelse(.curEval=="exp", "Cv", "Sd"))
                         if (.xf1 %in% c("sd", "var")) {
                           .curTheta$lower <- 0 # doesn't work with saem
                         }
                         .env$maxtheta <- .curTheta$ntheta <- .env$maxtheta + 1L
                         .env$thetas <- rbind(.env$thetas, .curTheta)

                         .env$extraThetas <- c(.env$extraThetas, .curTheta)
                         if (.m1 != "omega") {
                           ## "theta" mode: one unit-variance eta per
                           ## occasion, fixed, scaled by the magnitude
                           ## theta above.  Under "omega" the etas are
                           ## built occasion-major after this loop, since
                           ## each occasion's parameters form one block.
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
                         }
                         if (.xf1 == "logsd") {
                           str2lang(paste0("rx.", v, " <- exp(", v, ")*(",
                                           paste(paste0("rx.", v, ".", .lvls[[l1]],
                                                        "*(", l1,
                                                        " == ", .lvls[[l1]], ")"),
                                                 collapse="+"),
                                           ")"))
                         } else if (.xf1 == "logvar") {
                             str2lang(paste0("rx.", v, " <- sqrt(exp(", v, "))*(",
                                             paste(paste0("rx.", v, ".", .lvls[[l1]],
                                                          "*(", l1,
                                                          " == ", .lvls[[l1]], ")"),
                                                   collapse="+"),
                                             ")"))
                         } else if (.xf1 == "sd") {
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
                         } else if (.xf1 == "var") {
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
                       if (.m1 == "omega") {
                         .uiIovOmegaEtas(.var, .lvls[[l1]], l1, .iniDf,
                                         .eta1, .env)
                       }
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
    .ret <- list(ui = rxode2::rxUiDecompress(suppressWarnings(suppressMessages(.ui$fun()))))
    if (.fellBack) .ret$control <- control
    .ret
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
  # the collapsed sampler has its own restoration (.saemIovFinalizeCollapsed,
  # R/saemIov.R): it removed the user's line entirely and has no magnitude theta,
  # so almost none of the mechanism below applies
  if (!is.null(.uiIovEnv$ui) && is.null(.uiIovEnv$iovCollapsed)) {
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
        # A correlated occasion block also produced `(a,b)` covariance rows.
        # Those are not in `iovDrop` -- it names variables, and a covariance
        # is not one -- so drop whatever still points at an eta that just
        # went away, or it survives with a dangling `neta2` and the
        # regenerated `ini({})` no longer parses.
        .keptEta <- .etaDf$neta1[.etaDf$neta1 == .etaDf$neta2]
        .etaDf <- .etaDf[.etaDf$neta1 %in% .keptEta &
                           .etaDf$neta2 %in% .keptEta, , drop=FALSE]

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
        ## PER VARIABLE, since the expansion is chosen per level of
        ## variability: `iovMaster` holds an entry exactly for the
        ## occasion parameters that got the shared-block treatment
        .omegaVar <- function(v) !is.null(.uiIovEnv$iovMaster[[v]])
        # in "omega" mode the magnitude theta is FIXED AT ONE and the
        # variability lives in the shared (`same()`) block, so the estimate
        # is read off that block's master row rather than back-transformed
        # out of the theta
        .masterEst <- function(nm) {
          .wm <- which(.iniDf$name == nm)
          ## nocov start
          if (length(.wm) != 1L) {
            stop("cannot find the occasion master estimate '", nm, "'",
                 call.=FALSE)
          }
          ## nocov end
          .iniDf$est[.wm]
        }
        # a correlated occasion block puts off-diagonal `(a,b)` rows in
        # `.iovDf` as well; they can only be numbered once every diagonal
        # they point at has a new eta number, so walk the diagonals first
        .diagI <- which(.iovDf$neta1 == .iovDf$neta2)
        .offI <- which(.iovDf$neta1 != .iovDf$neta2)
        .newEta <- integer(0)
        .rmTheta <- function(w) {
          # `x[-integer(0)]` is EMPTY, not everything -- never let an
          # unmatched name silently drop every theta
          if (length(w) == 1L) .thetaDf <<- .thetaDf[-w, , drop=FALSE]
        }
        .fillRow <- function(i, cur) {
          cur$fix <- .iovDf$fix[i]
          cur$upper <- .iovDf$upper[i]
          cur$lower <- .iovDf$lower[i]
          cur$label <- .iovDf$label[i]
          cur$backTransform <- .iovDf$backTransform[i]
          cur$err <- .iovDf$err[i]
          cur$condition <- .iovDf$condition[i]
          # .etaTemplate is a COPY of the first remaining eta (or of the
          # first iniDf row): its prior belongs to that parameter.  This row
          # is the user's own `iov.x ~ v | occ` being restored, so it takes
          # the prior they wrote on it.
          if (any(names(cur) == "prior")) {
            cur$prior <- .iovDf$prior[i]
          }
          cur
        }
        for (i in .diagI) {
          .v <- .iovDf$name[i]
          .w <- which(.thetaDf$name == .v)
          if (.omegaVar(.v)) {
            # shared (`same()`) block: the magnitude theta is fixed at one,
            # so the variance is the block's own master row
            .est <- .masterEst(.uiIovEnv$iovMaster[[.v]])
          } else if (is.null(.uiIovEnv$iovTwoLevel)) {
            # shared rewrite: the variance is a magnitude theta on the
            # iovXform scale, converted back through its own back-transform
            .fun <- sub("Cv$", "Sd", .thetaDf[.w, "backTransform"])
            .fun <- get(.fun)
            .est <- .fun(.thetaDf[.w, "est"])^2
          } else {
            # two-level: the variance IS an omega entry already, shared by every
            # occasion level (poolOmegaGroups, src/saem.cpp), so read it off the
            # first one.  This runs over the FINAL frame and over iniDf0, and
            # iniDf0 is the user's own frame -- it never went through the
            # expansion, so its `iov.x ~ v | occ` row is already what we would
            # be rebuilding.  Leave it alone.
            .w <- integer(0)
            .poolEta <- .uiIovEnv$iovTwoLevel[[.v]][1]
            if (!(.poolEta %in% .iniDf$name)) next
            .est <- .iniDf$est[.iniDf$name == .poolEta]
          }
          .maxEta <- .maxEta + 1L
          .newEta[[.v]] <- .maxEta
          .cur <- .etaTemplate
          .cur$neta1 <- .cur$neta2 <- .maxEta
          .cur$est <- .est
          .cur$name <- paste0("rx.", .v) # Matches replacement
          .etaDf <- rbind(.etaDf, .fillRow(i, .cur))
          .rmTheta(.w)
        }
        for (i in .offI) {
          .nm <- .iovDf$name[i]
          .pair <- strsplit(substr(.nm, 2L, nchar(.nm) - 1L), ",",
                            fixed=TRUE)[[1]]
          ## nocov start
          if (length(.pair) != 2L || !all(.pair %in% names(.newEta))) {
            stop("cannot restore the occasion covariance '", .nm, "'",
                 call.=FALSE)
          }
          ## nocov end
          .ma <- .uiIovEnv$iovMaster[[.pair[1]]]
          .mb <- .uiIovEnv$iovMaster[[.pair[2]]]
          # `as.data.frame()` names an off diagonal by matrix position, so
          # accept either spelling rather than guessing the order
          .wo <- which(.iniDf$name %in% c(paste0("(", .ma, ",", .mb, ")"),
                                          paste0("(", .mb, ",", .ma, ")")))
          if (length(.wo) == 0L) next
          .e1 <- .newEta[[.pair[1]]]
          .e2 <- .newEta[[.pair[2]]]
          .cur <- .etaTemplate
          .cur$neta1 <- max(.e1, .e2)
          .cur$neta2 <- min(.e1, .e2)
          .cur$est <- .iniDf$est[.wo[1]]
          .cur$name <- paste0("(rx.", .pair[1], ",rx.", .pair[2], ")")
          .etaDf <- rbind(.etaDf, .fillRow(i, .cur))
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

      ## The per-occasion columns of an occasion parameter `d` are named
      ## `rx.<d>.<occ>`.  A SUBSTRING match also picks up `rx.<d>2.<occ>`
      ## whenever one parameter's name is a prefix of another's
      ## (`iov.v` and `iov.v2`), which silently mixed the two
      ## parameters' columns together and left the occasion NA.
      .iovOccCols <- function(d, nms) {
        which(grepl(paste0("^rx\\.",
                           gsub(".", "\\.", d, fixed=TRUE),
                           "\\.[^.]+$"), nms))
      }
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
                  .w <- .iovOccCols(d, names(.shrinkN))
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

      # the shared rewrite's occasion etas are unit-variance, so they have to be
      # rescaled by the fitted SD; the two-level ones already carry c_ik
      .sdIov <- if (is.null(.uiIovEnv$iovTwoLevel)) sqrt(.est) else NULL

      .dt <- NULL
      .omegaModeFin <- function(d) !is.null(.uiIovEnv$iovMaster[[d]])
      .iov <- lapply(.n, function(var) {
        .cur <- .omega[[var]]
        .dn <- dimnames(.cur)[[1]]
        # `.dt` is created on the FIRST dimname and later ones only add a
        # column, so its occasion-label column always carries that first
        # parameter's `rx.<name>.<occ>` spellings
        .dFirst <- .dn[1]
        for (d in .dn) {
          .w <- c(1L, .iovOccCols(d, names(.iov)))
          .curd <- data.table::data.table(.iov[,.w])
          .curd <- data.table::melt(.curd,
                                    id.vars=names(.curd)[1],
                                    measure.vars=names(.curd)[-1],
                                    variable.name = var,
                                    value.name = d)
          # rescale the derived eta (fixed to 1) by the IOV variable's sd.
          # Under the shared-block expansion the eta already carries the
          # whole magnitude -- its theta is fixed at one -- so there is
          # nothing to rescale; the two-level path sets `.sdIov` NULL for
          # the same reason.
          if (!is.null(.sdIov) && !.omegaModeFin(d)) {
            .curd[[d]] <- .curd[[d]] * .sdIov[d]
          }
          if (is.null(.dt)) {
            .dt <- .curd
          } else {
            .dt[[d]] <- .curd[[d]]
          }
        }
        # strip the FIRST parameter's prefix -- `d` here is whatever the
        # loop above left behind, and with two occasion parameters that is
        # the wrong one, which silently turned every occasion into NA
        .dt[[var]] <- as.integer(sub(paste0("rx.", .dFirst, "."), "",
                                     as.character(.dt[[var]]), fixed=TRUE))
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
      # `x[-integer(0)]` empties the vector, and the two-level path has no
      # magnitude theta in fixef to begin with
      if (length(.w) > 0L) {
        .fixef <- .fixef[-.w]
        assign("fixef", .fixef, envir = ret$env)
      }

      .parFixedDf <- ret$env$parFixedDf
      .bck <- which(grepl("Back",names(.parFixedDf)))
      .bsv <- which(grepl("BSV", names(.parFixedDf)))
      .est <- which(grepl("Est", names(.parFixedDf)))

      .hasIovTheta <- all(.uiIovEnv$iovVars %in% rownames(.parFixedDf))
      .valCharPrep <- NULL
      if (.hasIovTheta) {
        # For a shared-block occasion parameter the magnitude theta is fixed
        # at one, so its back-transformed cell is `nlmixr2iovSdCv(1)` -- a
        # constant 131%, not an estimate.  Take the CV from the variance
        # that WAS estimated: the restored occasion block on the final
        # `iniDf`.  Only those parameters; a level expanded the old way
        # still reports its magnitude theta.
        .omegaVars <- intersect(.uiIovEnv$iovVars, names(.uiIovEnv$iovMaster))
        if (length(.omegaVars) > 0L) {
          .finIni <- ret$env$ui$iniDf
          .parFixedDf[.omegaVars, .bck] <-
            vapply(.omegaVars, function(.v) {
              .wv <- which(.finIni$name == .v & !is.na(.finIni$neta1) &
                             .finIni$neta1 == .finIni$neta2)
              if (length(.wv) != 1L) return(NA_real_)  # nocov
              nlmixr2iovSdCv(sqrt(.finIni$est[.wv]))
            }, double(1), USE.NAMES=FALSE)
        }
        .valCharPrep <-
          .parFixedDf[.uiIovEnv$iovVars,.bsv] <-
          .parFixedDf[.uiIovEnv$iovVars, .bck]
        .parFixedDf[.uiIovEnv$iovVars,.bsv] <- NA_real_
        .parFixedDf[.uiIovEnv$iovVars,.est] <- NA_real_
      }

      .parFixedDf <- .parFixedDf[!grepl("^rx[.]", rownames(.parFixedDf)),]
      assign("parFixedDf", .parFixedDf, envir = ret$env)

      .parFixed <- ret$env$parFixed
      .bck2 <- which(grepl("Back", names(.parFixed)))
      .bsv2 <- which(grepl("BSV",  names(.parFixed)))
      .est2 <- which(grepl("Est", names(.parFixed)))

      .sigdig <- ret$control$sigdig
      if (.hasIovTheta) {
        .parFixed[.uiIovEnv$iovVars, .bck2] <- ""
        .parFixed[.uiIovEnv$iovVars, .est2] <- ""
        .parFixed[.uiIovEnv$iovVars, .bsv2] <- formatC(
          signif(.valCharPrep, digits = .sigdig),
          digits = .sigdig, format = "fg", flag = "#")
      }
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
