#' Generic for nlmixr2 estimation methods
#'
#' @param env Environment for the nlmixr2 estimation routines.
#'
#' This needs to have:
#'
#' - rxode2 ui object in `$ui`
#'
#' - data to fit in the estimation routine in `$data`
#'
#' - control for the estimation routine's control options in `$ui`
#'
#' @param ... Other arguments provided to `nlmixr2Est()` provided for
#'   flexibility but not currently used inside nlmixr
#'
#' @return nlmixr2 fit object
#'
#' @author Matthew Fidler
#'
#' @details
#'
#' This is a S3 generic that allows others to use the nlmixr2
#'   environment to do their own estimation routines
#'
#' @export
nlmixr2Est <- function(env, ...) {
  on.exit({
    .nlmixr2clearPipe()
    nlmixr2global$nlmixr2SimInfo <- NULL
  })
  if (!exists("ui", envir=env)) {
    stop("need 'ui' object", call.=FALSE)
  } else if (!inherits(get("ui", envir=env), "rxUi")) {
    stop("'ui' is not an rxode2 object", call.=FALSE)
  }
  if (!inherits(env, "output")) {
    nlmixr2global$nlmixr2EstEnv$iniDf0 <- data.frame(get("ui", envir=env)$iniDf)
  }
  if (!exists("data", envir=env)) {
    stop("need 'data' object", call.=FALSE)
  } else if (!inherits(get("data", envir=env), "data.frame")) {
    stop("'data' is not a data.frame", call.=FALSE)
  }
  assign("data", as.data.frame(get("data", envir=env)), envir=env)
  if (!exists("control", envir=env)) {
    stop("need 'control' object", call.=FALSE)
  } else if (is.null(get("control", envir=env))) {
  } else {
  }
  if (!exists("table", envir=env)) {
    stop("need 'table' object", call.=FALSE)
  } else if (is.null(get("table", envir=env))) {
  }
  UseMethod("nlmixr2Est")
}

#' Show all the current estimation methods
#'
#' @return List of supported nlmixr2 estimation options (est=...)
#' @examples
#' nlmixr2AllEst()
#' @export
nlmixr2AllEst <- function() {
  .ret <- vapply(as.character(utils::methods("nlmixr2Est")), function(est){substr(est,12,nchar(est))}, character(1), USE.NAMES=FALSE)
  .ret[!(.ret %fin% c("default", "output"))]
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.default <- function(env, ...) {
  .curEst <- class(env)[1]
  stop("nlmixr2 estimation '", .curEst, "' not supported\n can be one of '", paste(nlmixr2AllEst(), collapse="', '"), "'",
       call.=FALSE)
}

#' For models with zero omega values or fixed values, update the
#' original model
#'
#' @param ret input for updating
#' @return nothing, called for side effects
#' @noRd
#' @author Matthew L. Fidler
.nlmixrEstUpdatesOrigModel <- function(ret) {
  .ui <- try(ret$ui, silent=TRUE)
  if (inherits(.ui, "try-error")) return(ret)
  if (inherits(.ui, "rxUi")) {
    # this needs to be in reverse order of the changes, which means apply zero omegas then fixed
    if (!is.null(nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi)) {
      .final <- nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi
      .finalIni <- .final$iniDf
      .iniDf <- .ui$iniDf
      .theta <- .iniDf[is.na(.iniDf$neta1), ]
      for (.i in seq_along(.theta$name)) {
        .w <- which(.finalIni$name == .theta$name[.i])
        if (length(.w) == 1) {
          .finalIni$est[.w] <- .theta$est[.i]
        }
      }
      .etas <- .iniDf[!is.na(.iniDf$neta1),, drop = FALSE]
      if (length(.etas$name) > 0) {
        .etaNames <- .etas[.etas$neta1 == .etas$neta2, "name"]
        .etaFinal <- vapply(.etaNames, function(n) {
          .finalIni[which(.finalIni$name == n), "neta1"]
        }, double(1), USE.NAMES=TRUE)

        .etaCur <- vapply(.etaNames, function(n) {
          .etas[which(.etas$name == n), "neta1"]
        }, double(1), USE.NAMES=TRUE)
        for (.i in seq_along(.etas$name)) {
          .eta1 <- .etaFinal[names(.etaCur)[which(.etas$neta1[.i] == .etaCur)]]
          .eta2 <- .etaFinal[names(.etaCur)[which(.etas$neta2[.i] == .etaCur)]]
          .finalIni$est[which(.finalIni$neta1 == .eta1 & .finalIni$neta2 == .eta2)] <- .etas$est[.i]
        }
      }
      assign("iniDf", .finalIni, .final)
      assign("ui", .final, envir=ret$env)
      assign("omega", .final$omega, envir=ret$env)
      .minfo("initial model updated with final estimates, some zero etas are excluded from output")
    }
    if (!is.null(nlmixr2global$nlmixr2EstEnv$uiUnfix)) {
      # Adjust to original model without literal fix
      .final <- nlmixr2global$nlmixr2EstEnv$uiUnfix
      .iniDf0 <- ret$env$ui$iniDf
      .iniDf2 <- .final$iniDf
      .iniDf2$est <- vapply(.iniDf2$name,
                            function(n) {
                              .w <- which(.iniDf0$name == n)
                              if (length(.w) == 1L) return(.iniDf0$est[.w])
                              .iniDf2[.iniDf2$name == n, "est"]
                            }, double(1), USE.NAMES = FALSE)
      assign("iniDf", .iniDf2, envir=.final)
      assign("ui", .final, envir=ret$env)
      assign("fixef", .final$theta, envir=ret$env)
    }
  }
  nlmixr2global$nlmixr2EstEnv$uiUnfix <- NULL
  nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi <- NULL
}

.tablePassthrough <- c("addDosing", "subsetNonmem", "cores", "keep", "drop")

#' Call nlmixr2Est wrapped to collect the warnings
#'
#'
#' @param env nlmixr2 estimate call
#' @param ... Other parameters
#' @return nlmixr2 object
#' @author Matthew L. Fidler
#' @noRd
nlmixr2Est0 <- function(env, ...) {
  rxode2::rxUnloadAll()
  .ui <- rxode2::rxUiDecompress(get("ui", env))
  assign("ui", .ui, envir=env)
  if (!exists("missingTable", envir=env)) {
    assign("missingTable", FALSE, envir=env)
  }
  if (!exists("missingControl", envir=env)) {
    assign("missingControl", FALSE, envir=env)
  }
  if (!exists("missingEst", envir=env)) {
    assign("missingEst", FALSE, envir=env)
  }
  if (!exists("missingTable", envir=env)) {
    assign("missingTable", TRUE, envir=env)
  }
  .doIt <- TRUE
  if (is.null(get("missingTable", envir=env))) {
  } else if (get("missingTable", envir=env)) {
  } else {
    .doIt <- FALSE
  }
  if (.doIt) {
    .meta <- get("ui", envir=env)$meta
    if (is.null(get("table", envir=env))) {
      assign("table", tableControl(), envir=env)
    }
    if (!is.environment(.meta)) {
      .meta <- new.env(parent=emptyenv())
    }
    .table <- get("table", envir=env)
    for (.elt in .tablePassthrough) {
      if (exists(.elt, envir=.meta)) {
        .table[[.elt]] <- .meta[[.elt]]
      }
    }
    assign("table", .table, envir=env)
  }
  .envReset <- new.env(parent=emptyenv())
  .envReset$reset <- TRUE
  if (!getOption("nlmixr2.resetCache", TRUE)) {
    .envReset$ret <- .collectWarn(nlmixr2Est(env, ...), lst = TRUE)
  } else {
    .envReset$reset <- TRUE
    .envReset$env <- new.env(parent=emptyenv())
    lapply(ls(envir = env, all.names = TRUE), function(item) {
      assign(item, get(item, envir = env), envir = .envReset$env)
    })
    .envReset$cacheReset <- FALSE
    .envReset$unload <- FALSE
    class(.envReset) <- class(env)
    if (length(get("reset", envir=.envReset)) != 1) assign("reset", TRUE, envir=.envReset)
    while (get("reset", envir=.envReset)) {
      assign("reset", FALSE, envir=.envReset)
      assign("ret", try(.collectWarn(nlmixr2Est(env, ...), lst = TRUE)), envir=.envReset)
      if (inherits(get("ret", envir=.envReset), "try-error")) {
        .msg <- attr(get("ret", envir=.envReset), "condition")$message
        if (regexpr("not provided by package", .msg) != -1) {
          if (get("cacheReset", envir=.envReset)) {
            .malert("unsuccessful cache reset; try manual reset with 'rxode2::rxClean()'")
            stop(.msg, call.=FALSE)
          } else {
            # reset
            if (is.environment(.envReset)) {
              rm(list=ls(envir = env, all.names = TRUE), envir=env)
              lapply(ls(envir = .envReset, all.names = TRUE), function(item) {
                assign(item, get(item, envir = .envReset), envir = env)
              })
            } else if (is.environment(.envReset$env)) {
              rm(list=ls(envir = env, all.names = TRUE), envir=env)
              lapply(ls(envir = .envReset$env, all.names = TRUE), function(item) {
                assign(item, get(item, envir = .envReset$env), envir = env)
              })
            }

            gc()
            .minfo("try resetting cache")
            rxode2::rxClean()
            assign("cacheReset", TRUE, envir=.envReset)
            assign("reset", TRUE, envir=.envReset)
            .msuccess("done")
          }
        } else if (regexpr("maximal number of DLLs reached", .msg) != -1) {
          if (.envReset$unload) {
            .malert("Could not unload rxode2 models, try restarting R")
            stop(.msg, call.=FALSE)
          } else {
            # reset
            if (is.environment(.envReset)) {
              rm(list=ls(envir = env, all.names = TRUE), envir=env)
              lapply(ls(envir = .envReset, all.names = TRUE), function(item) {
                assign(item, get(item, envir = .envReset), envir = env)
              })
            } else if (is.environment(.envReset$env)) {
              rm(list=ls(envir = env, all.names = TRUE), envir=env)
              lapply(ls(envir = .envReset$env, all.names = TRUE), function(item) {
                assign(item, get(item, envir = .envReset$env), envir = env)
              })
            }
            gc()
            .minfo("try resetting cache and unloading all rxode2 models")
            rxode2::rxUnloadAll(TRUE) # make sure this is actually unloading models
            try(rxode2::rxUnloadAll())
            rxode2::rxClean()
            assign("unload", TRUE, envir=.envReset)
            assign("reset", TRUE, envir=.envReset)
            .msuccess("done")
          }
        } else {
          stop(.msg, call.=FALSE)
        }
      }
      if (length(get("reset", envir=.envReset)) != 1) assign("reset", TRUE, .envReset)
    }
  }
  .lst <- get("ret", envir=.envReset)
  .ret <- .lst[[1]]
  if (inherits(.ret, "nlmixr2FitCore") ||
        inherits(.ret, "nlmixr2Fit")) {
    if (is.environment(.ret)) {
      try(assign("runInfo", .lst[[2]], .ret), silent=TRUE)
    } else {
      try(assign("runInfo", .lst[[2]], .ret$env), silent=TRUE)
    }
  } else {
    .w <-.lst[[2]]
    lapply(seq_along(.w), function(i) {
      warning(.w[[i]])
    })
  }
  .nlmixrEstUpdatesOrigModel(.ret)
  .ret
}
