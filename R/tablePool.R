#' Set up the post-fit table models in the ODE-swap pool
#'
#' Builds one solve pool over the fit's translated events, sized for the larger
#' of `inner` and `predOnly`, and registers both as pool slots (#1356).
#'
#' @param fit focei style fit
#' @param inner the table's ipred/pred model (`fit$innerModel`), or `NULL`
#' @param predOnly the display-column model (`fit$predOnlyModel`), or `NULL`
#' @param pars parameters the pool is built with (one row per subject)
#' @param method rxode2 ODE method name
#' @param atol,rtol ODE tolerances
#' @return named logical: whether each model can be solved in the pool
#' @author Matthew L. Fidler
#' @noRd
.tablePoolSetup <- function(fit, inner, predOnly, pars, method, atol, rtol) {
  # Registering binds entry points by name, so each dll must be loaded: one
  # unloaded with an earlier model object leaves stale entry points behind.
  # Loaded BEFORE the pool is built, since a load rebinds rxode2's model
  # globals to that model.
  for (.m in list(predOnly, inner)) {
    if (!is.null(.m)) rxode2::rxLoad(.m)
  }
  .pool <- tablePoolDeclare_(inner, predOnly)
  if (is.null(.pool)) {
    return(c(inner = FALSE, predOnly = FALSE))
  }
  .covsi <- .residCovsInterpolation(fit)
  .rxSolveNonmem(
    .pool,
    pars,
    fit$dataSav,
    atol = atol,
    rtol = rtol,
    maxsteps = fit$maxstepsOde,
    hmin = fit$hmin,
    hmax = fit$hmax,
    hini = fit$hini,
    maxordn = fit$maxordn,
    maxords = fit$maxords,
    method = rxode2::odeMethodToInt(method),
    tolFactor = fit$env$tolFactor,
    covsInterpolation = .covsi$covsInterpolation,
    naInterpolation = .covsi$naInterpolation,
    nonmem = .covsi$nonmem,
    warnIdSort = FALSE,
    .setupOnly = 1L
  )
  tablePoolRegister_(inner, predOnly)
}

#' Solve one table model in the pool
#'
#' @param model the model, which must be the one registered for `slot`
#' @param slot `"inner"` or `"predOnly"`
#' @param pars named parameter vector (all subjects) or matrix (one row per
#'   subject)
#' @param lhs,state lhs/state names to output as columns
#' @param sens lhs names to return as sparse triplets instead
#' @param cov model parameters (covariates) to output, as interpolated
#' @param tbs output `rxLambda`/`rxYj`/`rxLow`/`rxHi`
#' @return list: `id` (subject index), `rec` (0-based record within the
#'   subject), `evid`, `time`, `vals` (data frame of the requested columns),
#'   `sens` (triplets) and `rowsOk`
#' @author Matthew L. Fidler
#' @noRd
.tablePoolSolve <- function(
  model,
  slot,
  pars,
  lhs = character(0),
  state = character(0),
  sens = character(0),
  cov = character(0),
  tbs = FALSE
) {
  .mv <- rxode2::rxModelVars(model)
  if (is.null(dim(pars))) {
    pars <- matrix(pars, nrow = 1L, dimnames = list(NULL, names(pars)))
  } else {
    pars <- as.matrix(pars)
  }
  .idx <- function(x, from) {
    .i <- match(x, from)
    if (anyNA(.i)) {
      stop("not in the table model: ", paste(x[is.na(.i)], collapse = ", "), call. = FALSE)
    }
    as.integer(.i - 1L)
  }
  .ret <- tablePoolSolve_(
    switch(slot, inner = 0L, predOnly = 1L),
    pars,
    .idx(colnames(pars), .mv$params),
    .idx(lhs, .mv$lhs),
    .idx(sens, .mv$lhs),
    .idx(state, .mv$state),
    .idx(cov, .mv$params),
    tbs
  )
  .vals <- as.data.frame(.ret$vals)
  names(.vals) <- c(lhs, state, cov, if (tbs) c("rxLambda", "rxYj", "rxLow", "rxHi"))
  .ret$vals <- .vals
  .ret
}

# ODE methods the pool solves; any other method takes the rxSolve() route
.tablePoolMethods <- c("liblsoda", "lsoda", "dop853")

.tablePoolEnv <- new.env(parent = emptyenv())
.tablePoolEnv$declineN <- 0L
.tablePoolEnv$rebuildN <- 0L
.tablePoolEnv$pooledN <- 0L
.tablePoolEnv$lastDecline <- NA_character_

#' Counters of the pooled table route
#'
#' @return list: `pooledN` (pooled table predictions), `declineN` and
#'   `lastDecline` (fell back to rxSolve()), `rebuildN` (pool rebuilt for an
#'   ODE method or tolerance change)
#' @noRd
.tablePoolInfo <- function() {
  as.list(.tablePoolEnv)
}

.tablePoolDecline <- function(why) {
  .tablePoolEnv$declineN <- .tablePoolEnv$declineN + 1L
  .tablePoolEnv$lastDecline <- why
  NULL
}

.tablePoolStop <- function(why) {
  stop(structure(class = c("tablePoolDecline", "error", "condition"), list(message = why, call = NULL)))
}

#' Why the pooled route cannot make this fit's table, or `NULL`
#'
#' @noRd
.tablePoolWhyNot <- function(fit, models, keep, addDosing) {
  if (!isFALSE(addDosing)) {
    return("addDosing")
  }
  .env <- fit
  if (!is.environment(.env) && is.environment(fit$env)) {
    .env <- fit$env
  }
  if (is.environment(.env) && exists("mixIcov", envir = .env, inherits = FALSE)) {
    return("mixture")
  }
  if (!(.residOdeMethodName(fit$methodOde) %in% .tablePoolMethods)) {
    return("method")
  }
  if (!is.null(fit$foceiControl$rxControl$scale)) {
    return("scale")
  }
  .data <- fit$dataSav
  if (!all(keep %in% names(.data))) {
    return("keep")
  }
  .evid <- .data[[which(toupper(names(.data)) == "EVID")[1]]]
  if (is.null(.evid) || any(.evid == 3L | (.evid >= 10L & .evid <= 99L))) {
    return("evid")
  }
  for (.m in models) {
    .mv <- rxode2::rxModelVars(.m)
    if (isTRUE(.mv$nMtime > 0L) || any(.mv$lhsStr != 0L) || length(.mv$strAssign) > 0L) {
      return("model")
    }
  }
  NULL
}

#' The ipred/pred list of `.foceiPredIpredList()`, solved in the ODE-swap pool
#'
#' Returns `NULL` (and counts it) when the pooled route does not cover the fit,
#' and the caller solves with rxSolve() instead.  The eta sensitivities are not
#' columns of the frames but their sparse `rxSens` attribute, read by
#' `cwresCalc`.
#'
#' @inheritParams .foceiPredIpredList
#' @param keep columns to carry from the data
#' @return list as `.foceiPredIpredList()`, or `NULL`
#' @author Matthew L. Fidler
#' @noRd
.tablePoolPredIpredList <- function(fit, thetaEtaParameters, keep, predOnly, addDosing) {
  if (!isTRUE(getOption("nlmixr2.tablePool", TRUE))) {
    return(NULL)
  }
  .po <- fit$predOnlyModel
  .inner <- if (predOnly) NULL else fit$innerModel
  if (is.null(.po)) {
    return(.tablePoolDecline("model"))
  }
  .why <- .tablePoolWhyNot(fit, c(list(.po), if (!is.null(.inner)) list(.inner)), keep, addDosing)
  if (!is.null(.why)) {
    return(.tablePoolDecline(.why))
  }
  .ret <- tryCatch(
    .tablePoolPredIpredList0(fit, thetaEtaParameters, keep, .inner, .po),
    tablePoolDecline = function(e) {
      .tablePoolDecline(conditionMessage(e))
    }
  )
  if (!is.null(.ret)) {
    .tablePoolEnv$pooledN <- .tablePoolEnv$pooledN + 1L
  }
  .ret
}

#' Seed state of rxode2 and of R, to put back after the pooled solves
#'
#' @noRd
.tablePoolSeedSave <- function() {
  list(rx = rxode2::rxGetSeed(), r = get0(".Random.seed", envir = globalenv(), inherits = FALSE))
}

#' Leave the seed where the rxSolve() route would have
#'
#' Each rxSolve() advances rxode2's running seed by its subjects and cores (or,
#' unseeded, draws one runif()), and the pooled solves take none.  Later
#' draws (the censored-DV imputation, a following simulation) must not depend
#' on which route made the table, so restore the state and advance it by
#' `nsolve` rxSolve() calls.
#'
#' @param st `.tablePoolSeedSave()` from before the pooled solves
#' @param nsolve number of rxSolve() calls the other route makes
#' @param step seed advance of one of them: subjects plus solve cores
#' @noRd
.tablePoolSeedAdvance <- function(st, nsolve, step) {
  if (st$rx != -1) {
    rxode2::rxSetSeed(st$rx + nsolve * step)
    return(invisible())
  }
  if (is.null(st$r)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", st$r, envir = globalenv())
  }
  for (.i in seq_len(nsolve)) {
    stats::runif(1L, 1.0, 4294967295)
  }
  invisible()
}

.tablePoolPredIpredList0 <- function(fit, thetaEtaParameters, keep, inner, po) {
  on.exit(tablePoolFree_(), add = TRUE)
  .seed <- .tablePoolSeedSave()
  # solve count and per-solve seed step, updated by the closures below
  .acc <- new.env(parent = emptyenv())
  .acc$nsolve <- 0L
  .acc$step <- 0L
  .done <- FALSE
  on.exit(
    .tablePoolSeedAdvance(.seed, if (.done) .acc$nsolve else 0L, .acc$step),
    add = TRUE
  )
  .data <- fit$dataSav
  .ipar <- thetaEtaParameters$ipred
  .ppar <- thetaEtaParameters$pred
  .need <- c(inner = !is.null(inner), predOnly = TRUE)
  on.exit(
    {
      rxode2::rxEventSensDeactivate()
      odeSwapEsNoteInstalled_(-1L)
    },
    add = TRUE
  )
  .st <- new.env(parent = emptyenv())
  .st$key <- NULL
  .ensure <- function(method, atol, rtol) {
    .key <- list(method, atol, rtol)
    if (identical(.st$key, .key)) {
      return(invisible())
    }
    if (!(method %in% .tablePoolMethods)) {
      .tablePoolStop("method")
    }
    if (!is.null(.st$key)) {
      .tablePoolEnv$rebuildN <- .tablePoolEnv$rebuildN + 1L
    }
    .st$key <- NULL
    # a set-up-only rxSolve() advances the seed by the cores a full solve
    # takes on top of its subjects
    .s0 <- rxode2::rxGetSeed()
    .ok <- .tablePoolSetup(fit, inner, po, .ipar, method, atol, rtol)
    if (.s0 != -1) {
      .acc$step <- length(.start) + (rxode2::rxGetSeed() - .s0)
    }
    if (!all(.ok[.need])) {
      .tablePoolStop("register")
    }
    # the inner model's event ("jump") sensitivity shape, installed after the
    # loads above (loading a model rebinds it), as .foceiFitInternal() does
    if (
      !is.null(inner) &&
        isTRUE(tryCatch(rxode2::rxEventSensLoadModel(inner), error = function(e) FALSE))
    ) {
      odeSwapEsNoteInstalled_(1L)
    } else {
      odeSwapEsNoteInstalled_(-1L)
    }
    .st$key <- .key
    invisible()
  }
  .solve <- function(what, model, slot, pars, ...) {
    .residSolveWithFallback(
      fit,
      what,
      function(method, atol, rtol) {
        .acc$nsolve <- .acc$nsolve + 1L
        .ensure(method, atol, rtol)
        .r <- .tablePoolSolve(model, slot, pars, ...)
        if (!isTRUE(.r$rowsOk)) {
          .tablePoolStop("rows")
        }
        .r
      },
      predOf = function(x) x$vals$rx_pred_
    )
  }
  # each subject's first row in the data: a solve's records map back onto it
  .start <- match(unique(.data$ID), .data$ID) - 1L
  .nsub <- length(.start)
  .hasEvid2 <- any(.data[[which(toupper(names(.data)) == "EVID")[1]]] == 2L)
  .frame <- function(r, keep = character(0), tbs = FALSE) {
    .df <- data.frame(id = factor(r$id, levels = seq_len(.nsub)))
    if (.hasEvid2) {
      .df$evid <- r$evid
    }
    .df$time <- r$time
    .v <- r$vals
    .tbs <- c("rxLambda", "rxYj", "rxLow", "rxHi")
    .df <- cbind(.df, .v[, setdiff(names(.v), .tbs), drop = FALSE])
    # as in rxSolve(), a kept column that is already output is not repeated
    keep <- setdiff(keep, names(.df))
    if (length(keep) > 0L) {
      .rows <- .start[r$id] + r$rec + 1L
      .df <- cbind(.df, .data[.rows, keep, drop = FALSE])
    }
    if (tbs) {
      .df <- cbind(.df, .v[, .tbs, drop = FALSE])
    }
    rownames(.df) <- NULL
    if (length(r$sens$x) > 0L || r$sens$ncol > 0L) {
      attr(.df, "rxSens") <- r$sens
    }
    .df
  }
  .mvP <- rxode2::rxModelVars(po)
  .covP <- setdiff(intersect(.mvP$params, names(.data)), c(.mvP$state, .mvP$lhs))
  if (is.null(inner)) {
    .ipred <- .solve(
      "ipred",
      po,
      "predOnly",
      .ipar,
      lhs = .mvP$lhs,
      state = .mvP$state,
      cov = .covP,
      tbs = TRUE
    )
    .pred <- .solve("pred", po, "predOnly", .ppar, lhs = c("rx_pred_", "rx_r_"))
    .done <- TRUE
    return(list(
      ipred = .residAdjustIpredNames(.frame(.ipred, keep, tbs = TRUE)),
      pred = .frame(.pred),
      etaLst = thetaEtaParameters$eta.lst
    ))
  }
  .mvI <- rxode2::rxModelVars(inner)
  .neta <- ncol(fit$omega)
  .sens <- paste0("rx__sens_rx_pred__BY_ETA_", seq_len(.neta), "___")
  if (!all(c("rx_pred_", "rx_r_", .sens) %in% .mvI$lhs)) {
    .tablePoolStop("model")
  }
  .stI <- intersect(.mvI$state, c(.mvP$state, .mvP$stateExtra))
  .ipred <- .solve(
    "ipred",
    inner,
    "inner",
    .ipar,
    lhs = c("rx_pred_", "rx_r_"),
    state = .stI,
    sens = .sens,
    tbs = TRUE
  )
  .pred <- .solve("pred", inner, "inner", .ppar, lhs = c("rx_pred_", "rx_r_"), sens = .sens)
  .ebe <- .solve("ebe", po, "predOnly", .ipar, lhs = .mvP$lhs, state = .mvP$state, cov = .covP)
  .done <- TRUE
  list(
    ipred = .residAdjustIpredNames(.frame(.ipred, keep, tbs = TRUE)),
    pred = .frame(.pred),
    etaLst = thetaEtaParameters$eta.lst,
    predOnly = .frame(.ebe, keep)
  )
}
