#' Likelihood penalty specification
#'
#' Describes one compiled likelihood penalty (see
#' `vignette("likelihood-penalties", package = "nlmixr2est")`) to apply
#' during a single fit.  A model extension attaches specs to a model with
#' [nlmixr2LikPenaltyAdd()]; they can also be passed to `nlmixr2()` as
#' `likPenalty=`.
#'
#' @param ptr the penalty: an external pointer to a C `nlmixrLikPenalty`, a
#'   function returning one (resolved when the fit starts, so the spec
#'   survives `saveRDS()`), or the name of a built-in penalty (`"thetaL2"`,
#'   `"etaL2"`).
#' @param thetas names of the population parameters the penalty acts on;
#'   handed to the C callbacks as 0-based theta indices (`ctx->sel`).
#' @param par numeric settings handed to the C callbacks (`ctx->par`).
#' @param data optional external pointer handed to the C callbacks
#'   (`ctx->data`); defaults to the bundle's own `data`.
#' @param inObjf should the penalty be counted in the reported objective
#'   function?  `NA` uses the bundle's `NLMIXR_PEN_IN_OBJF` flag.
#' @param name label used in `fit$likPenalty`; defaults to the bundle's name.
#' @return a `nlmixr2LikPenaltySpec` object
#' @author Matthew L. Fidler
#' @export
#' @examples
#' # ridge penalty of 10*tka^2 (on the -2LL scale), shaping the fit only
#' nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = 10)
nlmixr2LikPenaltySpec <- function(ptr, thetas = character(0), par = numeric(0), data = NULL, inObjf = NA, name = NULL) {
  if (!(inherits(ptr, "externalptr") || is.function(ptr) || checkmate::testString(ptr))) {
    stop("'ptr' must be an external pointer, a function or a built-in name", call. = FALSE)
  }
  checkmate::assertCharacter(thetas, any.missing = FALSE)
  checkmate::assertNumeric(par, any.missing = FALSE)
  checkmate::assertLogical(inObjf, len = 1)
  if (!is.null(data) && !inherits(data, "externalptr")) {
    stop("'data' must be NULL or an external pointer", call. = FALSE)
  }
  if (is.null(name)) {
    name <- if (is.character(ptr)) ptr else NA_character_
  }
  checkmate::assertString(name, na.ok = TRUE)
  structure(
    list(
      ptr = ptr,
      thetas = thetas,
      par = as.double(par),
      data = data,
      inObjf = inObjf,
      name = name
    ),
    class = "nlmixr2LikPenaltySpec"
  )
}

#' @export
print.nlmixr2LikPenaltySpec <- function(x, ...) {
  .nm <- if (is.na(x$name)) "<compiled>" else x$name
  cat("likelihood penalty", .nm)
  if (length(x$thetas) > 0) {
    cat(" on", paste(x$thetas, collapse = ", "))
  }
  if (length(x$par) > 0) {
    cat(" par:", paste(x$par, collapse = ", "))
  }
  cat(" inObjf:", x$inObjf, "\n")
  invisible(x)
}

#' Attach likelihood penalties to a model
#'
#' Used by model extensions (for example from an rxode2 user-function UI
#' handler) so the penalty is active only while this model is fit.
#'
#' @param ui rxode2 model (or function)
#' @param ... `nlmixr2LikPenaltySpec()` objects
#' @param replace replace the penalties already attached instead of adding to
#'   them
#' @return the rxode2 ui with the penalties stored in `ui$nlmixr2LikPenalty`
#' @author Matthew L. Fidler
#' @export
#' @examples
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45; tcl <- 1; tv <- 3.45
#'     eta.ka ~ 0.6
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl)
#'     v <- exp(tv)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#' ui <- nlmixr2LikPenaltyAdd(one.cmt,
#'   nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = 10))
#' ui$nlmixr2LikPenalty
nlmixr2LikPenaltyAdd <- function(ui, ..., replace = FALSE) {
  .specs <- list(...)
  for (.s in .specs) {
    if (!inherits(.s, "nlmixr2LikPenaltySpec")) {
      stop("'...' must be nlmixr2LikPenaltySpec() objects", call. = FALSE)
    }
  }
  .ui <- rxode2::rxUiDecompress(rxode2::rxode2(ui))
  .meta <- get("meta", envir = .ui)
  .old <- if (replace) NULL else .likPenaltyAsList(.ui$nlmixr2LikPenalty)
  assign("nlmixr2LikPenalty", c(.old, .specs), envir = .meta)
  .ui
}

#' Register a likelihood penalty from an rxode2 user-function expansion
#'
#' Called from an `rxUdfUi` method while it expands a model term (for
#' example a smoothing spline) against the data.  The method must return
#' `uiUseData = TRUE` from its data-free stage, so `nlmixr2()` re-parses the
#' model with the estimation data (`rxode2::rxUdfUiData()`); the specs
#' registered during that re-parse are attached to the expanded model with
#' [nlmixr2LikPenaltyAdd()], so the penalty is active only while that model
#' is estimated.  Outside of an nlmixr2 data expansion (for example in
#' `rxSolve()`) the call does nothing.  See the "A model term that brings its
#' own penalty" section of `vignette("likelihood-penalties")` and the
#' 'nlmixr2spline' package for a complete example.
#'
#' @param spec a [nlmixr2LikPenaltySpec()] with a unique `name`; a later
#'   registration with the same name replaces it
#' @return `TRUE` when the penalty was registered, invisibly
#' @author Matthew L. Fidler
#' @export
nlmixr2LikPenaltyUdf <- function(spec) {
  if (!inherits(spec, "nlmixr2LikPenaltySpec")) {
    stop("'spec' must be a nlmixr2LikPenaltySpec() object", call. = FALSE)
  }
  if (is.na(spec$name)) {
    stop("a user-function likelihood penalty needs a 'name'", call. = FALSE)
  }
  .cur <- nlmixr2global$likPenaltyUdf
  if (!is.list(.cur)) {
    return(invisible(FALSE))
  }
  .cur[[spec$name]] <- spec
  nlmixr2global$likPenaltyUdf <- .cur
  invisible(TRUE)
}

#' Attach the penalties registered while expanding a model
#' @param ui expanded ui
#' @param specs named list of specs from nlmixr2LikPenaltyUdf()
#' @return ui with the specs attached (same-named ones replaced)
#' @noRd
.likPenaltyUdfAttach <- function(ui, specs) {
  if (length(specs) == 0L) {
    return(ui)
  }
  .old <- .likPenaltyAsList(ui$nlmixr2LikPenalty)
  .keep <- Filter(function(s) !(s$name %in% names(specs)), .old)
  do.call(nlmixr2LikPenaltyAdd, c(list(ui), .keep, unname(specs), list(replace = TRUE)))
}

.likPenaltyAsList <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "nlmixr2LikPenaltySpec")) {
    return(list(x))
  }
  if (!is.list(x) || !all(vapply(x, inherits, logical(1), "nlmixr2LikPenaltySpec"))) {
    stop("likelihood penalties must be nlmixr2LikPenaltySpec() objects", call. = FALSE)
  }
  x
}

#' Penalty specs requested for this estimation (model + `likPenalty=`)
#' @param env estimation environment
#' @return list of specs
#' @noRd
.likPenaltySpecs <- function(env) {
  .ui <- get("ui", envir = env)
  .fromUi <- tryCatch(.ui$nlmixr2LikPenalty, error = function(e) NULL)
  .dots <- if (exists(".nlmixr2Dots", envir = env, inherits = FALSE)) {
    env$.nlmixr2Dots$likPenalty
  }
  c(.likPenaltyAsList(.fromUi), .likPenaltyAsList(.dots))
}

.likPenaltyResolvePtr <- function(ptr) {
  if (is.character(ptr)) {
    return(.Call(`_nlmixr2est_likPenaltyBuiltin`, ptr))
  }
  if (is.function(ptr)) {
    ptr <- ptr()
  }
  if (!inherits(ptr, "externalptr")) {
    stop("the likelihood penalty function must return an external pointer", call. = FALSE)
  }
  ptr
}

.likPenaltyLevels <- c("none", "ignore", "subject", "full")

#' Which likelihood penalties does the dispatched method support?
#' @param env estimation environment
#' @return `"none"` (refused), `"ignore"` (no estimation, e.g. simulation),
#'   `"subject"` (subject terms only; population terms are added to the
#'   reported objective) or `"full"`
#' @noRd
.likPenaltySupport <- function(env) {
  for (.cls in class(env)) {
    .fn <- utils::getS3method("nlmixr2Est", .cls, optional = TRUE)
    if (is.null(.fn)) {
      next
    }
    .a <- attr(.fn, "nlmixr2LikPenalty")
    if (is.null(.a)) {
      return("none")
    }
    return(match.arg(.a, .likPenaltyLevels))
  }
  "none"
}

#' Activate this fit's likelihood penalties
#'
#' Only the outermost estimation activates; nested estimations (a saem fit's
#' objective, setOfv(), ...) run under the same penalties.
#'
#' @param env estimation environment
#' @return TRUE when this call activated penalties (and must deactivate)
#' @noRd
.likPenaltyActivate <- function(env) {
  if (isTRUE(nlmixr2global$likPenaltyActive)) {
    return(FALSE)
  }
  .specs <- .likPenaltySpecs(env)
  if (length(.specs) == 0L) {
    return(FALSE)
  }
  .level <- .likPenaltySupport(env)
  if (.level == "none") {
    stop("this estimation method does not support likelihood penalties", call. = FALSE)
  }
  if (.level == "ignore") {
    return(FALSE)
  }
  .iniDf <- get("ui", envir = env)$iniDf
  .thetaNames <- .iniDf$name[!is.na(.iniDf$ntheta)][order(.iniDf$ntheta[!is.na(.iniDf$ntheta)])]
  .sel <- lapply(.specs, function(s) {
    .w <- match(s$thetas, .thetaNames)
    if (anyNA(.w)) {
      stop("likelihood penalty theta not in model: ", paste(s$thetas[is.na(.w)], collapse = ", "), call. = FALSE)
    }
    as.integer(.w - 1L)
  })
  .ptrs <- lapply(.specs, function(s) .likPenaltyResolvePtr(s$ptr))
  .Call(
    `_nlmixr2est_likPenaltyActivate`,
    .ptrs,
    .sel,
    lapply(.specs, function(s) s$par),
    lapply(.specs, function(s) s$data),
    vapply(.specs, function(s) if (is.na(s$inObjf)) NA_integer_ else as.integer(s$inObjf), integer(1))
  )
  nlmixr2global$likPenaltyActive <- TRUE
  nlmixr2global$likPenaltySpecs <- .specs
  nlmixr2global$likPenaltyLevel <- .level
  TRUE
}

.likPenaltyDeactivate <- function() {
  .Call(`_nlmixr2est_likPenaltyDeactivate`)
  nlmixr2global$likPenaltyActive <- NULL
  nlmixr2global$likPenaltySpecs <- NULL
  nlmixr2global$likPenaltyLevel <- NULL
  nlmixr2global$likPenaltySaemMap <- NULL
  nlmixr2global$likPenaltySaemTheta0 <- NULL
  invisible()
}

#' Are likelihood penalties active right now?
#' @return integer vector: total, subject, population, subject-changed
#' @noRd
.likPenaltyInfo <- function() {
  .r <- .Call(`_nlmixr2est_likPenaltyInfo`)
  names(.r) <- c("n", "subject", "population", "subjectChanged")
  .r
}

#' Evaluate the active penalties at a fit's estimates
#' @param fit nlmixr2 fit
#' @return data frame with one row per penalty
#' @noRd
.likPenaltyEvalFit <- function(fit) {
  .iniDf <- fit$iniDf
  .th <- .iniDf[!is.na(.iniDf$ntheta), ]
  .th <- .th[order(.th$ntheta), ]
  .theta <- fit$theta
  .thetaV <- setNames(.th$est, .th$name)
  if (!is.null(.theta)) {
    .w <- intersect(names(.theta), names(.thetaV))
    .thetaV[.w] <- .theta[.w]
  }
  .eta <- fit$eta
  .etaNames <- .iniDf$name[!is.na(.iniDf$neta1) & .iniDf$neta1 == .iniDf$neta2]
  .etaM <- matrix(numeric(0), 0, 0)
  .nobs <- integer(0)
  if (is.data.frame(.eta) && nrow(.eta) > 0) {
    .etaM <- as.matrix(.eta[, intersect(.etaNames, names(.eta)), drop = FALSE])
    .ids <- .eta$ID
  } else {
    .ids <- unique(fit$ID)
  }
  .tab <- table(factor(fit$ID, levels = unique(.ids)))
  .nobs <- as.integer(.tab)
  if (nrow(.etaM) == 0L) {
    .etaM <- matrix(numeric(0), length(.nobs), 0)
  }
  .om <- fit$omega
  if (!is.matrix(.om)) {
    .om <- matrix(numeric(0), 0, 0)
  }
  .r <- .Call(`_nlmixr2est_likPenaltyEval`, as.double(.thetaV), .etaM, .om, .nobs)
  .specs <- nlmixr2global$likPenaltySpecs
  .nm <- .r$name
  if (length(.specs) == length(.nm)) {
    .nm2 <- vapply(.specs, function(s) s$name, character(1))
    .nm[!is.na(.nm2)] <- .nm2[!is.na(.nm2)]
  }
  data.frame(
    name = .nm,
    inObjf = .r$inObjf,
    population = .r$population,
    subjectLL = .r$subjectLL,
    objf = .r$population - 2 * .r$subjectLL
  )
}

#' Record the penalties on a fit and remove optimization-only ones from its
#' reported objective
#' @param fit nlmixr2 fit
#' @return nothing, called for side effects
#' @noRd
.likPenaltyFinalize <- function(fit) {
  if (!isTRUE(nlmixr2global$likPenaltyActive)) {
    return(invisible())
  }
  if (!inherits(fit, "nlmixr2FitCore")) {
    return(invisible())
  }
  .env <- fit$env
  # a lazily computed objective (saem) must be computed while the penalties
  # are still on
  .objDf0 <- try(get("objDf", envir = .env), silent = TRUE)
  if (is.data.frame(.objDf0) && all(is.na(.objDf0$OBJF))) {
    try(fit$objf, silent = TRUE)
  }
  .df <- try(.likPenaltyEvalFit(fit), silent = TRUE)
  if (inherits(.df, "try-error")) {
    warning("could not evaluate the likelihood penalty at the estimates", call. = FALSE)
    return(invisible())
  }
  assign("likPenalty", .df, envir = .env)
  # the fit's objective holds every subject term, and every population term
  # unless the method only uses subject terms
  .sub <- -2 * .df$subjectLL
  .delta <- sum(.sub[!.df$inObjf])
  if (identical(nlmixr2global$likPenaltyLevel, "subject")) {
    .delta <- .delta - sum(.df$population[.df$inObjf])
    if (any(.df$population != 0)) {
      warning("population likelihood penalty added to objective after fit", call. = FALSE)
    }
  } else {
    .delta <- .delta + sum(.df$population[!.df$inObjf])
  }
  if (.delta == 0 || !is.finite(.delta)) {
    return(invisible())
  }
  .objDf <- get("objDf", envir = .env)
  if (is.data.frame(.objDf) && nrow(.objDf) > 0) {
    for (.c in c("OBJF", "AIC", "BIC")) {
      if (!is.null(.objDf[[.c]])) .objDf[[.c]] <- .objDf[[.c]] - .delta
    }
    if (!is.null(.objDf[["Log-likelihood"]])) {
      .objDf[["Log-likelihood"]] <- .objDf[["Log-likelihood"]] + .delta / 2
    }
    assign("objDf", .objDf, envir = .env)
  }
  for (.v in c("OBJF", "objf", "objective", "AIC", "BIC")) {
    if (exists(.v, envir = .env, inherits = FALSE)) {
      assign(.v, get(.v, envir = .env) - .delta, envir = .env)
    }
  }
  if (exists("logLik", envir = .env, inherits = FALSE)) {
    .ll <- get("logLik", envir = .env)
    .a <- attributes(.ll)
    .ll <- as.numeric(.ll) + .delta / 2
    attributes(.ll) <- .a
    assign("logLik", .ll, envir = .env)
  }
  if (any(!.df$inObjf & .df$objf != 0)) {
    warning("likelihood penalty excluded from objective (see $likPenalty)", call. = FALSE)
  }
  invisible()
}

#' Built-in likelihood penalty pointer
#'
#' @param name `"thetaL2"` (population: `par[1] * sum((theta - par[2])^2)`
#'   on the -2LL scale) or `"etaL2"` (subject: `par[2] - par[1]/2 *
#'   sum(eta^2)` added to each subject's log-likelihood)
#' @return external pointer to the C `nlmixrLikPenalty`
#' @export
#' @keywords internal
.nlmixr2LikPenaltyBuiltin <- function(name) {
  .Call(`_nlmixr2est_likPenaltyBuiltin`, name)
}

#' Add the likelihood-penalty theta map to a saem configuration
#'
#' @param a saem configuration list
#' @param estNames names of the saem fixed effects (`Plambda` order)
#' @param thetaNames ui theta names (theta order)
#' @param thetaEst ui theta estimates (theta order)
#' @return the configuration, with `penMap`/`penTheta0` when penalties are on
#' @noRd
.saemLikPenaltyCfg <- function(a, estNames, thetaNames, thetaEst) {
  if (.likPenaltyInfo()[["n"]] == 0L) {
    return(a)
  }
  if (length(estNames) != a$nlambda1 + a$nlambda0) {
    warning("saem fixed effects not matched to thetas; penalty not in M-step", call. = FALSE)
    return(a)
  }
  .m <- match(estNames, thetaNames) - 1L
  .m[is.na(.m)] <- -1L
  a$penMap <- as.integer(.m)
  a$penTheta0 <- as.double(thetaEst)
  # the fit's saem.cfg is rebuilt later; the objective reads these instead
  nlmixr2global$likPenaltySaemMap <- a$penMap
  nlmixr2global$likPenaltySaemTheta0 <- a$penTheta0
  a
}
