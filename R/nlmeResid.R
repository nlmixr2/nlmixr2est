#' Stop with an nlme residual-error refusal
#'
#' @param ... pieces of the message, after "nlme cannot estimate "
#' @return never returns
#' @noRd
.nlmeResidStop <- function(...) {
  stop("nlme cannot estimate ", ..., "; use est=\"focei\", \"saem\" or \"imp\"", call. = FALSE)
}

#' Is a call an assignment (`~` is a suppressed one)?
#'
#' @param e expression
#' @return logical
#' @noRd
.nlmeIsAssign <- function(e) {
  is.call(e) &&
    length(e) == 3L &&
    (identical(e[[1]], quote(`<-`)) || identical(e[[1]], quote(`=`)) || identical(e[[1]], quote(`~`)))
}

#' Names a statement assigns
#'
#' @param e expression
#' @return character vector of the names assigned anywhere in `e`
#' @noRd
.nlmeAssigned <- function(e) {
  if (!is.call(e)) {
    return(character(0))
  }
  if (.nlmeIsAssign(e)) {
    if (is.name(e[[2]])) {
      return(c(as.character(e[[2]]), .nlmeAssigned(e[[3]])))
    }
    return(.nlmeAssigned(e[[3]]))
  }
  unlist(lapply(as.list(e)[-1], .nlmeAssigned))
}

#' Is a statement a simple top-level assignment?
#'
#' @param e expression
#' @return the assigned name, or `NULL`
#' @noRd
.nlmeSimpleAssign <- function(e) {
  # a chained `a <- b <- x` is not simple
  if (.nlmeIsAssign(e) && is.name(e[[2]]) && !.nlmeIsAssign(e[[3]])) {
    return(as.character(e[[2]]))
  }
  NULL
}

#' Residual error structure of a model for nlme
#'
#' The residual SD of a single add/prop/pow endpoint may be an estimated
#' parameter times a function of time and covariates (intermediate variables
#' are resolved); anything else is refused with a specific message.
#'
#' @param ui rxode2 ui
#' @return `NULL` when every residual parameter is a plain theta, otherwise a
#'   list with the residual `slot`, its resolved `expr`, the estimated `theta`,
#'   the `fixed` theta values, the `covs` it uses and the residual-only
#'   thetas (`exclude`) to keep out of the structural parameters
#' @noRd
.nlmeResidInfo <- function(ui) {
  .predDf <- ui$predDf
  .slots <- intersect(c("a", "b", "c", "d", "e", "f", "lambda"), names(.predDf))
  .val <- vapply(.slots, function(s) as.character(.predDf[[s]][1]), character(1))
  .val <- .val[!is.na(.val)]
  if (length(.val) == 0L) {
    return(NULL)
  }
  .other <- setdiff(names(.val), c("a", "b"))
  if (length(.other) > 0L) {
    .nlmeResidStop("a power exponent or transform that depends on time, covariates or other parameters")
  }
  .errType <- as.character(.predDf$errType[1])
  if (!(.errType %in% c("add", "prop", "pow"))) {
    .nlmeResidStop("a combined residual error whose SD depends on time, covariates or other parameters")
  }
  .slot <- names(.val)[1]
  .iniDf <- ui$iniDf
  .thetaDf <- .iniDf[!is.na(.iniDf$ntheta), ]
  .etas <- .iniDf$name[!is.na(.iniDf$neta1)]
  .covs <- ui$allCovs
  .states <- ui$state
  # endpoint lines (`cp ~ add(...)`) assign nothing
  .lst <- ui$lstExpr
  .lst <- .lst[-.predDf$line]
  .defs <- list()
  .acc <- new.env(parent = emptyenv())
  .acc$cnt <- integer(0)
  .add <- function(v, n) {
    .acc$cnt[v] <- ifelse(is.na(.acc$cnt[v]), 0L, .acc$cnt[v]) + n
  }
  for (.e in .lst) {
    .v <- .nlmeSimpleAssign(.e)
    if (!is.null(.v)) {
      .defs[[.v]] <- .e[[3]]
      .add(.v, 1L)
      lapply(.nlmeAssigned(.e[[3]]), .add, n = 2L)
    } else {
      lapply(.nlmeAssigned(.e), .add, n = 2L)
    }
  }
  .resolve <- function(expr, chain = character(0)) {
    for (.i in 1:100) {
      .vars <- all.vars(expr)
      .sub <- list()
      for (.v in .vars) {
        if (.v %in% c(.thetaDf$name, .etas, .covs, "time", "t", "pi")) {
          next
        }
        if (.v %in% names(.defs) && !(.v %in% .states) && isTRUE(.acc$cnt[.v] == 1L)) {
          .sub[[.v]] <- call("(", .defs[[.v]])
          chain <- c(chain, .v)
        } else {
          .nlmeResidStop("a residual error that depends on the model-calculated '", .v, "'")
        }
      }
      if (length(.sub) == 0L) {
        return(list(expr = expr, chain = unique(chain)))
      }
      expr <- do.call(substitute, list(expr, .sub))
    }
    .nlmeResidStop("a residual error defined recursively")
  }
  .name <- .val[[1]]
  .r <- .resolve(str2lang(.name))
  .expr <- .r$expr
  .chain <- .r$chain
  # what the structural model uses; an intermediate it shares must be data-only
  .struct <- character(0)
  for (.e in .lst) {
    .v <- .nlmeSimpleAssign(.e)
    if (!is.null(.v) && .v %in% .chain) {
      next
    }
    .struct <- c(.struct, all.vars(.e))
  }
  # the endpoint's prediction and pow() exponent belong to the model too
  .struct <- c(.struct, as.character(.predDf$var), .iniDf$name[which(.iniDf$err == "pow2")])
  for (.v in intersect(.struct, .chain)) {
    if (any(all.vars(.resolve(.defs[[.v]])$expr) %in% c(.thetaDf$name, .etas))) {
      .nlmeResidStop("a residual error that depends on the model-calculated '", .v, "'")
    }
  }
  .eta <- intersect(all.vars(.expr), .etas)
  if (length(.eta) > 0L) {
    .nlmeResidStop("an eta on the residual error ('", .eta[1], "')")
  }
  .th <- intersect(all.vars(.expr), .thetaDf$name)
  .shared <- intersect(.th, .struct)
  if (length(.shared) > 0L) {
    .nlmeResidStop("a residual error sharing '", .shared[1], "' with the structural model")
  }
  .fix <- .thetaDf$name[.thetaDf$fix]
  .est <- setdiff(.th, .fix)
  if (length(.est) == 0L) {
    .nlmeResidStop("a residual SD without an estimated scale parameter")
  }
  if (length(.est) > 1L) {
    .nlmeResidStop(
      "a residual SD with more than one estimated parameter (",
      paste0("'", .est, "'", collapse = ", "),
      ")"
    )
  }
  .fixTh <- intersect(.th, .fix)
  .pw <- which(.iniDf$err == "pow2")
  list(
    slot = .slot,
    errType = .errType,
    expr = .expr,
    theta = .est,
    fixed = setNames(.thetaDf$est[match(.fixTh, .thetaDf$name)], .fixTh),
    covs = intersect(all.vars(.expr), .covs),
    exclude = .th,
    pow = if (length(.pw) == 1L) {
      list(name = .iniDf$name[.pw], ini = .iniDf$est[.pw], fix = .iniDf$fix[.pw])
    }
  )
}

#' @export
rxUiGet.nlmeResidInfo <- function(x, ...) {
  .nlmeResidInfo(x[[1]])
}
attr(rxUiGet.nlmeResidInfo, "rstudio") <- list()

#' Evaluate the residual SD expression on the observations
#'
#' @param info `.nlmeResidInfo()` result
#' @param data observation rows of the nlme data
#' @param theta named values of the estimated residual thetas
#' @return numeric vector, one per observation
#' @noRd
.nlmeResidEval <- function(info, data, theta) {
  .lst <- c(
    as.list(data[, info$covs, drop = FALSE]),
    list(time = data$TIME, t = data$TIME),
    as.list(info$fixed),
    as.list(theta)
  )
  .env <- list2env(.lst, parent = asNamespace("rxode2"))
  .ret <- try(eval(info$expr, .env), silent = TRUE)
  if (inherits(.ret, "try-error") || !is.numeric(.ret)) {
    .nlmeResidStop("the residual SD '", deparse1(info$expr), "' (it cannot be evaluated in R)")
  }
  rep_len(.ret, nrow(data))
}

#' nlme variance covariate for a time/covariate scaled residual SD
#'
#' Checks that the SD is the estimated theta times a data-only factor `g` and
#' returns `g^2`, the variance covariate of `nlme::varFixed()`.
#'
#' @param info `.nlmeResidInfo()` result
#' @param data observation rows of the nlme data
#' @return numeric vector, one per observation
#' @noRd
.nlmeResidVarCov <- function(info, data) {
  .g1 <- .nlmeResidEval(info, data, setNames(1, info$theta))
  .g2 <- .nlmeResidEval(info, data, setNames(2, info$theta))
  .ok <- all(is.finite(.g1)) &&
    all(is.finite(.g2)) &&
    all(.g1 != 0) &&
    isTRUE(all.equal(.g2, 2 * .g1, tolerance = 1e-8))
  if (!.ok) {
    .nlmeResidStop(
      "the residual SD '",
      deparse1(info$expr),
      "' (it must be '",
      info$theta,
      "' times a nonzero function of time and covariates)"
    )
  }
  .g1^2
}

#' nlme weights for a time/covariate scaled residual SD
#'
#' @param info `.nlmeResidInfo()` result
#' @return nlme `varFunc`
#' @noRd
.nlmeResidWeights <- function(info) {
  .fixed <- nlme::varFixed(~rxNlmeResidV)
  switch(
    info$errType,
    add = .fixed,
    prop = nlme::varComb(.fixed, nlme::varPower(fixed = list(power = 1))),
    pow = nlme::varComb(.fixed, .nlmeResidVarPower(info))
  )
}

#' nlme `varPower()` for the exponent of a `pow()` residual error
#'
#' @param info `.nlmeResidInfo()` result
#' @return nlme `varFunc`, holding a fixed exponent at its value
#' @noRd
.nlmeResidVarPower <- function(info) {
  if (isTRUE(info$pow$fix)) {
    return(nlme::varPower(fixed = list(power = info$pow$ini)))
  }
  nlme::varPower(value = info$pow$ini)
}

#' Thetas only the residual error uses, kept out of nlme's fixed effects
#'
#' @param ui rxode2 ui
#' @return character vector
#' @noRd
.nlmeResidExclude <- function(ui) {
  .info <- .nlmeResidInfo(ui)
  if (is.null(.info)) {
    return(character(0))
  }
  .info$exclude
}

#' Structural parameters nlme estimates
#'
#' @param x list holding the rxode2 ui
#' @param ... ignored
#' @return character vector
#' @noRd
.nlmeParamsToEstimate <- function(x, ...) {
  setdiff(rxUiGet.saemParamsToEstimate(x, ...), .nlmeResidExclude(x[[1]]))
}
