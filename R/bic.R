#' Bayesian information criterion variants for mixed-effects fits
#'
#' `BIC(fit)` charges every estimated parameter `log(nobs)`.  These variants
#' charge the parameters tied to random effects `log(nsub)` instead (Delattre
#' et al. 2014), which keeps BIC from favouring over-parameterised
#' random-effect structures.
#'
#' Parameters are classified from the fit's model:
#'
#' - Omega elements (including inter-occasion kappas) are always random-effect
#'   parameters.
#'
#' - A theta is a structural random-effect parameter when it enters the
#'   expression that defines a parameter carrying an eta (directly or through
#'   intermediate assignments, but not through ODE states), as with
#'   `tcl` and the covariate effect `cwt` in
#'   `cl <- exp(tcl + cwt * log(WT / 70) + eta.cl)`.
#'
#' - A residual-error theta is a random-effect parameter only when an eta sits
#'   on the residual error, as with `prop.sd` in
#'   `rsd <- prop.sd * exp(eta.rsd); cp ~ prop(rsd)`.
#'
#' Fixed parameters are not counted.  A theta whose linkage cannot be
#' determined is charged `log(nobs)`, as in `BIC(fit)`.
#'
#' @param fit nlmixr2 fit
#' @param type which parameters are charged `log(nsub)`:
#'
#' - `"fixed"`: none; equal to `BIC(fit)`.
#'
#' - `"mixed"`: omega elements and the thetas of eta-bearing structural
#'   parameters.
#'
#' - `"iiv"`: omega (and kappa) elements only.
#'
#' - `"random"`: every parameter tied to a random effect, which adds to
#'   `"mixed"` the residual-error thetas when an eta sits on the residual
#'   error.
#'
#' @return the BIC value.  It carries the attributes `nobs` and `nsub` (the
#'   counts the penalties use) and `dfObs` and `dfSub` (the number of
#'   parameters charged `log(nobs)` and `log(nsub)`).
#' @references Delattre M, Lavielle M, Poursat MA (2014). A note on BIC in
#'   mixed-effects models. Electronic Journal of Statistics, 8(1), 456-475.
#' @author Matthew L. Fidler
#' @export
#' @examples
#' \donttest{
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45
#'     tcl <- log(2.7)
#'     tv <- 3.45
#'     eta.cl ~ 0.3
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#' fit <- nlmixr2(one.cmt, theo_sd, est = "focei",
#'                control = foceiControl(print = 0))
#' bic(fit, "fixed")
#' bic(fit, "mixed")
#' }
bic <- function(fit, type = c("fixed", "mixed", "iiv", "random")) {
  assertNlmixrFit(fit)
  type <- match.arg(type)
  .ll <- fit$logLik
  .df <- attr(.ll, "df")
  .nobs <- attr(.ll, "nobs")
  if (is.null(.nobs)) {
    .nobs <- fit$nobs
  }
  .nsub <- fit$env$nsub
  if (is.null(.nsub)) {
    .nsub <- length(unique(as.character(fit$ID)))
  }
  .dfSub <- 0
  if (type != "fixed") {
    .cls <- .bicParClass(fit$ui)
    .dfSub <- .cls$omega - max(0, .cls$total - .df)
    if (type %in% c("mixed", "random")) {
      .dfSub <- .dfSub + .cls$structural
    }
    if (type == "random") {
      .dfSub <- .dfSub + .cls$residual
    }
    .dfSub <- min(max(.dfSub, 0), .df)
  }
  .dfObs <- .df - .dfSub
  structure(
    -2 * as.numeric(.ll) + log(.nobs) * .dfObs + log(.nsub) * .dfSub,
    nobs = .nobs,
    nsub = .nsub,
    dfObs = .dfObs,
    dfSub = .dfSub
  )
}

#' Count a ui's estimated parameters by BIC penalty class
#'
#' @param ui rxode2 ui of the fit
#' @return list with `total` (estimated parameters in `iniDf`), `omega`
#'   (estimated omega/kappa elements), `structural` (thetas of eta-bearing
#'   structural parameters) and `residual` (residual-error thetas with an eta
#'   on the residual error)
#' @noRd
.bicParClass <- function(ui) {
  .idf <- ui$iniDf
  .est <- !.idf$fix
  .isTheta <- !is.na(.idf$ntheta)
  .ret <- list(
    total = sum(.est),
    omega = sum(.est & !.isTheta),
    structural = 0,
    residual = 0
  )
  .thetas <- .idf$name[.isTheta & .est]
  if (length(.thetas) == 0L) {
    return(.ret)
  }
  .link <- tryCatch(
    .bicThetaLink(ui),
    error = function(e) {
      warning("BIC could not link thetas to etas; they use log(nobs)", call. = FALSE)
      list(structural = character(0), residual = character(0))
    }
  )
  .err <- .idf$name[.isTheta & !is.na(.idf$err)]
  .res <- intersect(.thetas, union(.link$residual, intersect(.link$structural, .err)))
  .str <- setdiff(intersect(.thetas, union(.link$structural, ui$muRefTable$theta)), .res)
  .ret$structural <- length(.str)
  .ret$residual <- length(.res)
  .ret
}

#' Thetas tied to an eta through the model's assignments
#'
#' @param ui rxode2 ui
#' @return list of theta names: `structural` (linked through a structural
#'   parameter) and `residual` (linked through the residual error)
#' @noRd
.bicThetaLink <- function(ui) {
  .idf <- ui$iniDf
  .thetas <- .idf$name[!is.na(.idf$ntheta)]
  .etas <- .idf$name[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2]
  .states <- ui$state
  # each assignment as list(lhs = name or NA, rhs = symbols, err = endpoint?)
  .lines <- list()
  .add <- function(lhs, rhs, err) {
    .lines[[length(.lines) + 1L]] <<- list(lhs = lhs, rhs = all.vars(rhs), err = err)
  }
  .walk <- function(x) {
    if (!is.call(x)) {
      return(invisible())
    }
    .f <- as.character(x[[1]])
    if (identical(.f, "{")) {
      lapply(as.list(x)[-1], .walk)
    } else if (identical(.f, "if")) {
      lapply(as.list(x)[-(1:2)], .walk)
    } else if (.f %in% c("<-", "=")) {
      .lhs <- if (is.name(x[[2]])) as.character(x[[2]]) else NA_character_
      .add(.lhs, x[[3]], FALSE)
    } else if (identical(.f, "~") && length(x) == 3L) {
      .add(NA_character_, x[[3]], TRUE)
    }
    invisible()
  }
  lapply(ui$lstExpr, .walk)
  .lhs <- vapply(.lines, function(l) l$lhs, character(1))
  .rhsOf <- function(v) {
    unique(unlist(lapply(.lines[which(.lhs == v)], function(l) l$rhs)))
  }
  .vars <- unique(.lhs[!is.na(.lhs)])
  # fixed points: variables that depend on an ODE state, and "random terms"
  # (variables built only from etas and covariates, e.g. ecl <- eta.cl)
  .stateDep <- character(0)
  .randTerm <- character(0)
  repeat {
    .n <- length(.stateDep) + length(.randTerm)
    for (.v in .vars) {
      .r <- .rhsOf(.v)
      if (any(.r %in% c(.states, .stateDep))) {
        .stateDep <- union(.stateDep, .v)
      }
      .eta <- intersect(.r, c(.etas, .randTerm))
      .other <- setdiff(.r, c(.etas, .randTerm))
      if (length(.eta) > 0L && !any(.other %in% c(.thetas, .vars))) {
        .randTerm <- union(.randTerm, .v)
      }
    }
    if (length(.stateDep) + length(.randTerm) == .n) break
  }
  .etaLike <- c(.etas, .randTerm)
  .isP <- vapply(.lines, function(l) any(l$rhs %in% .etaLike), logical(1)) &
    !(.lhs %in% .randTerm)
  .pVars <- unique(.lhs[.isP & !is.na(.lhs)])
  .errSyms <- unique(unlist(lapply(.lines[vapply(.lines, function(l) l$err, logical(1))], function(l) l$rhs)))
  .trace <- function(syms) {
    .seen <- character(0)
    .out <- character(0)
    .todo <- syms
    while (length(.todo) > 0L) {
      .s <- .todo[1]
      .todo <- .todo[-1]
      if (.s %in% .seen) {
        next
      }
      .seen <- c(.seen, .s)
      if (.s %in% .thetas) {
        .out <- c(.out, .s)
      } else if (.s %in% .vars && !(.s %in% c(.stateDep, .pVars, .randTerm))) {
        .todo <- c(.todo, .rhsOf(.s))
      }
    }
    unique(.out)
  }
  .str <- character(0)
  .res <- character(0)
  for (.i in which(.isP)) {
    .l <- .lines[[.i]]
    .t <- .trace(.l$rhs)
    if (.l$err || (!is.na(.l$lhs) && .l$lhs %in% .errSyms)) {
      .res <- union(.res, .t)
    } else {
      .str <- union(.str, .t)
    }
  }
  list(structural = .str, residual = .res)
}
