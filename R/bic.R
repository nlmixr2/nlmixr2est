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
    .dfSub <- max(0, .cls$omega - max(0, .cls$total - .df))
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
  .isTheta <- !is.na(.idf$ntheta)
  # a same() copy row is not separately estimated
  .cnd <- .idf$condition
  .copy <- !.isTheta & !is.na(.cnd) & .cnd != lotri::lotriBaseCondition(.cnd)
  .est <- !.idf$fix & !.copy
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
  .mu <- ui$muRefDataFrame$theta
  .str <- setdiff(intersect(.thetas, union(.link$structural, .mu)), .res)
  .ret$structural <- length(.str)
  .ret$residual <- length(.res)
  .ret
}

#' Flatten a chain of one associative operator family into its operands
#' @param x expression
#' @param ops operator names of the family
#' @return list of operand expressions
#' @noRd
.bicFlatten <- function(x, ops) {
  if (is.call(x) && as.character(x[[1]]) %in% ops) {
    return(do.call(c, lapply(as.list(x)[-1], .bicFlatten, ops = ops)))
  }
  if (is.call(x) && identical(as.character(x[[1]]), "(")) {
    return(.bicFlatten(x[[2]], ops))
  }
  list(x)
}

#' The parts of an expression that make up an eta's parameter
#'
#' The largest sub-expressions that hold an eta and no state-dependent
#' symbol, with the state-free factors of a product they sit in (other
#' than separately assigned variables), so
#' `ka * depot - exp(tcl + eta.cl) * center / v` gives `exp(tcl + eta.cl)`.
#' @param x expression
#' @param etaLike eta (and eta-only variable) names
#' @param dirty state and state-dependent names
#' @param own assigned variable names; a bare one is its own parameter
#' @return list of expressions
#' @noRd
.bicEtaPart <- function(x, etaLike, dirty, own) {
  .v <- all.vars(x)
  if (!any(.v %in% etaLike)) {
    return(list())
  }
  if (!any(.v %in% dirty) || !is.call(x)) {
    return(list(x))
  }
  .f <- as.character(x[[1]])
  .rec <- function(ops) {
    lapply(.bicFlatten(x, ops), function(o) {
      .ov <- all.vars(o)
      list(clean = !any(.ov %in% dirty), eta = any(.ov %in% etaLike), x = o)
    })
  }
  .sub <- function(lst) {
    do.call(c, lapply(lst, function(o) .bicEtaPart(o$x, etaLike, dirty, own)))
  }
  if (.f %in% c("*", "/")) {
    .o <- .rec(c("*", "/"))
    .clean <- Filter(function(o) o$clean && !(is.name(o$x) && as.character(o$x) %in% own), .o)
    return(c(lapply(.clean, function(o) o$x), .sub(Filter(function(o) !o$clean, .o))))
  }
  if (.f %in% c("+", "-", "(")) {
    return(.sub(.rec(c("+", "-"))))
  }
  do.call(c, lapply(as.list(x)[-1], .bicEtaPart, etaLike = etaLike, dirty = dirty, own = own))
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
  # each statement as list(lhs = name or NA, x = rhs, err = residual spec?)
  .lines <- list()
  .preds <- character(0)
  .add <- function(lhs, x, err) {
    .lines[[length(.lines) + 1L]] <<- list(lhs = lhs, x = x, rhs = all.vars(x), err = err)
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
      if (is.name(x[[2]])) {
        .preds <<- c(.preds, as.character(x[[2]]))
      }
      # a likelihood endpoint, ll(err) ~ expr, is structural, not an error spec
      .ll <- is.call(x[[2]]) && identical(as.character(x[[2]][[1]]), "ll")
      .add(NA_character_, x[[3]], !.ll)
    }
    invisible()
  }
  lapply(ui$lstExpr, .walk)
  .lhs <- vapply(.lines, function(l) l$lhs, character(1))
  .rhsOf <- function(v) {
    unique(unlist(lapply(.lines[which(.lhs == v)], function(l) l$rhs)))
  }
  .vars <- unique(.lhs[!is.na(.lhs)])
  # stopping points: predictions and anything depending on a state or
  # linCmt(), and "random terms" built only from etas and covariates
  .stateDep <- unique(c(.preds, .lhs[vapply(.lines, function(l) "linCmt" %in% all.names(l$x), logical(1))]))
  .stateDep <- .stateDep[!is.na(.stateDep)]
  .randTerm <- character(0)
  repeat {
    .n <- length(.stateDep) + length(.randTerm)
    for (.v in .vars) {
      .r <- .rhsOf(.v)
      if (any(.r %in% c(ui$state, .stateDep))) {
        .stateDep <- union(.stateDep, .v)
      }
      .eta <- intersect(.r, c(.etas, .randTerm))
      .other <- setdiff(.r, c(.etas, .randTerm))
      if (length(.eta) > 0L && !any(.other %in% c(.thetas, .vars))) {
        .randTerm <- union(.randTerm, .v)
      }
    }
    if (length(.stateDep) + length(.randTerm) == .n) {
      break
    }
  }
  .etaLike <- c(.etas, .randTerm)
  .dirty <- c(ui$state, .stateDep)
  .isP <- vapply(.lines, function(l) any(l$rhs %in% .etaLike), logical(1)) &
    !(.lhs %in% .randTerm)
  .pVars <- unique(.lhs[.isP & !is.na(.lhs)])
  # variables feeding a residual error spec: the spec's own symbols, then
  # through state-free assignments
  .errFeed <- unique(unlist(lapply(Filter(function(l) l$err, .lines), function(l) l$rhs)))
  .todo <- setdiff(intersect(.errFeed, .vars), .preds)
  while (length(.todo) > 0L) {
    .new <- setdiff(unlist(lapply(.todo, .rhsOf)), .errFeed)
    .errFeed <- c(.errFeed, .new)
    .todo <- setdiff(intersect(.new, .vars), .stateDep)
  }
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
    # an error spec's terms (add(), prop(), ...) are separate parameters
    .terms <- if (.l$err) .bicFlatten(.l$x, "+") else list(.l$x)
    .parts <- do.call(
      c,
      lapply(.terms, function(t) {
        if (.l$err && any(all.vars(t) %in% .etaLike)) list(t) else .bicEtaPart(t, .etaLike, .dirty, .vars)
      })
    )
    .t <- .trace(unique(unlist(lapply(.parts, all.vars))))
    if (.l$err || (!is.na(.l$lhs) && .l$lhs %in% .errFeed)) {
      .res <- union(.res, .t)
    } else {
      .str <- union(.str, .t)
    }
  }
  list(structural = .str, residual = .res)
}
