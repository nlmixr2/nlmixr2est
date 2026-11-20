#' Bayesian information criterion for nlmixr2 fits
#'
#' `BIC(fit)` charges every estimated parameter `log(nobs)`.  The `type`
#' variants
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
#' The linkage is read from the expression around each eta: the whole sum it
#' is added to (`tcl + eta.cl`), or the product it multiplies, where the other
#' factors are its typical value (`tvcl * exp(eta.cl)`) unless the eta's own
#' factor already has one; then only covariate effects such as
#' `(WT / 70)^cwt` join it (in `exp(tcl + eta.cl) / v`, `v` stays its own
#' parameter).  This is exact for mu-referenced parameters; for other forms
#' it is an approximation, so write each eta on its own parameter line for an
#' exact count.  Fixed parameters are not counted.  A theta whose linkage
#' cannot be determined is charged `log(nobs)`, as in `BIC(fit)`.
#'
#' @param object nlmixr2 fit
#' @param ... other fitted models to compare, giving a data frame with
#'   columns `df` and `BIC`, one row per model; a character value here is
#'   taken as `type`, so `BIC(fit, "mixed")` works
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
#' @return the BIC value, or a data frame when several models are given
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
#' BIC(fit)
#' BIC(fit, type = "mixed")
#' }
BIC.nlmixr2FitCore <- function(object, ..., type = c("fixed", "mixed", "iiv", "random")) {
  .dots <- list(...)
  .chr <- vapply(.dots, is.character, logical(1))
  if (any(.chr) && missing(type)) {
    type <- .dots[[which(.chr)[1]]]
  }
  type <- match.arg(type, c("fixed", "mixed", "iiv", "random"))
  .objs <- c(list(object), .dots[!.chr])
  if (length(.objs) == 1L) {
    return(as.numeric(.bicValue(object, type)))
  }
  .val <- lapply(.objs, function(.o) {
    if (inherits(.o, c("nlmixr2FitCore", "nlmixr2FitCoreSilent"))) {
      .b <- .bicValue(.o, type)
      c(attr(.b, "dfObs") + attr(.b, "dfSub"), as.numeric(.b))
    } else {
      c(attr(stats::logLik(.o), "df"), stats::BIC(.o))
    }
  })
  .ret <- data.frame(
    df = vapply(.val, `[`, numeric(1), 1L),
    BIC = vapply(.val, `[`, numeric(1), 2L)
  )
  .call <- match.call(expand.dots = TRUE)
  .call$type <- NULL
  .args <- as.list(.call[-1L])
  .args <- .args[!vapply(.args, is.character, logical(1))]
  row.names(.ret) <- make.unique(vapply(.args, deparse1, character(1)))
  .ret
}

#' @rdname BIC.nlmixr2FitCore
#' @export
BIC.nlmixr2FitCoreSilent <- BIC.nlmixr2FitCore

#' BIC with its penalty counts
#'
#' @param fit nlmixr2 fit
#' @param type BIC variant, see `BIC.nlmixr2FitCore()`
#' @return the BIC value with the attributes `nobs`, `nsub`, `dfObs` and
#'   `dfSub` (parameters charged `log(nobs)` and `log(nsub)`)
#' @noRd
.bicValue <- function(fit, type) {
  assertNlmixrFit(fit)
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
  .mu <- setdiff(ui$muRefDataFrame$theta, .link$residual)
  .str <- setdiff(intersect(.thetas, union(.link$structural, .mu)), .err)
  .res <- setdiff(intersect(.thetas, union(.link$residual, .link$structural)), .str)
  .ret$structural <- length(.str)
  .ret$residual <- length(.res)
  .ret
}

#' Flatten a chain of one associative operator family into its operands
#'
#' Parentheses and unary signs are transparent.
#' @param x expression
#' @param ops operator names of the family
#' @return list of operand expressions
#' @noRd
.bicFlatten <- function(x, ops) {
  if (is.call(x)) {
    .f <- as.character(x[[1]])
    if (.f == "(" || (.f %in% c("+", "-") && length(x) == 2L)) {
      return(.bicFlatten(x[[2]], ops))
    }
    if (.f %in% ops) {
      return(do.call(c, lapply(as.list(x)[-1], .bicFlatten, ops = ops)))
    }
  }
  list(x)
}

#' The parts of an expression that make up an eta's parameter
#'
#' A sum with a bare direct eta term (`tcl + eta.cl`) is the parameter;
#' otherwise only the terms holding an eta are.  In a product the other
#' factors are the eta's typical value (`tvcl * exp(eta.cl)`) unless the
#' eta's own factor already carries one (`exp(tcl + eta.cl) / v`); then only
#' covariate effects join.
#' @param x expression
#' @param etaLike eta and eta-carrying variable names
#' @param direct etas and variables built only from etas and covariates
#' @param typ thetas and eta-carrying variables that hold a theta
#' @param own assigned variables that carry no eta
#' @param covs covariate names
#' @return list of expressions
#' @noRd
.bicEtaPart <- function(x, etaLike, direct, typ, own, covs) {
  .has <- function(e) any(all.vars(e) %in% etaLike)
  if (!.has(x)) {
    return(list())
  }
  if (!is.call(x)) {
    return(list(x))
  }
  .rec <- function(lst) {
    do.call(c, lapply(lst, .bicEtaPart, etaLike = etaLike, direct = direct, typ = typ, own = own, covs = covs))
  }
  .f <- as.character(x[[1]])
  if (.f %in% c("*", "/")) {
    .o <- .bicFlatten(x, c("*", "/"))
    .eta <- vapply(.o, .has, logical(1))
    .typ <- any(vapply(.o[.eta], function(o) any(all.vars(o) %in% typ), logical(1)))
    .keep <- !.eta
    if (.typ) {
      # only a covariate effect joins an eta that has its typical value
      .keep <- .keep &
        vapply(
          .o,
          function(o) {
            .v <- all.vars(o)
            !any(.v %in% own) && (!any(.v %in% typ) || any(.v %in% covs))
          },
          logical(1)
        )
    }
    return(c(.o[.keep], .rec(.o[.eta])))
  }
  if (.f %in% c("+", "-", "(")) {
    .o <- .bicFlatten(x, c("+", "-"))
    .bare <- vapply(.o, function(o) is.name(o) && as.character(o) %in% direct, logical(1))
    if (any(.bare)) {
      return(list(x))
    }
    return(.rec(.o[vapply(.o, .has, logical(1))]))
  }
  .rec(as.list(x)[-1])
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
  .acc <- new.env(parent = emptyenv())
  .acc$lines <- list()
  .acc$preds <- character(0)
  .add <- function(lhs, x, err) {
    .acc$lines[[length(.acc$lines) + 1L]] <- list(lhs = lhs, x = x, rhs = all.vars(x), err = err)
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
        .acc$preds <- c(.acc$preds, as.character(x[[2]]))
      }
      # a likelihood endpoint, ll(err) ~ expr, is structural, not an error spec
      .ll <- is.call(x[[2]]) && identical(as.character(x[[2]][[1]]), "ll")
      .add(NA_character_, x[[3]], !.ll)
    }
    invisible()
  }
  lapply(ui$lstExpr, .walk)
  .lines <- .acc$lines
  .preds <- .acc$preds
  .lhs <- vapply(.lines, function(l) l$lhs, character(1))
  .rhsOf <- function(v) {
    unique(unlist(lapply(.lines[which(.lhs == v)], function(l) l$rhs)))
  }
  .vars <- unique(.lhs[!is.na(.lhs)])
  .grow <- function(set, test) {
    repeat {
      .n <- length(set)
      for (.v in setdiff(.vars, set)) {
        if (test(.rhsOf(.v), set)) {
          set <- c(set, .v)
        }
      }
      if (length(set) == .n) {
        return(set)
      }
    }
  }
  # predictions and what depends on a state or linCmt() stop the tracing
  .linCmt <- .lhs[vapply(.lines, function(l) "linCmt" %in% all.names(l$x), logical(1))]
  .stateDep <- .grow(
    unique(c(.preds, .linCmt[!is.na(.linCmt)])),
    function(r, set) any(r %in% c(ui$state, set))
  )
  .thetaFree <- .grow(character(0), function(r, set) {
    !any(r %in% .thetas) && all(intersect(r, .vars) %in% set)
  })
  .etaVar <- .grow(character(0), function(r, set) any(r %in% c(.etas, set)))
  .etaVar <- setdiff(.etaVar, .stateDep)
  .direct <- c(.etas, intersect(.etaVar, .thetaFree))
  .closure <- function(seed, through) {
    .out <- seed
    .todo <- intersect(seed, through)
    while (length(.todo) > 0L) {
      .new <- setdiff(unlist(lapply(.todo, .rhsOf)), .out)
      .out <- c(.out, .new)
      .todo <- intersect(.new, through)
    }
    .out
  }
  .predFeed <- .closure(.preds, .vars)
  .errSeed <- unique(unlist(lapply(Filter(function(l) l$err, .lines), function(l) l$rhs)))
  # an error spec's own variables, then state-free assignments feeding them
  .first <- unlist(lapply(setdiff(intersect(.errSeed, .vars), .preds), .rhsOf))
  .errFeed <- .closure(unique(c(.errSeed, .first)), setdiff(.vars, c(.preds, .stateDep)))
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
      } else if (.s %in% .vars && !(.s %in% c(.stateDep, .etaVar))) {
        .todo <- c(.todo, .rhsOf(.s))
      }
    }
    unique(.out)
  }
  .str <- character(0)
  .res <- character(0)
  for (.l in .lines) {
    if (!is.na(.l$lhs) && .l$lhs %in% .direct) {
      next
    }
    # a state-dependent statement takes only a direct eta, and a residual one
    # no eta that reaches it through the prediction
    .dirty <- any(.l$rhs %in% c(ui$state, .stateDep))
    .etaLike <- if (.dirty) .direct else c(.etas, .etaVar)
    .resLine <- .l$err || (!is.na(.l$lhs) && .l$lhs %in% setdiff(.errFeed, .predFeed))
    if (.resLine) {
      .etaLike <- setdiff(.etaLike, .predFeed)
    }
    if (!any(.l$rhs %in% .etaLike)) {
      next
    }
    .parts <- .bicEtaPart(
      .l$x,
      .etaLike,
      intersect(.direct, .etaLike),
      c(.thetas, setdiff(.etaVar, .thetaFree)),
      setdiff(.vars, .etaLike),
      ui$allCovs
    )
    .t <- .trace(unique(unlist(lapply(.parts, all.vars))))
    if (.resLine) {
      .res <- union(.res, .t)
    } else {
      .str <- union(.str, .t)
    }
  }
  list(structural = .str, residual = .res)
}
