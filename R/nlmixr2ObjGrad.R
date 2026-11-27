#' Objective function and analytic gradient of a FOCEi-family fit at new parameters
#'
#' Evaluates the fit's own objective function and its analytic outer gradient
#' (the `foceiControl(fast=TRUE)` C++ gradient) at user-supplied population
#' parameters, without running an outer optimization.  The fit's estimation
#' method (`"focei"`, `"foce"`, `"laplace"`, `"agq"`, ...) is reused, since it
#' determines the objective being differentiated.
#'
#' At the fit's own estimates the fit's empirical Bayes estimates are used as-is,
#' so the objective equals `fit$objective`.  At any other point the etas are
#' re-optimized (warm-started from the fit's) before the objective and gradient
#' are evaluated.
#'
#' @param fit A fit from a FOCEi-family estimation method.
#' @param theta Population fixed effects and residual parameters to evaluate
#'   at, on the scale they are written in `ini()`.  Either named (a subset of
#'   `names(fit$theta)`; the rest stay at the fit's values) or unnamed and as
#'   long as `fit$theta`.  `NULL` uses `fit$theta`.
#' @param omega Between-subject covariance matrix to evaluate at, with the same
#'   dimensions as `fit$omega`.  `NULL` uses `fit$omega`.
#' @param scale The scale of the returned gradient for the omega elements.
#'   `"untransformed"` differentiates with respect to the variances and
#'   covariances themselves, named as in `ini()` (`eta.cl`, `(eta.cl,eta.v)`).
#'   `"transformed"` differentiates with respect to the Cholesky-based
#'   estimation-scale omega parameters the optimizer uses (`om.chol.<k>`, with
#'   `foceiControl(diagXform=)` applied to the diagonal).  Thetas are always on
#'   their `ini()` scale.
#' @param individual Also return each subject's contribution to the objective
#'   and gradient.
#' @return A list with `objective` (the fit's objective at the point),
#'   `gradient` (the derivative of `objective` for every estimated parameter,
#'   named), `par` (the parameter values on the same scale as `gradient`) and
#'   `scale`.  With `individual = TRUE` it also has `individual`, a data frame of
#'   each subject's `ID`, objective contribution (`OBJI`) and gradient (one
#'   column per parameter), and `population`, the part of the objective and
#'   gradient that belongs to no subject (such as an `ini()` prior), so the
#'   subjects plus `population` sum to the totals.  The attribute `nAnalyticGradDirect` counts the analytic gradient
#'   evaluations behind the result.
#' @export
#' @examples
#' \donttest{
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45
#'     tcl <- 1
#'     tv <- 3.45
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     add.sd <- 0.7
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv)
#'     d / dt(depot) <- -ka * depot
#'     d / dt(center) <- ka * depot - cl / v * center
#'     cp <- center / v
#'     cp ~ add(add.sd)
#'   })
#' }
#' fit <- nlmixr2(one.cmt, theo_sd, "focei", control = foceiControl(print = 0))
#' nlmixr2ObjGrad(fit)
#' nlmixr2ObjGrad(fit, theta = c(tka = 0.5))
#' }
nlmixr2ObjGrad <- function(
  fit,
  theta = NULL,
  omega = NULL,
  scale = c("untransformed", "transformed"),
  individual = FALSE
) {
  checkmate::assertLogical(individual, len = 1, any.missing = FALSE)
  scale <- match.arg(scale)
  if (!inherits(fit, "nlmixr2FitCore")) {
    stop("'fit' must be an nlmixr2 fit", call. = FALSE)
  }
  .env <- fit$env
  if (is.null(.env$foceiControl) || is.null(.env$est)) {
    stop("'fit' must come from a FOCEi-family estimation method", call. = FALSE)
  }
  .th0 <- fit$theta
  .om0 <- fit$omega
  .th <- .objGradTheta(theta, .th0)
  .om <- .objGradOmega(omega, .om0)
  .atFit <- isTRUE(all.equal(unname(.th), unname(.th0), tolerance = 0)) &&
    (is.null(.om0) || isTRUE(all.equal(unname(.om), unname(.om0), tolerance = 0)))
  .ui <- .fitUiCopy(fit)
  if (is.null(.ui)) {
    stop("could not copy the fit's model", call. = FALSE)
  }
  .ui <- .objGradPinUi(.ui, .th, .om)
  .r <- .foceiObjGradRefit(fit, .ui, forceFast = TRUE, reoptimize = !.atFit)
  if (is.null(.r) || is.null(.r$gradient)) {
    stop("no analytic gradient for this fit; see foceiControl(fast=)", call. = FALSE)
  }
  if (is.null(names(.r$gradient))) {
    stop("could not name the analytic gradient for this fit", call. = FALSE)
  }
  .tr <- .objGradOmegaScale(.r$gradient, .ui, .om, scale, .r$gradInd)
  .par <- c(.th[intersect(names(.tr$gradient), names(.th))], .tr$par)
  .ret <- list(
    objective = .r$objective,
    gradient = .tr$gradient,
    par = .par[names(.tr$gradient)],
    scale = scale
  )
  if (individual) {
    .ret <- c(.ret, .objGradIndividual(.r, .tr))
  }
  structure(
    .ret,
    nAnalyticGradDirect = .r$nAnalyticGradDirect
  )
}

#' Fill the supplied thetas into the fit's theta vector
#' @noRd
.objGradTheta <- function(theta, th0) {
  if (is.null(theta)) {
    return(th0)
  }
  if (!is.numeric(theta) || anyNA(theta)) {
    stop("'theta' must be a numeric vector without missing values", call. = FALSE)
  }
  if (is.null(names(theta))) {
    if (length(theta) != length(th0)) {
      stop("an unnamed 'theta' must have length ", length(th0), call. = FALSE)
    }
    return(stats::setNames(as.numeric(theta), names(th0)))
  }
  .bad <- setdiff(names(theta), names(th0))
  if (length(.bad) > 0L) {
    stop("unknown 'theta' names: ", paste(.bad, collapse = ", "), call. = FALSE)
  }
  th0[names(theta)] <- as.numeric(theta)
  th0
}

#' Validate the supplied omega against the fit's
#' @noRd
.objGradOmega <- function(omega, om0) {
  if (is.null(omega)) {
    return(om0)
  }
  omega <- as.matrix(omega)
  if (is.null(om0) || !identical(dim(omega), dim(om0))) {
    stop("'omega' must have the same dimensions as 'fit$omega'", call. = FALSE)
  }
  if (!isSymmetric(unname(omega)) || anyNA(omega)) {
    stop("'omega' must be a symmetric matrix", call. = FALSE)
  }
  dimnames(omega) <- dimnames(om0)
  omega
}

#' Pin theta and omega into a ui's `iniDf`
#' @noRd
.objGradPinUi <- function(ui, theta, omega) {
  .ini <- ui$iniDf
  .w <- match(names(theta), .ini$name)
  .ini$est[.w[!is.na(.w)]] <- as.numeric(theta)[!is.na(.w)]
  if (!is.null(omega)) {
    .eta <- which(!is.na(.ini$neta1))
    for (.i in .eta) {
      .r <- .ini$neta1[.i]
      .c <- .ini$neta2[.i]
      if (.r <= nrow(omega) && .c <= ncol(omega)) {
        .ini$est[.i] <- omega[.r, .c]
      }
    }
  }
  ui$iniDf <- .ini
  ui
}

#' Re-enter estimation at a pinned ui and read back the objective and gradient
#'
#' The gradient is the one `analyticOuterGrad()` stashes in the fit's own C++ path
#' (natural scale: thetas, sigmas, estimation-scale omega), named by
#' `.foceiGradNames()`.  `NULL` gradient when the analytic route declined.
#' @param fit completed fit (source of data, etas and control)
#' @param ui ui with the evaluation point pinned in `iniDf`
#' @param forceFast turn on `fast=TRUE` for the evaluation
#' @param reoptimize re-optimize the etas (warm-started at the fit's)
#' @return list(objective, gradient), or `NULL` if the re-entry failed
#' @noRd
.foceiObjGradRefit <- function(fit, ui, forceFast = FALSE, reoptimize = FALSE) {
  tryCatch(
    {
      .env <- if (rxode2::rxIs(fit, "nlmixr2FitData")) fit$env else fit
      .est <- .env$est
      if (is.null(.est) || !nzchar(.est)) {
        return(NULL)
      }
      .control <- .env$foceiControl
      if (!reoptimize) {
        .control$maxInnerIterations <- 0L # evaluate at the fit's EBEs
      }
      .control$maxOuterIterations <- 0L # no outer step: the gradient is at THIS theta
      .control$calcTables <- FALSE
      .control$covMethod <- 0L
      .control$skipCov <- fit$skipCov
      if (forceFast) {
        .control$fast <- TRUE
        # no outer step is taken, but a derivative-free optimizer would revert fast
        if (.control$outerOptTxt %in% c("bobyqa", "uobyqa", "newuoa")) {
          .control$outerOpt <- -1L
          .control$outerOptFun <- .nlminb
          .control$outerOptTxt <- "nlminb"
        }
      }
      # Re-run under the fit's OWN est, not est="none": the gradient SHAPE is the
      # estimation method (FOCE freezes the residual variance, AGQ adds quadrature).
      .etaMat <- .fitEtaMat(fit)
      if (!is.null(.etaMat)) {
        .control$etaMat <- .etaMat
      }
      # the nested re-fit resets mu-referencing global state (.muRefTrans$cur); restore it
      .savedMuRef <- .muRefTrans$cur
      on.exit(.muRefTrans$cur <- .savedMuRef, add = TRUE)
      .f2 <- suppressMessages(suppressWarnings(
        nlmixr2(ui, data = getData(fit), est = .est, control = .control)
      ))
      .src <- tryCatch(.f2$env, error = function(.) NULL)
      if (is.null(.src)) {
        return(NULL)
      }
      .g <- NULL
      if (exists(".gradDirectFirst", .src, inherits = FALSE)) {
        .uiG <- tryCatch(.f2$ui, error = function(e) NULL)
        if (is.null(.uiG)) {
          .uiG <- ui
        }
        .g <- .foceiGradNames(as.numeric(get(".gradDirectFirst", .src)), .uiG, fit, .src)
      }
      .nDirect <- if (exists("nAnalyticGradDirect", .src, inherits = FALSE)) {
        as.integer(get("nAnalyticGradDirect", .src))
      } else {
        0L
      }
      .gInd <- NULL
      if (!is.null(.g) && exists(".gradDirectFirstInd", .src, inherits = FALSE)) {
        .gInd <- t(get(".gradDirectFirstInd", .src))
        if (ncol(.gInd) == length(.g)) {
          colnames(.gInd) <- names(.g)
        } else {
          .gInd <- NULL
        }
      }
      .etaObf <- tryCatch(.f2$etaObf, error = function(e) NULL)
      list(
        objective = .f2$objective,
        gradient = .g,
        nAnalyticGradDirect = .nDirect,
        gradInd = .gInd,
        etaObf = .etaObf
      )
    },
    error = function(e) NULL
  )
}

#' Map the estimation-scale omega part of a named gradient to the requested scale
#' @return list(gradient, par) where `par` holds the omega values on that scale
#' @noRd
.objGradOmegaScale <- function(g, ui, omega, scale, gInd = NULL) {
  .isOm <- grepl("^om[.]chol[.][0-9]+$", names(g))
  if (!any(.isOm)) {
    return(list(gradient = g, par = numeric(0), individual = gInd))
  }
  .k <- as.integer(sub("^om[.]chol[.]", "", names(g)[.isOm]))
  .oe <- .foceiEstOmegaDeriv(ui, omega)
  .rxInv <- rxode2::rxSymInvCholCreate(
    mat = omega,
    diag.xform = rxode2::rxGetControl(ui, "diagXform", "sqrt")
  )
  if (is.null(.oe) || max(.k) > length(.oe$dOi)) {
    stop("could not map the omega gradient to the estimation scale", call. = FALSE)
  }
  if (scale == "transformed") {
    return(list(
      gradient = g,
      par = stats::setNames(.rxInv$theta[.k], names(g)[.isOm]),
      individual = gInd
    ))
  }
  # the estimated (not fixed) omega elements, as written in ini()
  .ini <- ui$iniDf
  .ini <- .ini[!is.na(.ini$neta1) & !.ini$fix, , drop = FALSE]
  if (nrow(.ini) != length(.k)) {
    stop("cannot map this omega structure to the untransformed scale; use scale=\"transformed\"", call. = FALSE)
  }
  # d(Omega)/d(chol_k) = -Omega d(Omega^-1)/d(chol_k) Omega
  .jac <- vapply(
    .k,
    function(k) {
      .d <- -omega %*% .oe$dOi[[k]] %*% omega
      .d[cbind(.ini$neta1, .ini$neta2)]
    },
    numeric(nrow(.ini))
  )
  .jac <- matrix(.jac, nrow = nrow(.ini))
  .gn <- tryCatch(solve(t(.jac), as.numeric(g[.isOm])), error = function(e) NULL)
  if (is.null(.gn)) {
    stop("singular omega Jacobian; use scale=\"transformed\"", call. = FALSE)
  }
  .gn <- stats::setNames(.gn, .ini$name)
  if (!is.null(gInd)) {
    .gi <- t(solve(t(.jac), t(gInd[, .isOm, drop = FALSE])))
    colnames(.gi) <- .ini$name
    gInd <- cbind(gInd[, !.isOm, drop = FALSE], .gi)
  }
  list(
    gradient = c(g[!.isOm], .gn),
    par = stats::setNames(omega[cbind(.ini$neta1, .ini$neta2)], .ini$name),
    individual = gInd
  )
}

#' Per-subject objective and gradient, and the population remainder
#' @param r result of `.foceiObjGradRefit()`
#' @param tr result of `.objGradOmegaScale()`
#' @return list(individual, population)
#' @noRd
.objGradIndividual <- function(r, tr) {
  .gi <- tr$individual
  .eo <- r$etaObf
  if (is.null(.gi) || is.null(.eo) || nrow(.eo) != nrow(.gi) || !("OBJI" %in% names(.eo))) {
    stop("per-subject contributions are not available for this fit", call. = FALSE)
  }
  .gi <- .gi[, names(tr$gradient), drop = FALSE]
  .ind <- data.frame(ID = .eo$ID, OBJI = .eo$OBJI, check.names = FALSE)
  .ind <- cbind(.ind, as.data.frame(.gi, optional = TRUE))
  names(.ind) <- c("ID", "OBJI", colnames(.gi))
  list(
    individual = .ind,
    population = list(
      objective = r$objective - sum(.eo$OBJI),
      gradient = tr$gradient - colSums(.gi)
    )
  )
}
