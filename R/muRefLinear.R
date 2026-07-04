#' Quick linear-model regression step for the mu-referenced FOCEI family
#'
#' Shared implementation for `muModel="lin"` (OLS) and `muModel="irls"`
#' (reweighted): regresses per-subject `phi` on covariate(s) to get a new
#' theta (intercept), covariate coefficients (slopes), and per-subject eta
#' residuals. Fixed covariate coefficients are subtracted from `phi` as a
#' known offset and not re-estimated.
#'
#' @param phi numeric vector, one back-calculated individual value per subject
#' @param cov data.frame/matrix of covariate columns, one row per subject;
#'   column names are the raw covariate names (e.g. `"logWT"`), not the
#'   theta names -- callers map the returned `coef` names back to
#'   covariate-coefficient theta names themselves
#' @param fixedCoef named numeric vector of user-fixed covariate
#'   coefficients (not re-estimated); names must be a subset of
#'   `colnames(cov)`
#' @param weights per-subject regression weight; `NULL` (the `"lin"` case)
#'   is equivalent to all-equal weights
#' @return `list(theta=, coef=, eta=, fit=)`
#' @author Matthew L. Fidler
#' @noRd
.muRefFit <- function(phi, cov, fixedCoef = NULL, weights = NULL) {
  covNames <- colnames(cov)
  fixedNames <- intersect(names(fixedCoef), covNames)
  freeNames <- setdiff(covNames, fixedNames)

  offset <- rep(0, length(phi))
  if (length(fixedNames) > 0) {
    fixedMat <- as.matrix(cov[, fixedNames, drop = FALSE])
    offset <- as.vector(fixedMat %*% fixedCoef[fixedNames])
  }
  y <- phi - offset

  if (is.null(weights)) weights <- rep(1, length(phi))

  if (length(freeNames) == 0L) {
    # only a population intercept left to estimate
    theta <- sum(weights * y) / sum(weights)
    eta <- y - theta
    coefAll <- fixedCoef[fixedNames]
    return(list(theta = theta, coef = coefAll, eta = eta, fit = NULL))
  }

  df <- as.data.frame(cov[, freeNames, drop = FALSE])
  df$y <- y
  fitForm <- stats::as.formula(paste0("y ~ ", paste(freeNames, collapse = " + ")))
  fit <- stats::lm(fitForm, data = df, weights = weights)

  theta <- unname(stats::coef(fit)[1])
  betaFree <- stats::coef(fit)[-1]
  names(betaFree) <- freeNames
  coefAll <- c(betaFree, fixedCoef[fixedNames])
  coefAll <- coefAll[covNames]
  eta <- stats::residuals(fit)

  list(theta = theta, coef = coefAll, eta = unname(eta), fit = fit)
}

#' Plain OLS mu-ref regression step (`muModel="lin"`)
#'
#' @inheritParams .muRefFit
#' @return `list(theta=, coef=, eta=, fit=)`
#' @author Matthew L. Fidler
#' @noRd
.muRefLin <- function(phi, cov, fixedCoef = NULL) {
  .muRefFit(phi, cov, fixedCoef = fixedCoef, weights = NULL)
}

#' Reweighted mu-ref regression step (`muModel="irls"`)
#'
#' @inheritParams .muRefFit
#' @return `list(theta=, coef=, eta=, fit=)`
#' @author Matthew L. Fidler
#' @noRd
.muRefIrls <- function(phi, cov, fixedCoef = NULL, weights = NULL) {
  .muRefFit(phi, cov, fixedCoef = fixedCoef, weights = weights)
}
