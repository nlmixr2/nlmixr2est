# vaeData.R -- VAE data preparation. Builds, from the rxode2 ui + data:
#  - encoder inputs: standardized (time, DV) sequences padded to Tmax + lengths
#  - per-subject decoder inputs: event table, obs times, observed DV
#  - eta<->theta (z_pop) mapping and initial z_pop / omega / a from ini()
# Encoder-input standardization follows the reference (time/max, (DV-mean)/sd).

# Data columns never treated as covariate candidates by the VAE search
.vaeReservedCols <- c("ID", "TIME", "DV", "EVID", "AMT", "CMT", "MDV", "SS", "II",
                      "ADDL", "RATE", "DUR", "DVID", "CENS", "LIMIT", "OCC")

# Finite fallback bound (+/-) for an unbounded model-declared covariate coefficient
# regressed by the VAE M-step; keeps the 1-D optimize() step from running away.
# Generous on a transformed/log scale (an effect past this is implausible); a user
# ini() bound overrides it.
.vaeCovCoefBound <- 10

# Max plausible log-scale effect used to derive a SCALE-AWARE fallback bound for a
# raw linear covariate coefficient (beta*COV): |beta*max|COV|| <= this, so the
# bound shrinks as the covariate magnitude grows (a raw WT coefficient is ~1/WT).
.vaeCovCoefEffect <- 5

#' Scale-aware finite fallback bounds for unbounded covariate coefficients.
#'
#' For a linear mu-ref covariate coefficient (in `muRefCovariateDataFrame`) the
#' sane magnitude scales like `1/max|COV|`, so a fixed wide bound makes the 1-D
#' optimize() overshoot a shallow interior optimum.  Scale that coefficient's
#' bound by the covariate's data magnitude; leave transformed/unrecognized
#' coefficients (already O(1) regressors) at the flat `.vaeCovCoefBound`.
#' @param ui rxode2 ui
#' @param data raw (un-normalized) modeling data
#' @param names covariate-coefficient theta names needing a bound
#' @return named numeric of positive half-widths (one per `names`)
#' @noRd
.vaeCovCoefBoundVec <- function(ui, data, names) {
  .b <- setNames(rep(.vaeCovCoefBound, length(names)), names)
  .mrc <- ui$muRefCovariateDataFrame
  if (is.null(.mrc) || nrow(.mrc) == 0L) return(.b)
  .dd <- as.data.frame(data); colnames(.dd) <- toupper(colnames(.dd))
  for (.nm in names) {
    .w <- which(.mrc$covariateParameter == .nm)
    if (length(.w) == 0L) next
    .cov <- toupper(.mrc$covariate[.w[1]])
    if (is.null(.dd[[.cov]])) next
    .mx <- max(abs(as.numeric(.dd[[.cov]])), na.rm = TRUE)
    if (is.finite(.mx) && .mx > 0) .b[.nm] <- min(.vaeCovCoefBound, .vaeCovCoefEffect / .mx)
  }
  .b
}

#' Discover + encode subject-level covariates for the VAE search
#' @param d normalized data.frame (upper-case column names)
#' @param ids unique subject ids in estimation order
#' @return list(covNames, covMat, covType, covPop, tvExcl)
#' @noRd
.vaeCovariateSearch <- function(d, ids) {
  N <- length(ids)
  ## auto-discover subject-level covariate candidates (constant within ID),
  ## excluding reserved data columns. Paper encoding: continuous -> log(v/mean),
  ## categorical (<=2 levels) -> centered.
  .cand <- setdiff(names(d), .vaeReservedCols)
  .isSubjConst <- vapply(.cand, function(nm) {
    all(vapply(ids, function(id) {
      v <- d[[nm]][d$ID == id]; length(unique(v)) == 1L
    }, logical(1))) && is.numeric(d[[nm]])
  }, logical(1))
  # The VAE covariate search absorbs covariates as subject-level (time-invariant)
  # effects, so a covariate that varies within a subject (time-varying) cannot be
  # searched.  Unlike saem/mu-focei (whose covariates are declared in the model,
  # detectable via .nlmixrTimeVaryingCovariates), the VAE search scans every
  # numeric data column, so time-varying ones are those that are not
  # subject-constant.  They are reported back (tvExcl) so callers can warn.
  .covNames <- .cand[.isSubjConst]
  .numCand <- .cand[vapply(.cand, function(nm) is.numeric(d[[nm]]), logical(1))]
  .tvExcl <- setdiff(.numCand, .covNames)
  .covVal <- vapply(.covNames, function(nm) {
    vapply(ids, function(id) d[[nm]][d$ID == id][1], numeric(1))
  }, numeric(N))
  if (length(.covNames) == 1L) .covVal <- matrix(.covVal, N, 1L, dimnames = list(NULL, .covNames))
  .covType <- character(length(.covNames)); .covPop <- numeric(length(.covNames))
  .covMat <- matrix(0, N, length(.covNames), dimnames = list(NULL, .covNames))
  for (j in seq_along(.covNames)) {
    v <- .covVal[, j]
    if (length(unique(v)) > 2L && all(v > 0)) {
      .covType[j] <- "continuous"; .covPop[j] <- mean(v); .covMat[, j] <- log(v / .covPop[j])
    } else {
      .covType[j] <- "categorical"; .covPop[j] <- mean(v); .covMat[, j] <- v - .covPop[j]
    }
  }
  list(covNames = .covNames, covMat = .covMat, covType = .covType,
       covPop = .covPop, tvExcl = .tvExcl)
}

#' Covariates explored by the VAE covariate search
#'
#' Returns the subject-level covariates that `nlmixr2(..., est = "vae")` would
#' explore during automated covariate selection, using the same discovery rules
#' as the fit: every non-reserved numeric data column that is constant within
#' each subject is a candidate; a candidate with more than two unique values
#' (all positive) is treated as continuous (encoded `log(value/mean)`),
#' anything else as categorical (mean-centered).  Time-varying numeric columns
#' cannot be searched and are excluded with a warning.
#'
#' @param data estimation dataset containing at least an `ID` column; column
#'   names are matched case-insensitively, as in the VAE fit
#' @param warn when `TRUE` (default) warn about time-varying numeric columns
#'   excluded from the search; when `FALSE` exclude them silently
#' @return a data frame with one row per explored covariate and columns
#'   `covariate` (upper-cased column name), `type` (`"continuous"` or
#'   `"categorical"`) and `center` (the population value the covariate is
#'   centered at); zero rows when no covariates qualify
#' @export
#' @author Matthew L. Fidler
#' @examples
#' d <- data.frame(id = rep(1:3, each = 2), time = rep(0:1, 3), dv = rnorm(6),
#'                 wt = rep(c(70, 80, 60), each = 2),
#'                 sex = rep(c(0, 1, 0), each = 2))
#' vaeCovariates(d)
vaeCovariates <- function(data, warn = TRUE) {
  checkmate::assertLogical(warn, len = 1, any.missing = FALSE)
  d <- as.data.frame(data)
  names(d) <- toupper(names(d))
  if (is.null(d$ID)) {
    stop("'data' must contain an ID column", call. = FALSE)
  }
  .cov <- .vaeCovariateSearch(d, unique(d$ID))
  if (warn && length(.cov$tvExcl) > 0L) {
    warning("time-varying covariate(s) were excluded from automatic covariate search: ",
            paste(.cov$tvExcl, collapse = ", "), call. = FALSE)
  }
  data.frame(covariate = .cov$covNames, type = .cov$covType,
             center = .cov$covPop, row.names = NULL)
}

#' Detect a clean `log(cov/center)` form for `cov` inside an expression.
#'
#' Walks the parse tree; returns the divisor `center` when `cov` occurs as
#' `log(cov)` (center 1) or `log(cov/<numeric literal>)`.  Any other form
#' (raw `cov`, `log(cov/expr)`, `log(a*cov)`) yields `inLog=FALSE` or a `NA`
#' center so the caller can fall back to the in-place regress M-step.
#' @noRd
.vaeLogCenter <- function(e, cov) {
  if (is.call(e)) {
    if (identical(e[[1L]], as.name("log")) && length(e) == 2L &&
          cov %in% all.vars(e[[2L]])) {
      .a <- e[[2L]]
      if (is.name(.a) && identical(as.character(.a), cov)) {
        return(list(inLog = TRUE, center = 1))
      }
      if (is.call(.a) && identical(.a[[1L]], as.name("/")) && length(.a) == 3L &&
            is.name(.a[[2L]]) && identical(as.character(.a[[2L]]), cov) &&
            is.numeric(.a[[3L]])) {
        return(list(inLog = TRUE, center = as.numeric(.a[[3L]])))
      }
      return(list(inLog = TRUE, center = NA_real_))
    }
    for (.i in seq_along(e)[-1L]) {
      .r <- .vaeLogCenter(e[[.i]], cov)
      if (isTRUE(.r$inLog)) return(.r)
    }
  }
  list(inLog = FALSE, center = NA_real_)
}

#' Model-declared covariate/parameter pairs for pinned VAE selection.
#'
#' One row per model-written covariate coefficient, resolving which latent dim
#' (eta) `k` it modifies, the data covariate, the user's coefficient theta, and
#' whether the pair can be handled by the restricted branch-and-bound search
#' (`inPool`) -- i.e. the covariate is in the subject-level search pool AND its
#' written functional form matches the VAE encoding so the estimated slope
#' transfers directly (continuous->`log(cov/center)`, categorical->linear).
#' Pairs that are not `inPool` (out-of-pool covariate, or a form whose slope
#' would not transfer) are estimated in place by the regress M-step instead.
#' @param ui rxode2 ui
#' @param covNames upper-cased search-pool covariate names (`prep$covNames`)
#' @param covType per-pool encoding (`"continuous"`/`"categorical"`)
#' @return data frame (k, covName, coefName, thetaName, covType, userCenter,
#'   inPool), or `NULL` when the model declares no covariate effects
#' @noRd
.vaeModelCovariatePairs <- function(ui, covNames, covType) {
  .coefThetas <- .vaeCovariateCoefThetas(ui)
  if (length(.coefThetas) == 0L) return(NULL)
  .thetaForEta <- .foceiEtaThetaMap(ui)$thetaForEta
  .thetaPool <- .thetaForEta[!is.na(.thetaForEta)]
  .allCov <- ui$allCovs
  if (is.null(.allCov)) .allCov <- character(0)
  .mrc <- ui$muRefCovariateDataFrame
  .lst <- ui$lstExpr
  .rows <- vector("list", 0L)
  for (.coef in .coefThetas) {
    .thName <- NA_character_; .covTok <- NA_character_; .linear <- FALSE
    if (!is.null(.mrc) && nrow(.mrc) > 0L && .coef %in% .mrc$covariateParameter) {
      .r <- .mrc[.mrc$covariateParameter == .coef, , drop = FALSE][1L, ]
      .thName <- as.character(.r$theta); .covTok <- as.character(.r$covariate)
      .linear <- TRUE
    } else {
      .lines <- Filter(function(e) .coef %in% all.vars(e), .lst)
      if (length(.lines) == 0L) next
      .vars <- all.vars(.lines[[1L]])
      .thHit <- intersect(.thetaPool, .vars)
      .covHit <- intersect(.allCov, .vars)
      if (length(.thHit) != 1L || length(.covHit) != 1L) next
      .thName <- .thHit; .covTok <- .covHit
    }
    .k <- match(.thName, .thetaForEta)
    if (is.na(.k)) next
    .j <- match(toupper(.covTok), covNames)
    .inPool <- !is.na(.j)
    .ct <- if (.inPool) covType[.j] else NA_character_
    .userCenter <- NA_real_
    if (.inPool) {
      if (identical(.ct, "categorical")) {
        ## VAE centers a categorical covariate; a plain linear beta*cov transfers
        ## (slope invariant to the shift).  A transformed categorical is unusual
        ## and not slope-transferable -- route it to the regress M-step.
        if (.linear) .userCenter <- 0 else .inPool <- FALSE
      } else {
        ## continuous: VAE uses log(cov/mean); only a written log(cov/center)
        ## transfers.  A raw linear beta*cov on a continuous covariate does not.
        .lc <- .vaeLogCenter(
          if (length(.lst)) Filter(function(e) .coef %in% all.vars(e), .lst)[[1L]] else NULL,
          .covTok)
        if (isTRUE(.lc$inLog) && is.finite(.lc$center)) .userCenter <- .lc$center
        else .inPool <- FALSE
      }
    }
    .rows[[length(.rows) + 1L]] <- data.frame(
      k = as.integer(.k),
      covName = if (.inPool) covNames[.j] else toupper(.covTok),
      coefName = .coef, thetaName = .thName,
      covType = if (is.na(.ct)) NA_character_ else .ct,
      userCenter = .userCenter, inPool = .inPool,
      stringsAsFactors = FALSE)
  }
  if (length(.rows) == 0L) return(NULL)
  do.call(rbind, .rows)
}

#' Prepare VAE inputs from a ui + data
#' @param ui rxode2 ui object
#' @param data estimation data (ID/TIME/DV/EVID/... columns)
#' @return list of prepared VAE inputs
#' @noRd
.vaeDataPrep <- function(ui, data, control = NULL) {
  .idf <- ui$iniDf
  .map <- .foceiEtaThetaMap(ui)
  .etaNames <- .map$etaNames
  .neta <- length(.etaNames)
  if (.neta == 0L) stop("est=\"vae\" requires at least one random effect", call. = FALSE)

  ## full theta vector (THETA_i_ in ntheta order), from ini estimates
  .thRows <- .idf[!is.na(.idf$ntheta), , drop = FALSE]
  .thRows <- .thRows[order(.thRows$ntheta), , drop = FALSE]
  .th <- setNames(as.numeric(.thRows$est), paste0("THETA_", seq_len(nrow(.thRows)), "_"))
  ## structural theta index (in the full theta vector) paired with each eta.
  ## A random effect that is not mu-referenced to a single theta -- a mixture eta
  ## (mix(exp(lke1+eta.ke),p,exp(lke2+eta.ke))), an eta on a fixed (literalFix-ed)
  ## theta, or a genuinely free eta -- is modeled as theta+eta with theta forced
  ## to 0: it centers at 0, is held there by the M-step, and is excluded from
  ## covariate selection; the rest of the model (literal / component thetas /
  ## covariate expression) carries its structure.
  .zPopThetaIdx <- match(.map$thetaForEta, .thRows$name)
  .isFree <- is.na(.zPopThetaIdx)
  .zPop <- numeric(.neta)                                      # structural population means (transformed)
  .zPop[!.isFree] <- as.numeric(.th[.zPopThetaIdx[!.isFree]])
  ## latent dims whose backing structural theta is FIXED (e.g. nonMuTheta="fix", or
  ## a user-fixed theta carrying an eta): the M-step holds their typical value at
  ## the ini() value, and they are dropped from the iteration print.
  .zPopFix <- logical(.neta)
  .zPopFix[!.isFree] <- isTRUE2(.thRows$fix[.zPopThetaIdx[!.isFree]])

  ## omega init (diagonal) for the etas + which variances are FIXED (held by the
  ## M-step, not estimated)
  .omega <- vapply(.etaNames, function(nm) {
    .r <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2 & .idf$name == nm, , drop = FALSE]
    as.numeric(.r$est[1])
  }, numeric(1))
  .omegaFix <- vapply(.etaNames, function(nm) {
    .r <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2 & .idf$name == nm, , drop = FALSE]
    isTRUE(as.logical(.r$fix[1]))
  }, logical(1))
  ## structural-theta bounds per eta (Inf/-Inf when unbounded or free): the M-step
  ## clamps the population estimate to [lower, upper], giving the constrained
  ## estimate (at the bound when the unconstrained optimum is outside).
  .zPopLower <- rep(-Inf, .neta); .zPopUpper <- rep(Inf, .neta)
  .zPopLower[!.isFree] <- as.numeric(.thRows$lower[.zPopThetaIdx[!.isFree]])
  .zPopUpper[!.isFree] <- as.numeric(.thRows$upper[.zPopThetaIdx[!.isFree]])

  ## normalize data columns (needed early for covariate discovery + pinning)
  d <- as.data.frame(data)
  names(d) <- toupper(names(d))
  ## no EVID: derive from AMT; with no AMT column either (dose-free data), all
  ## rows are observations (d$AMT is NULL -> the ifelse would yield length 0).
  if (is.null(d$EVID)) {
    d$EVID <- if (is.null(d$AMT)) rep(0L, nrow(d)) else ifelse(is.na(d$AMT) | d$AMT == 0, 0L, 1L)
  }
  .ids <- unique(d$ID)
  N <- length(.ids)

  ## subject-level covariate discovery + encoding (shared with vaeCovariates())
  .cov <- .vaeCovariateSearch(d, .ids)
  if (length(.cov$tvExcl) > 0L) {
    warning("time-varying covariate(s) were excluded from automatic covariate search: ",
            paste(.cov$tvExcl, collapse = ", "), call. = FALSE)
  }

  ## pinCovariates=FALSE with a model that declares covariates: turn OFF the
  ## automatic search and estimate every declared covariate in place by the
  ## regress M-step (the covariateSelection=FALSE treatment).  Emptying the
  ## search pool makes the C++ M-step skip covariate selection entirely.  With no
  ## model-declared covariates there is nothing to switch off -- the full search
  ## runs.  (Explicit covariateSelection=FALSE keeps its own path below.)
  .declaredCoefs <- .vaeCovariateCoefThetas(ui)
  .searchOff <- isFALSE(control$pinCovariates) && length(.declaredCoefs) > 0L &&
    !isFALSE(control$covariateSelection)
  if (.searchOff) {
    warning("pinCovariates=FALSE: model covariates estimated in place", call. = FALSE)
    .cov$covNames <- character(0)
    .cov$covMat <- matrix(0, N, 0L)
    .cov$covType <- character(0)
    .cov$covPop <- numeric(0)
  }

  ## pinned covariate selection: restrict the search to model-declared covariate
  ## /parameter pairs.  Build the per-(eta k x covariate j) allow-mask from the
  ## `inPool` declared pairs; zero those coefficients in the training theta so the
  ## decoder stays covariate-free (the effect is recovered by the M-step prior
  ## regression, exactly as in the unconstrained search) and injected back into
  ## the model afterward.  Declared pairs that cannot be searched are routed to
  ## the regress M-step below (`.pinCovCoef`).
  .pinActive <- FALSE
  .covAllow <- NULL
  .pinPairs <- NULL
  .pinCovCoef <- character(0)
  if (isTRUE(control$pinCovariates) && !isFALSE(control$covariateSelection)) {
    .pinPairs <- .vaeModelCovariatePairs(ui, .cov$covNames, .cov$covType)
    if (!is.null(.pinPairs) && nrow(.pinPairs) > 0L) {
      .pinActive <- TRUE
      .inRows <- .pinPairs[.pinPairs$inPool, , drop = FALSE]
      .nCov <- length(.cov$covNames)
      if (.nCov > 0L && nrow(.inRows) > 0L) {
        .covAllow <- matrix(0L, .neta, .nCov)
        for (.r in seq_len(nrow(.inRows))) {
          .j <- match(.inRows$covName[.r], .cov$covNames)
          if (!is.na(.j)) .covAllow[.inRows$k[.r], .j] <- 1L
        }
        ## decoder covariate-free during training: hold the declared (in-pool)
        ## coefficient at 0 so the covariate enters only through the prior.
        for (.cn in unique(.inRows$coefName)) {
          .ti <- match(.cn, .thRows$name)
          if (!is.na(.ti)) .th[.ti] <- 0
        }
      }
      .pinCovCoef <- unique(.pinPairs$coefName[!.pinPairs$inPool])
      warning("covariate selection pinned to model-specified covariates", call. = FALSE)
      if (length(.pinCovCoef) > 0L) {
        warning("pinned covariate(s) outside search pool estimated in place", call. = FALSE)
      }
    }
  }

  ## Fixed-effect thetas estimated directly by a bounded bobyqa regression in the
  ## M-step (vs. the latent-space zPop update).  Two sources, unioned:
  ##  * nonMuTheta="regress": non-mu-referenced structural thetas (no eta); and
  ##  * covariateSelection=FALSE: model-declared covariate coefficients -- always
  ##    estimated in place, independent of nonMuTheta, so the mu-referenced
  ##    covariate expression the user wrote is fit rather than held fixed.
  ## Carry each 0-based index into the full theta vector (`.th`, ntheta order)
  ## plus the ini() bounds (NA -> +-Inf).
  .regressNames <- character(0)
  .regressThetaIdx0 <- integer(0)
  .regressLower <- numeric(0); .regressUpper <- numeric(0)
  if (identical(control$nonMuTheta, "regress")) {
    .regressNames <- .vaeNonMuThetas(ui)
  }
  .covCoefNames <- character(0)
  if (isFALSE(control$covariateSelection)) {
    .covCoefNames <- .declaredCoefs
  } else if (.searchOff) {
    ## pinCovariates=FALSE with model-declared covariates: all of them regress.
    .covCoefNames <- .declaredCoefs
  } else if (.pinActive && length(.pinCovCoef) > 0L) {
    ## pinned selection: a declared covariate that cannot be handled by the
    ## restricted search (out-of-pool, or a form whose slope will not transfer)
    ## is estimated in place by the regress M-step, like covariateSelection=FALSE.
    .covCoefNames <- .pinCovCoef
  }
  .regressNames <- c(.regressNames, .covCoefNames)
  .regressNames <- unique(.regressNames)
  if (length(.regressNames) > 0L) {
    .ri <- match(.regressNames, .thRows$name)
    .regressThetaIdx0 <- as.integer(.ri - 1L)
    .lo <- as.numeric(.thRows$lower[.ri]); .hi <- as.numeric(.thRows$upper[.ri])
    .regressLower <- ifelse(is.na(.lo), -Inf, .lo)
    .regressUpper <- ifelse(is.na(.hi), Inf, .hi)
    ## An UNBOUNDED covariate coefficient regressed alone routes through the 1-D
    ## optimize() branch of .boundedResidOpt, which searches the whole interval and
    ## overshoots a shallow interior optimum on a too-wide interval.  Give an
    ## unbounded coefficient a finite, scale-aware fallback interval (user ini()
    ## bounds still win); structural regress thetas keep their bounds (bobyqa is
    ## stable, so they are left alone).
    .isCov <- .regressNames %in% .covCoefNames
    .noLo <- .isCov & !is.finite(.regressLower)
    .noHi <- .isCov & !is.finite(.regressUpper)
    if (any(.noLo | .noHi)) {
      .bnd <- .vaeCovCoefBoundVec(ui, data, .regressNames[.isCov])
      .regressLower[.noLo] <- -.bnd[.regressNames[.noLo]]
      .regressUpper[.noHi] <- .bnd[.regressNames[.noHi]]
    }
  }

  ## residual error params (all of them, in theta order): value, theta index,
  ## type (add/prop/...), and bounds. Combined models have >1 row; log-likelihood
  ## models may have none. `a` is the (named) error-param vector.
  .errRow <- .idf[!is.na(.idf$err) & !is.na(.idf$ntheta), , drop = FALSE]
  .errRow <- .errRow[order(.errRow$ntheta), , drop = FALSE]
  .a <- if (nrow(.errRow) > 0) setNames(as.numeric(.errRow$est), .errRow$name) else numeric(0)
  .errThetaIdx <- as.integer(.errRow$ntheta)
  .errType <- as.character(.errRow$err)
  .errLower <- as.numeric(.errRow$lower); .errUpper <- as.numeric(.errRow$upper)

  ## per-subject decoder inputs + gather all obs for standardization
  subj <- vector("list", N)
  .allTime <- numeric(0); .allDv <- numeric(0)
  for (i in seq_len(N)) {
    .di <- d[d$ID == .ids[i], , drop = FALSE]
    .obs <- .di[.di$EVID == 0, , drop = FALSE]
    .times <- .obs$TIME
    .y <- .obs$DV
    ## M2/M3/M4 censoring columns (0 / NA when absent)
    .cens <- if (is.null(.obs$CENS)) integer(length(.y)) else as.integer(.obs$CENS)
    .limit <- if (is.null(.obs$LIMIT)) rep(NA_real_, length(.y)) else as.numeric(.obs$LIMIT)
    subj[[i]] <- list(ev = .di, times = .times, y = .y, n = length(.times),
                      cens = .cens, limit = .limit)
    .allTime <- c(.allTime, .times); .allDv <- c(.allDv, .y)
  }
  .tMax <- max(.allTime); .dvMean <- mean(.allDv); .dvSd <- stats::sd(.allDv)

  ## encoder inputs: [N, Tmax, 2] standardized (time, DV), padded; lengths
  Tmax <- max(vapply(subj, function(s) s$n, integer(1)))
  dataIn <- array(0, c(N, Tmax, 2L))
  lengths <- integer(N)
  for (i in seq_len(N)) {
    s <- subj[[i]]; ni <- s$n; lengths[i] <- ni
    dataIn[i, seq_len(ni), 1L] <- s$times / .tMax
    dataIn[i, seq_len(ni), 2L] <- (s$y - .dvMean) / .dvSd
  }
  covIn <- matrix(0, N, 0L)                     # encoder-head covariates (unused for now)

  list(N = N, neta = .neta, zDim = .neta, etaNames = .etaNames,
       th = .th, zPopThetaIdx = .zPopThetaIdx, isFree = .isFree, omegaFix = .omegaFix,
       zPopFix = .zPopFix,
       zPopLower = .zPopLower, zPopUpper = .zPopUpper,
       errThetaIdx = .errThetaIdx, errType = .errType,
       errLower = .errLower, errUpper = .errUpper,
       regressNames = .regressNames, regressThetaIdx0 = .regressThetaIdx0,
       regressLower = .regressLower, regressUpper = .regressUpper,
       zPop = .zPop, omega = .omega, a = .a,
       subj = subj, dataIn = dataIn, lengths = lengths, covIn = covIn,
       covNames = .cov$covNames, covMat = .cov$covMat, covType = .cov$covType,
       covPop = .cov$covPop,
       pinActive = .pinActive, pinPairs = .pinPairs, covAllow = .covAllow,
       tMax = .tMax, dvMean = .dvMean, dvSd = .dvSd, Nobs = length(.allDv))
}
