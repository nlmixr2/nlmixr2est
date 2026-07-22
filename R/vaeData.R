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

# Fallback bound half-width for an UNBOUNDED structural non-mu theta estimated by
# the M-step: `ini() estimate +/- max(.vaeNonMuThetaBound, |est| * .vaeNonMuThetaRel)`.
# Purely a divergence guard -- with +-Inf bounds a flat direction runs the estimate
# away (~1e68 on an unbounded theo_sd `tv`).  Deliberately generous so it does not
# bind at a sane optimum; a user `ini()` bound always wins.
.vaeNonMuThetaBound <- 10
.vaeNonMuThetaRel <- 3

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
    ## mu2/mu3 covariates are pre-transformed by the hook into a linear
    ## nlmixrMuDerCov# data column (the centering/transform is already baked in),
    ## so encode them LINEARLY (mean-centered) -- never re-apply a log transform.
    .isMuDer <- grepl("^NLMIXRMUDERCOV[0-9]+$", .covNames[j], ignore.case = TRUE)
    ## a 0/1 indicator column (e.g. SEXF) is already in its natural
    ## parameterization: leave it RAW so its coefficient is the level-1 shift and
    ## the structural theta stays the reference (0) value.  (`%in%` yields FALSE
    ## for NA, so this is already a strict TRUE/FALSE; isTRUE makes that explicit.)
    .isInd <- isTRUE(all(v %in% c(0, 1)))
    if (!.isMuDer && !.isInd && length(unique(v)) > 2L && all(v > 0)) {
      .covType[j] <- "continuous"; .covPop[j] <- mean(v); .covMat[, j] <- log(v / .covPop[j])
    } else if (.isInd) {
      .covType[j] <- "categorical"; .covPop[j] <- 0; .covMat[, j] <- v
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
            is.numeric(.a[[3L]]) && length(.a[[3L]]) == 1L && is.finite(.a[[3L]])) {
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

#' Covariate that a coefficient multiplies, within an expression.
#'
#' Disambiguates which data covariate `coef` pairs with when a model line carries
#' several covariate effects (e.g. `wt.cl*log(WT/70) + sex.cl*SEX`): walks to the
#' `*` term containing `coef` and returns the single covariate on the other side.
#' `NULL` when it cannot be resolved to exactly one covariate.
#' @noRd
.vaeCoefCov <- function(e, coef, covs) {
  if (is.call(e)) {
    if (identical(e[[1L]], as.name("*")) && length(e) == 3L) {
      .lv <- all.vars(e[[2L]]); .rv <- all.vars(e[[3L]])
      if (coef %in% .lv) { .c <- intersect(.rv, covs); if (length(.c) == 1L) return(.c) }
      if (coef %in% .rv) { .c <- intersect(.lv, covs); if (length(.c) == 1L) return(.c) }
    }
    for (.i in seq_along(e)[-1L]) {
      .r <- .vaeCoefCov(e[[.i]], coef, covs)
      if (!is.null(.r)) return(.r)
    }
  }
  NULL
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
      .thName <- as.character(.r$theta)
      ## rxode2 may record an ALGEBRAIC covariate expression here (mu2-style,
      ## e.g. "log(0.0142857 * WT)") rather than a bare data column.  Only take
      ## it as a plain linear effect when it names a pool covariate directly;
      ## otherwise fall through to the model-line scan below.
      .cand <- as.character(.r$covariate)
      if (!is.na(match(toupper(.cand), covNames))) {
        .covTok <- .cand
        .linear <- TRUE
      }
    }
    if (is.na(.covTok)) {
      .lines <- Filter(function(e) .coef %in% all.vars(e), .lst)
      if (length(.lines) == 0L) next
      .vars <- all.vars(.lines[[1L]])
      if (is.na(.thName)) {
        .thHit <- intersect(.thetaPool, .vars)
        if (length(.thHit) == 1L) .thName <- .thHit
      }
      ## a line may carry several covariate effects (e.g.
      ## wt.cl*log(WT/70) + sex.cl*SEX): pick the covariate THIS coefficient
      ## multiplies rather than skipping the coefficient (skipping could drop
      ## pinning to the unrestricted full search).  Unresolved -> not pinnable.
      .cc <- .vaeCoefCov(.lines[[1L]], .coef, .allCov)
      if (!is.null(.cc)) .covTok <- .cc
    }
    ## Always emit a row for a detected coefficient so pinning stays restrictive;
    ## a pair that cannot be resolved/transferred is marked not `inPool` and
    ## estimated in place by the regress M-step.
    .k <- match(.thName, .thetaForEta)
    .j <- if (!is.na(.covTok)) match(toupper(.covTok), covNames) else NA_integer_
    .inPool <- !is.na(.k) && !is.na(.j)
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
        .cl <- Filter(function(e) .coef %in% all.vars(e), .lst)
        .lc <- if (length(.cl)) .vaeLogCenter(.cl[[1L]], .covTok)
               else list(inLog = FALSE, center = NA_real_)
        if (isTRUE(.lc$inLog) && is.finite(.lc$center)) .userCenter <- .lc$center
        else .inPool <- FALSE
      }
    }
    .rows[[length(.rows) + 1L]] <- data.frame(
      k = if (is.na(.k)) NA_integer_ else as.integer(.k),
      covName = if (.inPool) covNames[.j] else if (is.na(.covTok)) NA_character_ else toupper(.covTok),
      coefName = .coef, thetaName = if (is.na(.thName)) NA_character_ else .thName,
      covType = if (is.na(.ct)) NA_character_ else .ct,
      userCenter = .userCenter, inPool = .inPool,
      stringsAsFactors = FALSE)
  }
  if (length(.rows) == 0L) return(NULL)
  do.call(rbind, .rows)
}

#' Which regressed thetas may be optimized in `residOptimize="twoStage"` stage 2.
#'
#' Stage 2 pins the ODE states from stage 1 and re-optimizes with only the
#' candidate moving, so a parameter is eligible exactly when the state trajectory
#' cannot depend on it.  The rule is PER PARAMETER:
#'
#'   eligible = it is an `err` parameter, OR no solve-defining expression reaches it
#'
#' The `err` half is the historic rule and keeps every error parameter in stage 2
#' exactly as before.  The second half is what an `ll()`/generalized endpoint
#' needs: its residual-like parameters are plain thetas with no `err` row, so on
#' the old rule stage 2 was empty and `"twoStage"` silently degraded to the joint
#' `"optimize"` solve.
#'
#' A model-level "does this model have error parameters" short-circuit would NOT
#' do -- a multi-endpoint model with one Gaussian and one `ll()` endpoint has
#' `err` rows AND log-density-only thetas, and both belong in stage 2.
#'
#' Deliberately conservative in the risky direction: an unparsable model, or one
#' with no ODE states to pin, contributes nothing beyond the `err` set.
#' @param ui rxode2 ui object
#' @param regressNames regressed theta names, in `regressThetaIdx0` order
#' @param regressErrIdx0 0-based slot in `a` per regressed name, -1 when not one
#' @return integer 0/1 per regressed name
#' @noRd
.vaeRegressStage2 <- function(ui, regressNames, regressErrIdx0) {
  if (length(regressNames) == 0L) return(integer(0))
  .isErr <- regressErrIdx0 >= 0L
  .odeFree <- tryCatch(.vaeOdeFreeThetas(ui, regressNames),
                       error = function(e) rep(FALSE, length(regressNames)))
  as.integer(.isErr | .odeFree)
}

#' Names the ODE state trajectory provably cannot depend on.
#'
#' Fixpoint over the model assignments, seeded from every expression that feeds
#' the solve.  Assignment ORDER is ignored, which can only over-collect symbols
#' and therefore only ever answers FALSE where the truth is TRUE -- the safe
#' direction (the theta stays in stage 1, as it is today).
#' @param ui rxode2 ui object
#' @param names candidate names
#' @return logical vector, `TRUE` when no solve-defining expression reaches the name
#' @noRd
.vaeOdeFreeThetas <- function(ui, names) {
  .no <- rep(FALSE, length(names))
  .lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  if (is.null(.lst) || length(.lst) == 0L) return(.no)
  if (length(tryCatch(rxode2::rxState(ui), error = function(e) character(0))) == 0L) return(.no)
  .isAssign <- function(.ex) is.call(.ex) && length(.ex) == 3L &&
    (identical(.ex[[1]], as.name("<-")) || identical(.ex[[1]], as.name("=")) ||
       identical(.ex[[1]], as.name("~")))
  .syms <- function(.e) {
    if (is.name(.e)) return(as.character(.e))
    if (is.call(.e)) return(unlist(lapply(as.list(.e)[-1L], .syms), use.names = FALSE))
    character(0)
  }
  ## a solve-defining left-hand side: d/dt(x), x(0), and the dosing modifiers
  .solveLhs <- function(.txt) {
    grepl("^d */ *dt *\\(", .txt) ||
      grepl("^[A-Za-z._][A-Za-z0-9._]* *\\( *0 *\\)$", .txt) ||
      grepl("^(f|alag|lag|rate|dur) *\\(", .txt)
  }
  ## the observation model, which the SOLVE does not read: `ll(x) ~ <density>`
  ## and the error-model endpoint lines (`x ~ add(...)`, `x ~ pois(...)`).  ONLY a
  ## `~` line qualifies -- `cp <- center / v` defines a variable a `d/dt()` may
  ## well read, and dropping it would hide a real dependency.
  .endpointVars <- tryCatch(as.character(ui$predDf$var), error = function(e) character(0))
  .seed <- character(0); .map <- list()
  for (.ex in .lst) {
    if (!.isAssign(.ex)) next
    .rhs <- .syms(.ex[[3]])
    .txt <- paste(deparse(.ex[[2]]), collapse = "")
    .tilde <- identical(.ex[[1]], as.name("~"))
    if (.solveLhs(.txt)) {
      .seed <- c(.seed, .rhs)
    } else if (.tilde && (grepl("^ll *\\(", .txt) || .txt %in% .endpointVars)) {
      next                                   # likelihood only -- never reaches the solve
    } else if (is.name(.ex[[2]])) {
      .map[[.txt]] <- unique(c(.map[[.txt]], .rhs))
    } else {
      ## an lhs shape not recognized above may still feed the solve -- take its
      ## rhs as reachable rather than guess (the safe direction: stage 1)
      .seed <- c(.seed, .rhs)
    }
  }
  .seen <- unique(.seed); .todo <- .seen
  while (length(.todo) > 0L) {
    .nxt <- unique(unlist(.map[intersect(.todo, names(.map))], use.names = FALSE))
    .todo <- setdiff(.nxt, .seen)
    .seen <- c(.seen, .todo)
  }
  !(as.character(names) %in% .seen)
}

#' Prepare VAE inputs from a ui + data
#' @param ui rxode2 ui object
#' @param data estimation data (ID/TIME/DV/EVID/... columns)
#' @return list of prepared VAE inputs
#' @noRd
.vaeDataPrep <- function(ui, data, control = NULL) {
  .inputScale <- if (is.null(control$inputScale)) "reference" else control$inputScale
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
    ## keep the $runInfo note single-line even with many covariates
    .tvPre <- "time-varying covariate(s) not searched: "
    warning(.tvPre, .vaeTruncList(.cov$tvExcl, prefix = .tvPre), call. = FALSE)
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
      .nCov <- length(.cov$covNames)
      if (.nCov > 0L) {
        ## A covariate column carries ONE encoding.  Claim each column for the
        ## first declared pair's center; if the same covariate is declared again
        ## with a DIFFERENT center (e.g. log(WT/70) on CL and log(WT/80) on KA)
        ## that pair cannot share the column, so demote it to the regress M-step.
        .claim <- rep(NA_real_, .nCov)
        for (.r in seq_len(nrow(.pinPairs))) {
          if (!.pinPairs$inPool[.r]) next
          .j <- match(.pinPairs$covName[.r], .cov$covNames)
          if (is.na(.j)) {
            .pinPairs$inPool[.r] <- FALSE
          } else if (is.na(.claim[.j])) {
            .claim[.j] <- .pinPairs$userCenter[.r]
          } else if (!isTRUE(all.equal(.claim[.j], .pinPairs$userCenter[.r]))) {
            .pinPairs$inPool[.r] <- FALSE
          }
        }
        .inRows <- .pinPairs[.pinPairs$inPool, , drop = FALSE]
        ## restrict the search to the declared in-pool cells.  An all-zero row
        ## means "no covariate may be selected on this dim" -- crucial when every
        ## declared pair is out-of-pool, so a non-declared (or the out-of-pool)
        ## covariate is never auto-selected under pinning.
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
        ## Retain ONLY the model's own centering (already carried by the mu2/mu3
        ## nlmixrMuDerCov# column, or by the written log(cov/center)); do NOT add
        ## the VAE's mean-centering.  Use each pinned covariate at its MODEL value
        ## so zPop is the model intercept and no post-hoc correction is needed.
        ## Centering a predictor in a regression WITH an intercept leaves the slope
        ## and the selection unchanged -- this only relocates the intercept.
        ## Each claimed column is adjusted EXACTLY once, off the original covPop.
        for (.j in which(!is.na(.claim))) {
          if (identical(.cov$covType[.j], "continuous")) {
            ## log(v/mean) -> log(v/userCenter)
            .cov$covMat[, .j] <- .cov$covMat[, .j] + log(.cov$covPop[.j]) - log(.claim[.j])
          } else {
            ## (v - mean) -> raw v (mu2/mu3 already applied the model transform)
            .cov$covMat[, .j] <- .cov$covMat[, .j] + .cov$covPop[.j]
          }
          .cov$covPop[.j] <- 0
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
  ##  * nonMuTheta="regress"/"grad": non-mu-referenced structural thetas (no eta); and
  ##  * covariateSelection=FALSE: model-declared covariate coefficients -- always
  ##    estimated in place, independent of nonMuTheta, so the mu-referenced
  ##    covariate expression the user wrote is fit rather than held fixed.
  ## Carry each 0-based index into the full theta vector (`.th`, ntheta order)
  ## plus the ini() bounds (NA -> +-Inf).
  .regressNames <- character(0)
  .regressThetaIdx0 <- integer(0)
  .regressLower <- numeric(0); .regressUpper <- numeric(0)
  if (.vaeNonMuIsRegress(control$nonMuTheta)) {
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
  ## residOptimize="optimize": the residual-error thetas join the SAME optimizer
  ## as the non-mu structural thetas, against the same full outer objective.  A
  ## FIXED error parameter is excluded (nothing to estimate), as is one already
  ## regressed for another reason.
  .errRegressNames <- character(0)
  if (!identical(control$residOptimize, "moment")) {
    .errAll <- .idf[!is.na(.idf$err) & !is.na(.idf$ntheta), , drop = FALSE]
    if (nrow(.errAll) > 0L) {
      .errFree <- .errAll[!(!is.na(.errAll$fix) & .errAll$fix), , drop = FALSE]
      ## Every free residual parameter enters the optimizer.  The stage-2
      ## objective evaluates the ordinary likelihood with the ODE frozen, so `r`
      ## comes from the model's own rx_r_ and ANY error form is scored correctly
      ## -- there is nothing left for a form-specific filter to protect against.
      ## (It was needed only while the objective recomputed `r` itself from a
      ## hardcoded per-form expression, which could silently ignore a parameter
      ## and let the optimizer move it on noise.)
      .errRegressNames <- setdiff(as.character(.errFree$name), .regressNames)
    }
  }
  .regressNames <- c(.regressNames, .errRegressNames)
  if (length(.regressNames) > 0L) {
    .ri <- match(.regressNames, .thRows$name)
    .regressThetaIdx0 <- as.integer(.ri - 1L)
    .lo <- as.numeric(.thRows$lower[.ri]); .hi <- as.numeric(.thRows$upper[.ri])
    .regressLower <- ifelse(is.na(.lo), -Inf, .lo)
    .regressUpper <- ifelse(is.na(.hi), Inf, .hi)
    ## A residual SCALE parameter must not be allowed to reach zero.  The
    ## likelihood floors a zero variance (r == 0 -> r = 1) to stay finite, which
    ## makes a collapsed residual look attractive rather than forbidden -- a
    ## boxCox fit converged to add.err = 0 exactly this way.  Floor the scale
    ## parameters strictly above zero: absolute forms (add, lnorm) relative to
    ## the spread of the data, relative forms (prop, pow) at a small constant.
    ## Exponents and lambdas are NOT scales and are left alone.
    .errScaleAbs <- as.character(.idf$name[!is.na(.idf$err) & .idf$err %in% c("add", "lnorm")])
    .errScaleRel <- as.character(.idf$name[!is.na(.idf$err) & .idf$err %in% c("prop", "pow")])
    .dvObs <- suppressWarnings(as.numeric(d$DV[d$EVID == 0]))
    .dvSpread <- stats::sd(.dvObs[is.finite(.dvObs)])
    if (!is.finite(.dvSpread) || .dvSpread <= 0) .dvSpread <- 1
    .isAbs <- .regressNames %in% .errScaleAbs
    .isRel <- .regressNames %in% .errScaleRel
    .regressLower[.isAbs] <- pmax(.regressLower[.isAbs], 1e-4 * .dvSpread)
    .regressLower[.isRel] <- pmax(.regressLower[.isRel], 1e-6)
    ## A transform-both-sides lambda (Box-Cox / Yeo-Johnson) is only meaningful
    ## on a narrow interval, and it is unbounded in the ini() block, so the
    ## optimizer would otherwise search a meaningless range.  Constrain it to
    ## (-2, 2); SAEM does the same thing by mapping lambda through a bounded
    ## transform (`toLambda`).  A tighter user bound still wins.
    .lamNames <- as.character(.idf$name[!is.na(.idf$err) &
                                        .idf$err %in% c("boxCox", "yeoJohnson")])
    if (length(.lamNames) > 0L) {
      .isLam <- .regressNames %in% .lamNames
      .regressLower[.isLam] <- pmax(.regressLower[.isLam], -2)
      .regressUpper[.isLam] <- pmin(.regressUpper[.isLam], 2)
    }
    ## An UNBOUNDED covariate coefficient regressed alone routes through the 1-D
    ## optimize() branch of .boundedResidOpt, which searches the whole interval and
    ## overshoots a shallow interior optimum on a too-wide interval.  Give an
    ## unbounded coefficient a finite, scale-aware fallback interval (user ini()
    ## bounds still win).
    .isCov <- .regressNames %in% .covCoefNames
    .noLo <- .isCov & !is.finite(.regressLower)
    .noHi <- .isCov & !is.finite(.regressUpper)
    if (any(.noLo | .noHi)) {
      .bnd <- .vaeCovCoefBoundVec(ui, data, .regressNames[.isCov])
      .regressLower[.noLo] <- -.bnd[.regressNames[.noLo]]
      .regressUpper[.noHi] <- .bnd[.regressNames[.noHi]]
    }
    ## A STRUCTURAL non-mu theta needs the same guard.  With +-Inf bounds nothing
    ## constrains the M-step (bobyqa's interval, or the "grad" Adam projection), and
    ## a theta whose likelihood is flat in one direction runs away: an unbounded
    ## `tv <- 3.45` on theo_sd reaches ~1e68 (an lnorm fit there reports an OFV of
    ## 359315 against focei's 686), while the same model with `tv <- c(2, 3.45, 5)`
    ## converges.  Fall back to a generous window around the ini() ESTIMATE, wide
    ## enough not to bind at a sane optimum but finite so the search cannot diverge.
    .isStruct <- !.isCov
    .sLo <- .isStruct & !is.finite(.regressLower)
    .sHi <- .isStruct & !is.finite(.regressUpper)
    if (any(.sLo | .sHi)) {
      .init <- as.numeric(.thRows$est[.ri])
      .init[!is.finite(.init)] <- 0
      ## scale-aware half-width: the absolute floor covers a log-scale parameter
      ## (init ~ 3.45 -> +-10 is exp(+-10), ample), the relative term keeps a
      ## large-magnitude natural-scale init (say 1000) from being over-constrained
      .hw <- pmax(.vaeNonMuThetaBound, abs(.init) * .vaeNonMuThetaRel)
      .regressLower[.sLo] <- (.init - .hw)[.sLo]
      .regressUpper[.sHi] <- (.init + .hw)[.sHi]
    }
  }

  ## residual error params (all of them, in theta order): value, theta index,
  ## type (add/prop/...), and bounds. Combined models have >1 row; log-likelihood
  ## models may have none. `a` is the (named) error-param vector.
  .errRow <- .idf[!is.na(.idf$err) & !is.na(.idf$ntheta), , drop = FALSE]
  .errRow <- .errRow[order(.errRow$ntheta), , drop = FALSE]
  .a <- if (nrow(.errRow) > 0) setNames(as.numeric(.errRow$est), .errRow$name) else numeric(0)
  ## For each regressed parameter, its 0-based slot in `a` (the error-parameter
  ## vector), or -1 when it is not an error parameter.  vaeBuildTh writes `a`
  ## over the error theta positions, so the objective must substitute the
  ## CANDIDATE value into `a`; without this map it would be flat in every error
  ## parameter and the optimizer would never move one.
  .regressErrIdx0 <- if (length(.regressNames) > 0L) {
    .m <- match(.regressNames, names(.a))
    as.integer(ifelse(is.na(.m), 0L, .m) - 1L)
  } else integer(0)
  .errThetaIdx <- as.integer(.errRow$ntheta)
  .errType <- as.character(.errRow$err)
  .errLower <- as.numeric(.errRow$lower); .errUpper <- as.numeric(.errRow$upper)
  ## residOptimize="twoStage" stage-2 eligibility (see .vaeRegressStage2)
  .regressStage2 <- .vaeRegressStage2(ui, .regressNames, .regressErrIdx0)

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
  Tmax <- max(vapply(subj, function(s) s$n, integer(1)))
  .tMax <- max(.allTime)
  ## inputScale: which DV values the encoder-input centering/scaling is computed
  ## over.  The reference takes mean/sd across the WHOLE padded [N, Tmax] matrix,
  ## so the zero padding of short subjects enters both -- on a ragged dataset
  ## that is a materially different scale from the observed-only statistics
  ## (neonatal: sd 1582 vs 506).  "reference" reproduces it; "observed" uses the
  ## observed values only.
  if (identical(.inputScale, "reference")) {
    .padded <- c(.allDv, rep(0, N * Tmax - length(.allDv)))
    .dvMean <- mean(.padded); .dvSd <- stats::sd(.padded)
  } else {
    .dvMean <- mean(.allDv); .dvSd <- stats::sd(.allDv)
  }
  if (!is.finite(.dvSd) || .dvSd <= 0) .dvSd <- 1

  ## encoder inputs: [N, Tmax, 2] standardized (time, DV), padded; lengths
  dataIn <- array(0, c(N, Tmax, 2L))
  lengths <- integer(N)
  for (i in seq_len(N)) {
    s <- subj[[i]]; ni <- s$n; lengths[i] <- ni
    dataIn[i, seq_len(ni), 1L] <- s$times / .tMax
    dataIn[i, seq_len(ni), 2L] <- (s$y - .dvMean) / .dvSd
  }
  ## Encoder-head covariates.  The reference concatenates them to the LSTM's
  ## FINAL HIDDEN STATE before the linear head that emits (mu, logSigma, L)
  ## -- `torch.cat((hidden[-1], covariates), dim=1)` in its encoder -- so the
  ## approximate posterior q(z|x) is conditioned on the covariates.  That is
  ## central to the method: it is how the encoder can express a covariate
  ## relationship at all, and the M-step then only has to read it off the
  ## posterior means.  Feeding zero columns here leaves the posterior
  ## unconditioned, which shows up as less between-subject spread aligned with
  ## the covariates (smaller omega, more variance pushed into residual error) and
  ## weaker covariate effects.  Use the same encoded matrix the selection step
  ## uses (continuous -> log(cov/mean), categorical -> linear).
  ## Only the covariates actually in play are supplied.  With no covariate
  ## search there is nothing for the encoder to condition on (the declared
  ## coefficients are estimated in place by the regress M-step instead), and
  ## under `pinCovariates` only the pinned candidates are; conditioning on a
  ## covariate the search cannot select would let the posterior encode a
  ## relationship the model never reports.
  covIn <- if (isFALSE(control$covariateSelection) || .searchOff) {
    matrix(0, N, 0L)
  } else if (!is.null(.covAllow) && ncol(.cov$covMat) > 0L) {
    .keep <- which(colSums(.covAllow) > 0L)
    .cov$covMat[, .keep, drop = FALSE]
  } else {
    .cov$covMat
  }
  if (!is.matrix(covIn) || nrow(covIn) != N) covIn <- matrix(0, N, 0L)

  list(N = N, neta = .neta, zDim = .neta, etaNames = .etaNames,
       th = .th, zPopThetaIdx = .zPopThetaIdx, isFree = .isFree, omegaFix = .omegaFix,
       zPopFix = .zPopFix,
       zPopLower = .zPopLower, zPopUpper = .zPopUpper,
       errThetaIdx = .errThetaIdx, errType = .errType,
       errLower = .errLower, errUpper = .errUpper,
       regressNames = .regressNames, regressThetaIdx0 = .regressThetaIdx0,
       regressErrIdx0 = .regressErrIdx0, regressStage2 = .regressStage2,
       regressLower = .regressLower, regressUpper = .regressUpper,
       zPop = .zPop, omega = .omega, a = .a,
       subj = subj, dataIn = dataIn, lengths = lengths, covIn = covIn,
       covNames = .cov$covNames, covMat = .cov$covMat, covType = .cov$covType,
       covPop = .cov$covPop,
       pinActive = .pinActive, pinPairs = .pinPairs, covAllow = .covAllow,
       tMax = .tMax, dvMean = .dvMean, dvSd = .dvSd, Nobs = length(.allDv))
}
