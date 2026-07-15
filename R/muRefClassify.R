#' Classify theta/eta pairs for the mu-referenced FOCEI family
#'
#' Thetas/etas with a mu-ref covariate relationship
#' (`ui$muRefCovariateDataFrame`) get the in-C++ regression machinery
#' (`updateMuGroups()`, `src/inner.cpp`). Plain (covariate-free) mu-ref
#' theta+eta pairs are also eligible, as intercept-only regression groups
#' (`muPlainThetas`/`muPlainEtas`), when the theta<->eta mapping is a clean
#' bijection and the eta enters the model in a single additive position.
#'
#' @param ui rxode2 ui object (post mu2-hook, if applicable)
#' @return A list with:
#' \itemize{
#'  \item{\code{muCovThetas}}{population thetas with at least one mu-ref
#'    covariate relationship}
#'  \item{\code{muCovEtas}}{etas whose associated theta is in
#'    \code{muCovThetas}}
#'  \item{\code{muCovCovariateParams}}{the covariate-coefficient thetas
#'    themselves (the "slopes")}
#'  \item{\code{muPlainThetas}}{plain mu-ref population thetas eligible for
#'    an intercept-only regression group}
#'  \item{\code{muPlainEtas}}{their etas, parallel to \code{muPlainThetas}}
#'  \item{\code{standardThetas}}{population thetas left untouched by the
#'    covariate machinery; excludes sigma/residual-error parameters}
#'  \item{\code{standardEtas}}{etas left untouched by the covariate
#'    machinery}
#' }
#' @author Matthew L. Fidler
#' @noRd
.muRefClassify <- function(ui) {
  .muRefDf <- ui$muRefDataFrame
  .muRefCovDf <- ui$muRefCovariateDataFrame
  .nonMuEtas <- ui$nonMuEtas
  .iniDf <- ui$iniDf

  if (is.null(.muRefDf) || nrow(.muRefDf) == 0L) {
    .muRefDf <- data.frame(theta=character(0), eta=character(0),
                            level=character(0), stringsAsFactors=FALSE)
  }
  if (is.null(.muRefCovDf) || nrow(.muRefCovDf) == 0L) {
    .muRefCovDf <- data.frame(theta=character(0), covariate=character(0),
                               covariateParameter=character(0),
                               stringsAsFactors=FALSE)
  }

  .muCovThetas <- unique(as.character(.muRefCovDf$theta))
  .muCovCovariateParams <- unique(as.character(.muRefCovDf$covariateParameter))

  .muCovRows <- .muRefDf[as.character(.muRefDf$theta) %in% .muCovThetas, , drop=FALSE]
  .muCovEtas <- unique(as.character(.muCovRows$eta))

  .stdMuRefRows <- .muRefDf[!(as.character(.muRefDf$theta) %in% .muCovThetas), , drop=FALSE]
  .standardEtas <- unique(as.character(.stdMuRefRows$eta))
  if (!is.null(.nonMuEtas) && length(.nonMuEtas) > 0L) {
    .standardEtas <- unique(c(.standardEtas, as.character(.nonMuEtas)))
  }

  .allThetas <- character(0)
  if (!is.null(.iniDf)) {
    .thetaRows <- .iniDf[!is.na(.iniDf$ntheta) & is.na(.iniDf$err), , drop=FALSE]
    .allThetas <- as.character(.thetaRows$name)
  }
  .standardThetas <- setdiff(.allThetas, c(.muCovThetas, .muCovCovariateParams))

  # Plain (covariate-free) mu-ref pairs eligible for an intercept-only
  # regression group: id-level, theta<->eta bijection in muRefDataFrame, a
  # diagonal non-IOV eta that enters the model in a single additive position
  # (.foceiEtaOccurrence; a shared eta makes the theta+eta rewrite invalid).
  if ("level" %in% names(.muRefDf)) {
    .idRows <- .muRefDf[is.na(.muRefDf$level) | .muRefDf$level == "id", , drop=FALSE]
  } else {
    .idRows <- .muRefDf
  }
  .iovVars <- as.character(.uiIovEnv$iovVars)
  .diagEtas <- character(0)
  if (!is.null(.iniDf)) {
    .e2 <- .iniDf[is.na(.iniDf$ntheta), , drop=FALSE]
    .diagEtas <- as.character(.e2$name[.e2$neta1 == .e2$neta2])
  }
  .thetaCount <- table(as.character(.muRefDf$theta))
  .etaCount <- table(as.character(.muRefDf$eta))
  .etaOcc <- .foceiEtaOccurrence(ui)
  .muPlainThetas <- character(0)
  .muPlainEtas <- character(0)
  for (.k in seq_len(nrow(.idRows))) {
    .th <- as.character(.idRows$theta[.k])
    .et <- as.character(.idRows$eta[.k])
    if (.th %in% .muCovThetas || .th %in% .muCovCovariateParams) next
    if (!(.th %in% .allThetas)) next
    if (.th %in% .iovVars || .et %in% .iovVars) next
    if (.thetaCount[[.th]] != 1L || .etaCount[[.et]] != 1L) next
    if (!(.et %in% .diagEtas)) next
    if (!is.null(.nonMuEtas) && .et %in% as.character(.nonMuEtas)) next
    if (!(.et %in% names(.etaOcc)) || .etaOcc[[.et]] != 1L) next
    .muPlainThetas <- c(.muPlainThetas, .th)
    .muPlainEtas <- c(.muPlainEtas, .et)
  }

  list(
    muCovThetas = .muCovThetas,
    muCovEtas = .muCovEtas,
    muCovCovariateParams = .muCovCovariateParams,
    muPlainThetas = .muPlainThetas,
    muPlainEtas = .muPlainEtas,
    standardThetas = .standardThetas,
    standardEtas = .standardEtas
  )
}

#' Build per-group theta/eta/covariate structures for the in-C++ regression
#' (`updateMuGroups()`, `src/inner.cpp`)
#'
#' Each group is one mu-ref-covariate population theta with its associated
#' eta (e.g. `cl <- exp(tcl + eta.cl + allo.cl*logWT)` groups
#' `theta="tcl"`/`eta="eta.cl"` with covariate `allo.cl` on `logWT`).
#' Covariate thetas with no associated eta are excluded (no per-subject
#' residual to regress against) and stay ordinary outer-optimized thetas.
#'
#' A user-fixed (`fix()`, `literalFix=FALSE`) population theta excludes the
#' group (the regression would move the fixed value); a user-fixed covariate
#' coefficient stays out of the regression (outer-supplied, not estimated by
#' it). Bounded parameters follow one of two styles, selected by `clamp`:
#'
#' - `clamp=TRUE` (the clamped mu-referenced FOCEI family: `mfocei`/
#'   `ifocei`/`mfoce`/`ifoce`/`magq`/`iagq`/`mlaplace`/`ilaplace`, i.e.
#'   `muModel != "none"`): bounded population thetas and covariate
#'   coefficients ARE grouped; the regression update is clamped to the
#'   bounds (box-constrained least squares, `updateMuGroups()`), and each
#'   group carries `thetaLower`/`thetaUpper` plus per-covariate
#'   `lower`/`upper` columns for the clamp.
#' - `clamp=FALSE` (the default for every other method): bounded parameters
#'   are rejected from the mu-referencing, since a plain (unclamped)
#'   regression cannot respect a box constraint. A finite bound on the
#'   group's population theta excludes the whole group (warns); a bound on
#'   a covariate coefficient excludes only that covariate (warns, marked
#'   `bounded=TRUE`, estimated as an ordinary bounded parameter by the
#'   outer optimizer) while the rest of the group keeps the speed-up.
#'
#' With `plain=TRUE`, plain (covariate-free) mu-ref pairs
#' (`muPlainThetas`, see `.muRefClassify`) are appended as intercept-only
#' groups (empty `covariates`); fixed plain thetas are excluded quietly, and
#' with `clamp=FALSE` bounded plain thetas are too (with an `.minfo`). If
#' profiling the plain groups would leave the outer optimizer with no free
#' parameter at all, the plain groups are dropped (covariate groups kept).
#'
#' @param ui rxode2 ui object (post mu2-hook, if applicable)
#' @param plain also form intercept-only groups for plain mu-ref pairs
#' @param clamp accept bounded parameters and carry their bounds for the
#'   clamped regression; defaults from the ui's `muModel` control so the
#'   clamped family clamps and everything else rejects
#' @return list of `list(theta=, eta=, thetaLower=, thetaUpper=,
#'   covariates=data.frame(covariate=, covariateParameter=, lower=,
#'   upper=, bounded=))`, one element per group whose population theta is
#'   not user-fixed
#' @author Matthew L. Fidler
#' @noRd
.muRefGroups <- function(ui, plain = FALSE,
                         clamp = !identical(rxode2::rxGetControl(ui, "muModel", "none"), "none")) {
  .cls <- .muRefClassify(ui)
  .muRefDf <- ui$muRefDataFrame
  .muRefCovDf <- ui$muRefCovariateDataFrame
  .iniDf <- ui$iniDf
  .thetaRows <- .iniDf[!is.na(.iniDf$ntheta), , drop = FALSE]
  .fixedThetas <- .thetaRows$name[which(.thetaRows$fix)]
  .boundedThetas <- .thetaRows$name[.thetaRows$lower > -Inf | .thetaRows$upper < Inf]
  .lowerOf <- setNames(.thetaRows$lower, .thetaRows$name)
  .upperOf <- setNames(.thetaRows$upper, .thetaRows$name)
  lapply(.cls$muCovThetas, function(.theta) {
    .w <- which(as.character(.muRefDf$theta) == .theta)
    if (length(.w) != 1L) return(NULL)
    .eta <- as.character(.muRefDf$eta[.w])
    if (!(.eta %in% .cls$muCovEtas)) return(NULL)
    if (!clamp && .theta %in% .boundedThetas) {
      warning(
        "mu-referenced theta '", .theta, "' has boundaries and cannot ",
        "benefit from the mu-referenced speed gains (this method's ",
        "regression cannot respect a boundary on its population theta); ",
        "estimated as an ordinary (bounded) parameter by the outer ",
        "optimizer instead",
        call. = FALSE
      )
      return(NULL)
    }
    if (.theta %in% .fixedThetas) {
      .minfo(paste0("mu-referenced theta '", .theta,
                    "' is fixed; kept out of the mu-referenced regression"))
      return(NULL)
    }
    .wc <- which(as.character(.muRefCovDf$theta) == .theta)
    .covParams <- as.character(.muRefCovDf$covariateParameter[.wc])
    .covNames <- as.character(.muRefCovDf$covariate[.wc])
    .boundedCov <- !clamp & (.covParams %in% .boundedThetas)
    if (any(.boundedCov)) {
      warning(
        "mu-referenced covariate coefficient(s) ",
        paste0("'", .covParams[.boundedCov], "'", collapse = ", "),
        " have boundaries and cannot benefit from the mu-referenced ",
        "speed gains (this method's regression cannot respect a ",
        "boundary); treated as if time-varying and excluded from the ",
        "mu-referencing -- estimated as ordinary (bounded) parameters ",
        "by the outer optimizer instead, while '", .theta,
        "' and any other covariate(s) on it still benefit",
        call. = FALSE
      )
    }
    list(
      theta = .theta,
      eta = .eta,
      thetaLower = unname(.lowerOf[.theta]),
      thetaUpper = unname(.upperOf[.theta]),
      covariates = data.frame(
        covariate = .covNames,
        covariateParameter = .covParams,
        lower = unname(.lowerOf[.covParams]),
        upper = unname(.upperOf[.covParams]),
        bounded = .boundedCov,
        stringsAsFactors = FALSE
      )
    )
  }) -> .lst
  .lst <- .lst[!vapply(.lst, is.null, logical(1))]
  if (!plain || length(.cls$muPlainThetas) == 0L) return(.lst)
  .emptyCov <- data.frame(covariate = character(0),
                          covariateParameter = character(0),
                          lower = numeric(0), upper = numeric(0),
                          bounded = logical(0),
                          stringsAsFactors = FALSE)
  .boundedPlain <- character(0)
  .plainLst <- list()
  for (.k in seq_along(.cls$muPlainThetas)) {
    .th <- .cls$muPlainThetas[.k]
    if (.th %in% .fixedThetas) next
    if (!clamp && .th %in% .boundedThetas) {
      .boundedPlain <- c(.boundedPlain, .th)
      next
    }
    .plainLst[[length(.plainLst) + 1L]] <-
      list(theta = .th, eta = .cls$muPlainEtas[.k],
           thetaLower = unname(.lowerOf[.th]),
           thetaUpper = unname(.upperOf[.th]),
           covariates = .emptyCov)
  }
  if (length(.boundedPlain) > 0L) {
    .minfo(paste0("mu-referenced theta(s) ",
                  paste0("'", .boundedPlain, "'", collapse = ", "),
                  " have boundaries; estimated by the outer optimizer ",
                  "instead of the mu-referenced regression"))
  }
  if (length(.plainLst) == 0L) return(.lst)
  # Guard: profiling every theta out with nothing left free (no non-grouped
  # estimated theta, no estimated omega) would give the outer optimizer an
  # empty problem; drop the plain groups (covariate groups kept) instead.
  .grouped <- unlist(lapply(c(.lst, .plainLst), function(g) {
    c(g$theta, g$covariates$covariateParameter[!g$covariates$bounded])
  }))
  .freeTh <- .thetaRows$name[!.thetaRows$fix & !(.thetaRows$name %in% .grouped)]
  .omRows <- .iniDf[is.na(.iniDf$ntheta), , drop = FALSE]
  .freeOm <- any(!.omRows$fix)
  if (length(.freeTh) == 0L && !.freeOm) {
    .minfo("profiling all mu-referenced thetas would leave the outer optimizer with no parameters; keeping them as ordinary parameters")
    return(.lst)
  }
  c(.lst, .plainLst)
}

#' Flatten `.muRefGroups()` into the 0-based index arrays the C++
#' mu-referenced-FOCEI-family regression (`updateMuGroups()`,
#' `src/inner.cpp`) expects
#'
#' UI-derived only; `.muRefCppCovData()` builds the covariate matrix
#' separately once the dataset is available.
#'
#' With `clamp=FALSE` (see `.muRefGroups`), covariates carved out for their
#' boundaries (`bounded=TRUE`) are left out of the flattened arrays
#' entirely, so they stay ordinary outer-optimized parameters (not skipped,
#' not regressed).
#'
#' @param ui rxode2 ui object (post mu2-hook, if applicable)
#' @param plain also include plain mu-ref pairs (see `.muRefGroups`)
#' @param clamp accept bounded parameters for the clamped regression (see
#'   `.muRefGroups`)
#' @return list with `muGroupTheta`, `muGroupEta`, `muGroupCovStart`,
#'   `muGroupCovCount`, `muGroupCovTheta`, `muGroupCovUserFixed` (integer
#'   vectors), `muGroupThetaLower`/`muGroupThetaUpper` (numeric, per group),
#'   `muGroupCovLower`/`muGroupCovUpper` (numeric, parallel to
#'   `muGroupCovTheta`) and `muGroupCovNames` (character, parallel to
#'   `muGroupCovTheta`, used by `.muRefCppCovData()` to pull the right
#'   dataset columns)
#' @author Matthew L. Fidler
#' @noRd
.muRefCppGroupSetup <- function(ui, plain = FALSE,
                                clamp = !identical(rxode2::rxGetControl(ui, "muModel", "none"), "none")) {
  .groups <- .muRefGroups(ui, plain = plain, clamp = clamp)
  if (length(.groups) == 0L) {
    return(list(
      muGroupTheta = integer(0), muGroupEta = integer(0),
      muGroupCovStart = integer(0), muGroupCovCount = integer(0),
      muGroupCovTheta = integer(0), muGroupCovUserFixed = integer(0),
      muGroupThetaLower = numeric(0), muGroupThetaUpper = numeric(0),
      muGroupCovLower = numeric(0), muGroupCovUpper = numeric(0),
      muGroupCovNames = character(0)
    ))
  }
  .iniDf <- ui$iniDf
  .thetaNames <- .iniDf$name[!is.na(.iniDf$ntheta)]
  .thetaIdxOf <- function(nm) match(nm, .thetaNames) - 1L

  .w <- which(!is.na(.iniDf$ntheta))
  .i2 <- .iniDf[-.w, ]
  .i2 <- .i2[.i2$neta1 == .i2$neta2, ]
  .i2 <- .i2[order(.i2$neta1), ]
  .etaNames <- .i2$name
  .etaIdxOf <- function(nm) match(nm, .etaNames) - 1L

  .userFixed <- setNames(.iniDf$fix[!is.na(.iniDf$ntheta)], .thetaNames)

  .muGroupTheta <- integer(0)
  .muGroupEta <- integer(0)
  .muGroupCovStart <- integer(0)
  .muGroupCovCount <- integer(0)
  .muGroupCovTheta <- integer(0)
  .muGroupCovUserFixed <- integer(0)
  .muGroupThetaLower <- numeric(0)
  .muGroupThetaUpper <- numeric(0)
  .muGroupCovLower <- numeric(0)
  .muGroupCovUpper <- numeric(0)
  .muGroupCovNames <- character(0)

  for (g in .groups) {
    .cov <- g$covariates[!g$covariates$bounded, , drop = FALSE]
    .muGroupTheta <- c(.muGroupTheta, .thetaIdxOf(g$theta))
    .muGroupEta <- c(.muGroupEta, .etaIdxOf(g$eta))
    .muGroupThetaLower <- c(.muGroupThetaLower, as.numeric(g$thetaLower))
    .muGroupThetaUpper <- c(.muGroupThetaUpper, as.numeric(g$thetaUpper))
    .muGroupCovStart <- c(.muGroupCovStart, length(.muGroupCovTheta))
    .muGroupCovCount <- c(.muGroupCovCount, nrow(.cov))
    for (.k in seq_len(nrow(.cov))) {
      .cp <- .cov$covariateParameter[.k]
      .muGroupCovTheta <- c(.muGroupCovTheta, .thetaIdxOf(.cp))
      .muGroupCovUserFixed <- c(.muGroupCovUserFixed,
                                 as.integer(isTRUE(.userFixed[[.cp]])))
      .muGroupCovLower <- c(.muGroupCovLower, as.numeric(.cov$lower[.k]))
      .muGroupCovUpper <- c(.muGroupCovUpper, as.numeric(.cov$upper[.k]))
      .muGroupCovNames <- c(.muGroupCovNames, .cov$covariate[.k])
    }
  }

  list(
    muGroupTheta = as.integer(.muGroupTheta),
    muGroupEta = as.integer(.muGroupEta),
    muGroupCovStart = as.integer(.muGroupCovStart),
    muGroupCovCount = as.integer(.muGroupCovCount),
    muGroupCovTheta = as.integer(.muGroupCovTheta),
    muGroupCovUserFixed = as.integer(.muGroupCovUserFixed),
    muGroupThetaLower = .muGroupThetaLower,
    muGroupThetaUpper = .muGroupThetaUpper,
    muGroupCovLower = .muGroupCovLower,
    muGroupCovUpper = .muGroupCovUpper,
    muGroupCovNames = .muGroupCovNames
  )
}

#' Build the baseline (first-observation) covariate-value matrix for the
#' mu-referenced-FOCEI-family C++ regression
#'
#' One row per subject (sorted by ID, matching `etaMat`'s ordering), one
#' column per flattened covariate entry from `.muRefCppGroupSetup()` in the
#' same order (duplicate covariate names across groups get their own
#' column, since `updateMuGroups()` slices by flattened index, not name).
#'
#' @param covNames character vector, `muGroupCovNames` from
#'   `.muRefCppGroupSetup()`
#' @param dataSav processed dataset with an `ID` column and the needed
#'   covariate columns (`env$dataSav` inside `.foceiFamilyReturn()`)
#' @return numeric matrix, `nsub x length(covNames)`
#' @author Matthew L. Fidler
#' @noRd
.muRefCppCovData <- function(covNames, dataSav) {
  if (length(covNames) == 0L) return(matrix(numeric(0), nrow = 0, ncol = 0))
  .byId <- dataSav[!duplicated(dataSav$ID), , drop = FALSE]
  .byId <- .byId[order(.byId$ID), , drop = FALSE]
  # a name can be a covariate expression (e.g. log(WT/70)) rather than a
  # bare column; evaluate those against the baseline rows
  .cols <- lapply(covNames, function(.cn) {
    if (.cn %in% names(.byId)) return(as.numeric(.byId[[.cn]]))
    as.numeric(eval(parse(text = .cn), envir = .byId))
  })
  matrix(unlist(.cols, use.names = FALSE),
         nrow = nrow(.byId), ncol = length(covNames))
}
