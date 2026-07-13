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
#' coefficient stays out of the design matrix as a live offset. Bounded
#' population thetas and covariate coefficients ARE grouped: the regression
#' update is clamped to the bounds (box-constrained least squares,
#' `updateMuGroups()`), and each group carries `thetaLower`/`thetaUpper` plus
#' per-covariate `lower`/`upper` columns for the clamp. Only affects this
#' family (`mfocei`/`ifocei`/`mfoce`/`ifoce`/`magq`/`iagq`/
#' `mlaplace`/`ilaplace`).
#'
#' With `plain=TRUE`, plain (covariate-free) mu-ref pairs
#' (`muPlainThetas`, see `.muRefClassify`) are appended as intercept-only
#' groups (empty `covariates`); fixed plain thetas are excluded quietly. If
#' profiling the plain groups would leave the outer optimizer with no free
#' parameter at all, the plain groups are dropped (covariate groups kept).
#'
#' @param ui rxode2 ui object (post mu2-hook, if applicable)
#' @param plain also form intercept-only groups for plain mu-ref pairs
#' @return list of `list(theta=, eta=, thetaLower=, thetaUpper=,
#'   covariates=data.frame(covariate=, covariateParameter=, lower=,
#'   upper=))`, one element per group whose population theta is not
#'   user-fixed
#' @author Matthew L. Fidler
#' @noRd
.muRefGroups <- function(ui, plain = FALSE) {
  .cls <- .muRefClassify(ui)
  .muRefDf <- ui$muRefDataFrame
  .muRefCovDf <- ui$muRefCovariateDataFrame
  .iniDf <- ui$iniDf
  .thetaRows <- .iniDf[!is.na(.iniDf$ntheta), , drop = FALSE]
  .fixedThetas <- .thetaRows$name[which(.thetaRows$fix)]
  .lowerOf <- setNames(.thetaRows$lower, .thetaRows$name)
  .upperOf <- setNames(.thetaRows$upper, .thetaRows$name)
  lapply(.cls$muCovThetas, function(.theta) {
    .w <- which(as.character(.muRefDf$theta) == .theta)
    if (length(.w) != 1L) return(NULL)
    .eta <- as.character(.muRefDf$eta[.w])
    if (!(.eta %in% .cls$muCovEtas)) return(NULL)
    if (.theta %in% .fixedThetas) {
      .minfo(paste0("mu-referenced theta '", .theta,
                    "' is fixed; kept out of the mu-referenced regression"))
      return(NULL)
    }
    .wc <- which(as.character(.muRefCovDf$theta) == .theta)
    .covParams <- as.character(.muRefCovDf$covariateParameter[.wc])
    .covNames <- as.character(.muRefCovDf$covariate[.wc])
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
        stringsAsFactors = FALSE
      )
    )
  }) -> .lst
  .lst <- .lst[!vapply(.lst, is.null, logical(1))]
  if (!plain || length(.cls$muPlainThetas) == 0L) return(.lst)
  .emptyCov <- data.frame(covariate = character(0),
                          covariateParameter = character(0),
                          lower = numeric(0), upper = numeric(0),
                          stringsAsFactors = FALSE)
  .plainLst <- list()
  for (.k in seq_along(.cls$muPlainThetas)) {
    .th <- .cls$muPlainThetas[.k]
    if (.th %in% .fixedThetas) next
    .plainLst[[length(.plainLst) + 1L]] <-
      list(theta = .th, eta = .cls$muPlainEtas[.k],
           thetaLower = unname(.lowerOf[.th]),
           thetaUpper = unname(.upperOf[.th]),
           covariates = .emptyCov)
  }
  if (length(.plainLst) == 0L) return(.lst)
  # Guard: profiling every theta out with nothing left free (no non-grouped
  # estimated theta, no estimated omega) would give the outer optimizer an
  # empty problem; drop the plain groups (covariate groups kept) instead.
  .grouped <- unlist(lapply(c(.lst, .plainLst), function(g) {
    c(g$theta, g$covariates$covariateParameter)
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
#' @param ui rxode2 ui object (post mu2-hook, if applicable)
#' @param plain also include plain mu-ref pairs (see `.muRefGroups`)
#' @return list with `muGroupTheta`, `muGroupEta`, `muGroupCovStart`,
#'   `muGroupCovCount`, `muGroupCovTheta`, `muGroupCovUserFixed` (integer
#'   vectors), `muGroupThetaLower`/`muGroupThetaUpper` (numeric, per group),
#'   `muGroupCovLower`/`muGroupCovUpper` (numeric, parallel to
#'   `muGroupCovTheta`) and `muGroupCovNames` (character, parallel to
#'   `muGroupCovTheta`, used by `.muRefCppCovData()` to pull the right
#'   dataset columns)
#' @author Matthew L. Fidler
#' @noRd
.muRefCppGroupSetup <- function(ui, plain = FALSE) {
  .groups <- .muRefGroups(ui, plain = plain)
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
    .muGroupTheta <- c(.muGroupTheta, .thetaIdxOf(g$theta))
    .muGroupEta <- c(.muGroupEta, .etaIdxOf(g$eta))
    .muGroupThetaLower <- c(.muGroupThetaLower, as.numeric(g$thetaLower))
    .muGroupThetaUpper <- c(.muGroupThetaUpper, as.numeric(g$thetaUpper))
    .muGroupCovStart <- c(.muGroupCovStart, length(.muGroupCovTheta))
    .muGroupCovCount <- c(.muGroupCovCount, nrow(g$covariates))
    for (.k in seq_len(nrow(g$covariates))) {
      .cp <- g$covariates$covariateParameter[.k]
      .muGroupCovTheta <- c(.muGroupCovTheta, .thetaIdxOf(.cp))
      .muGroupCovUserFixed <- c(.muGroupCovUserFixed,
                                 as.integer(isTRUE(.userFixed[[.cp]])))
      .muGroupCovLower <- c(.muGroupCovLower, as.numeric(g$covariates$lower[.k]))
      .muGroupCovUpper <- c(.muGroupCovUpper, as.numeric(g$covariates$upper[.k]))
      .muGroupCovNames <- c(.muGroupCovNames, g$covariates$covariate[.k])
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
