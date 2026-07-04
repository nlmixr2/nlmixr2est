#' Classify theta/eta pairs for the mu-referenced FOCEI family
#'
#' Only thetas/etas with a mu-ref covariate relationship
#' (`ui$muRefCovariateDataFrame`) get the in-C++ regression machinery
#' (`updateMuGroups()`, `src/inner.cpp`); plain pairs use standard FOCEI
#' inner/outer optimization (unlike SAEM, which mu-refs the plain case too).
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
#'  \item{\code{standardThetas}}{population thetas left untouched by this
#'    family; excludes sigma/residual-error parameters}
#'  \item{\code{standardEtas}}{etas left untouched by this family}
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

  list(
    muCovThetas = .muCovThetas,
    muCovEtas = .muCovEtas,
    muCovCovariateParams = .muCovCovariateParams,
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
#' A finite bound on the group's population theta excludes the whole group
#' (warns) since OLS/IRLS can't respect a box constraint on its intercept.
#' A bound on one covariate-coefficient theta excludes only that covariate
#' (marked `bounded=TRUE`, treated like a time-varying covariate), while
#' the population theta and other unbounded covariates in the group still
#' get the mu-ref speed-up. Only affects this family (`mufocei`/`irlsfocei`/
#' `mufoce`/`irlsfoce`/`muagq`/`irlsagq`/`mulaplace`/`irlslaplace`).
#'
#' @param ui rxode2 ui object (post mu2-hook, if applicable)
#' @return list of `list(theta=, eta=, covariates=data.frame(covariate=,
#'   covariateParameter=, bounded=))`, one element per mu-ref-covariate-
#'   with-eta group whose population theta is unbounded
#' @author Matthew L. Fidler
#' @noRd
.muRefGroups <- function(ui) {
  .cls <- .muRefClassify(ui)
  .muRefDf <- ui$muRefDataFrame
  .muRefCovDf <- ui$muRefCovariateDataFrame
  .iniDf <- ui$iniDf
  .thetaRows <- .iniDf[!is.na(.iniDf$ntheta), , drop = FALSE]
  .boundedThetas <- .thetaRows$name[.thetaRows$lower > -Inf | .thetaRows$upper < Inf]
  lapply(.cls$muCovThetas, function(.theta) {
    .w <- which(as.character(.muRefDf$theta) == .theta)
    if (length(.w) != 1L) return(NULL)
    .eta <- as.character(.muRefDf$eta[.w])
    if (!(.eta %in% .cls$muCovEtas)) return(NULL)
    if (.theta %in% .boundedThetas) {
      warning(
        "mu-referenced theta '", .theta, "' has boundaries and cannot ",
        "benefit from the mu-referenced FOCEI-family speed gains (the ",
        "closed-form/IRLS regression cannot respect a boundary on its ",
        "population theta); estimated as an ordinary (bounded) ",
        "parameter by the outer optimizer instead",
        call. = FALSE
      )
      return(NULL)
    }
    .wc <- which(as.character(.muRefCovDf$theta) == .theta)
    .covParams <- as.character(.muRefCovDf$covariateParameter[.wc])
    .covNames <- as.character(.muRefCovDf$covariate[.wc])
    .boundedCov <- .covParams %in% .boundedThetas
    if (any(.boundedCov)) {
      warning(
        "mu-referenced covariate coefficient(s) ",
        paste0("'", .covParams[.boundedCov], "'", collapse = ", "),
        " have boundaries and cannot benefit from the mu-referenced ",
        "FOCEI-family speed gains (the closed-form/IRLS regression ",
        "cannot respect a boundary); treated as if time-varying and ",
        "excluded from the mu-referencing -- estimated as ordinary ",
        "(bounded) parameters by the outer optimizer instead, while '",
        .theta, "' and any other covariate(s) on it still benefit",
        call. = FALSE
      )
    }
    list(
      theta = .theta,
      eta = .eta,
      covariates = data.frame(
        covariate = .covNames,
        covariateParameter = .covParams,
        bounded = .boundedCov,
        stringsAsFactors = FALSE
      )
    )
  }) -> .lst
  .lst[!vapply(.lst, is.null, logical(1))]
}

#' Flatten `.muRefGroups()` into the 0-based index arrays the C++
#' mu-referenced-FOCEI-family regression (`updateMuGroups()`,
#' `src/inner.cpp`) expects
#'
#' UI-derived only; `.muRefCppCovData()` builds the covariate matrix
#' separately once the dataset is available.
#'
#' @param ui rxode2 ui object (post mu2-hook, if applicable)
#' @return list with `muGroupTheta`, `muGroupEta`, `muGroupCovStart`,
#'   `muGroupCovCount`, `muGroupCovTheta`, `muGroupCovUserFixed`,
#'   `muGroupCovBounded` (all integer vectors) and `muGroupCovNames`
#'   (character vector, parallel to `muGroupCovTheta`, used by
#'   `.muRefCppCovData()` to pull the right dataset columns)
#' @author Matthew L. Fidler
#' @noRd
.muRefCppGroupSetup <- function(ui) {
  .groups <- .muRefGroups(ui)
  if (length(.groups) == 0L) {
    return(list(
      muGroupTheta = integer(0), muGroupEta = integer(0),
      muGroupCovStart = integer(0), muGroupCovCount = integer(0),
      muGroupCovTheta = integer(0), muGroupCovUserFixed = integer(0),
      muGroupCovBounded = integer(0), muGroupCovNames = character(0)
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
  .muGroupCovBounded <- integer(0)
  .muGroupCovNames <- character(0)

  for (g in .groups) {
    .muGroupTheta <- c(.muGroupTheta, .thetaIdxOf(g$theta))
    .muGroupEta <- c(.muGroupEta, .etaIdxOf(g$eta))
    .muGroupCovStart <- c(.muGroupCovStart, length(.muGroupCovTheta))
    .muGroupCovCount <- c(.muGroupCovCount, nrow(g$covariates))
    for (.k in seq_len(nrow(g$covariates))) {
      .cp <- g$covariates$covariateParameter[.k]
      .muGroupCovTheta <- c(.muGroupCovTheta, .thetaIdxOf(.cp))
      .muGroupCovUserFixed <- c(.muGroupCovUserFixed,
                                 as.integer(isTRUE(.userFixed[[.cp]])))
      .muGroupCovBounded <- c(.muGroupCovBounded,
                               as.integer(isTRUE(g$covariates$bounded[.k])))
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
    muGroupCovBounded = as.integer(.muGroupCovBounded),
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
  as.matrix(.byId[, covNames, drop = FALSE])
}
