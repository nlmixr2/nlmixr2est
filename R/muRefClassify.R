#' Classify theta/eta pairs for the mu-referenced FOCEI family
#'
#' The mu-referenced FOCEI family (`mufocei`/`irlsfocei`/etc.) only applies
#' its restart-loop/linear-model machinery to thetas and etas that
#' participate in a **mu-ref covariate** relationship
#' (`ui$muRefCovariateDataFrame`). A theta+eta pair with no covariate ("eta
#' by itself") and a theta with no covariate and no eta ("theta by itself")
#' are left completely alone, handled by standard FOCEI inner/outer
#' optimization exactly as today -- this is a deliberate difference from
#' SAEM's mu-referencing, which mu-refs even the plain theta+eta case.
#'
#' Complex covariate expressions (mu2/mu3/mu4-eligible, see `R/mu2.R`) are
#' expected to have already been algebraically rewritten into simple linear
#' `muRefCovariateDataFrame` terms by the existing `.uiApplyMu2hook()`
#' pre-processing hook before this function is called -- that hook is
#' triggered generically by `.isMuMethod()` checking the `"mu"` S3 attribute
#' on `nlmixr2Est.<method>`, exactly as SAEM/nlme already do (see
#' `R/mu2.R`). This function does not need to know anything about mu2/mu3/mu4
#' itself; it only reads the already-simplified `ui$muRefCovariateDataFrame`.
#'
#' @param ui rxode2 ui object (post mu2-hook, if applicable)
#' @return A list with:
#' \itemize{
#'  \item{\code{muCovThetas}}{character vector of population thetas that
#'    have at least one mu-ref covariate relationship (with or without an
#'    associated eta)}
#'  \item{\code{muCovEtas}}{character vector of etas whose associated theta
#'    is in \code{muCovThetas}}
#'  \item{\code{muCovCovariateParams}}{character vector of the
#'    covariate-coefficient thetas themselves (the "slopes")}
#'  \item{\code{standardThetas}}{character vector of population thetas left
#'    untouched by this family (no covariate involvement); excludes sigma/
#'    residual-error parameters}
#'  \item{\code{standardEtas}}{character vector of etas left untouched by
#'    this family (no covariate involvement on their associated theta, or
#'    not mu-referenced to any theta at all)}
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

#' Build per-group theta/eta/covariate structures for the restart-loop engine
#'
#' Each group is one mu-ref-covariate population theta that also has an
#' associated eta (e.g. `cl <- exp(tcl + eta.cl + allo.cl*logWT)` is one
#' group: `theta="tcl"`, `eta="eta.cl"`, one covariate row `covariate=
#' "logWT", covariateParameter="allo.cl"`; a theta can have more than one
#' covariate).
#'
#' Mu-ref-covariate thetas with **no** associated eta (a pure fixed-effect
#' covariate relationship, e.g. bioavailability with no random effect) are
#' intentionally **not** included here -- there is no per-subject residual
#' to regress against, so the restart-loop's linear-model step has nothing
#' to update them with. They are left out of the "fixed for the outer
#' optimizer" set entirely and behave as ordinary, outer-optimized thetas.
#' This is a deliberate scoping decision for the initial engine, not an
#' oversight.
#'
#' A group's *population* theta having a finite bound
#' (`ini(...~c(lower, est, upper))`) excludes the **whole group** (with a
#' `warning()`, captured into the fit's `runInfo` the same way every
#' other pre-fit warning is) -- the population theta is the regression's
#' intercept, and an unconstrained OLS/IRLS solve has no way to respect a
#' box constraint on it.
#'
#' A bound on one of the group's *covariate-coefficient* thetas (the
#' "slopes") is handled more narrowly: only that one covariate is
#' excluded from the regression (marked `bounded=TRUE` in the returned
#' `covariates` data frame) -- it is treated exactly like a time-varying
#' covariate would be (see `.muRefClassify()`'s docs: a plain covariate
#' relationship with no eta is never mu-ref-eligible to begin with), i.e.
#' left to ordinary, bounded outer-optimizer handling while its *current*
#' (not fixed -- gradient-updated every real outer iteration) contribution
#' is still subtracted from the regression's target as a live offset,
#' exactly like `.muRefLin()`'s `fixedCoef` contract but recomputed each
#' iteration instead of held constant (`src/inner.cpp`'s
#' `updateMuGroups()`). The population theta and any of the group's
#' *other*, unbounded covariates still get the full mu-ref speed-up.
#'
#' Both cases warn. Both only ever apply to this function (and therefore
#' only to the `mufocei`/`irlsfocei`/`mufoce`/`irlsfoce`/`muagq`/
#' `irlsagq`/`mulaplace`/`irlslaplace` family) -- `.muRefGroups()` has
#' exactly one caller (`.muRefCppGroupSetup()`, itself only reached from
#' `rxUiGet.foceiOptEnv` when `muModel != "none"`), so SAEM's own,
#' unrelated mu-referencing mechanism and every other estimator are
#' untouched.
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
#' Purely UI-derived (no dataset needed) -- the covariate *values* matrix
#' is built separately by `.muRefCppCovData()` once the dataset is
#' available (see `R/focei.R`'s `.foceiFamilyReturn()`).
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
#' One row per subject, sorted by ID (matching the subject ordering
#' `foceiSetup_()` already uses for `etaMat`), one column per flattened
#' covariate entry from `.muRefCppGroupSetup()` (in the same order --
#' duplicate covariate names across groups get their own column, which is
#' fine, `src/inner.cpp`'s `updateMuGroups()` slices by flattened index,
#' not by name).
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
