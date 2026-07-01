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
#' @param ui rxode2 ui object (post mu2-hook, if applicable)
#' @return list of `list(theta=, eta=, covariates=data.frame(covariate=,
#'   covariateParameter=))`, one element per mu-ref-covariate-with-eta group
#' @author Matthew L. Fidler
#' @noRd
.muRefGroups <- function(ui) {
  .cls <- .muRefClassify(ui)
  .muRefDf <- ui$muRefDataFrame
  .muRefCovDf <- ui$muRefCovariateDataFrame
  lapply(.cls$muCovThetas, function(.theta) {
    .w <- which(as.character(.muRefDf$theta) == .theta)
    if (length(.w) != 1L) return(NULL)
    .eta <- as.character(.muRefDf$eta[.w])
    if (!(.eta %in% .cls$muCovEtas)) return(NULL)
    .wc <- which(as.character(.muRefCovDf$theta) == .theta)
    list(
      theta = .theta,
      eta = .eta,
      covariates = data.frame(
        covariate = as.character(.muRefCovDf$covariate[.wc]),
        covariateParameter = as.character(.muRefCovDf$covariateParameter[.wc]),
        stringsAsFactors = FALSE
      )
    )
  }) -> .lst
  .lst[!vapply(.lst, is.null, logical(1))]
}
