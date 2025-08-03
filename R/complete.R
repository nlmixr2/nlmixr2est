
.nmObjGetEnvInfo <- list(
  ui="rxode2 user interface",
  conditionNumberCor="Condition Number (Correlation)",
  conditionNumberCov="Condition Number (Covariance)",
  cov="Covariance of fixed effects",
  covMethod="Covariance Method for fixed effects",
  etaObf="ETAs and their individual objective function contribution (if applicable)",
  objDf="Objective Function DF",
  omega="Omega Matrix",
  origData="Original Data",
  phiC="covariance matrix of each individual's eta (if present)",
  parFixed="Formatted Parameter Values for Fixed effects",
  parFixedDf="Parameter Values for Fixed Effects (data frame)",
  parHistData="Parameter History (including gradients)",
  scaleInfo="Scaling Information",
  shrink="Shrinkage data frame",
  table="Table Control Value",
  fixef="Fixed effects",
  time="Timing data frame"
)

.nmObjGetSupportedDollars <- function() {
  .v <- as.character(utils::methods("nmObjGet"))
  .v <- .v[.v != "nmObjGet.default"]
  .cls <- vapply(.v, function(methodStr){
    substr(methodStr,10,nchar(methodStr))
  }, character(1), USE.NAMES=FALSE)
  .v <- vapply(.cls, function(cls){
    .desc <- attr(utils::getS3method("nmObjGet", cls), "desc")
    if (is.null(.desc)) .desc <- ""
    .desc
  }, character(1), USE.NAMES=TRUE)
  # Take out any "hidden methods"
  .w <- which(.v != "")
  .v <- c(.v[.w], .nmObjGetEnvInfo)
  .v
}

##' @export
.DollarNames.nlmixr2FitCore <- function(x, pattern) {
  .env <- x$env
  .cmp <- c(
    names(x),
    names(.nmObjGetSupportedDollars()))
  .cmp <- c(.cmp, "env")
  grep(pattern, .cmp, value = TRUE)
}

.nmObjGetDataSupportedDollars <- function() {
  .v <- as.character(utils::methods("nmObjGetData"))
  .v <- .v[.v != "nmObjGetData.default"]
  .cls <- vapply(.v, function(methodStr){
    substr(methodStr, 14, nchar(methodStr))
  }, character(1), USE.NAMES=FALSE)
  .v <- vapply(.cls, function(cls){
    .desc <- attr(utils::getS3method("nmObjGetData", cls), "desc")
    if (is.null(.desc)) .desc <- ""
    .desc
  }, character(1), USE.NAMES=TRUE)
  # Take out any "hidden methods"
  .w <- which(.v != "")
  .v <- c(.v[.w], .nmObjGetEnvInfo)
  c(.nmObjGetSupportedDollars(), .v)
}

#' @export
.DollarNames.nlmixr2FitData <- function(x, pattern) {
  ##FIXME
  .env <- x$env
  .cmp <- c(
    names(x),
    names(.nmObjGetDataSupportedDollars()))
  .cmp <- c(.cmp, "env")
  grep(pattern, .cmp, value = TRUE)
}
