#' Get the mixture probabilities from the estimated log-scale parameters
#'
#' @param val numeric vector of the full parameter set in focei
#'
#' @param idx integer vector of the indices of the mixture log-scale
#'   parameters
#'
#' @return A numeric vector of the mixture probabilities
#'
#' @noRd
#'
#' @author Matthew L. Fidler
.getMixFromLog <- function(val, idx) {
  v <- do.call("rxode2::mexpit", list(val[idx]))
  c(v, 1-sum(v))
}
#' Get the mixture gradients of the estimated log-scale parameters
#'
#'
#' @param val numeric vector of the full parameter set in focei
#'
#' @param idx vector of the indices of the mixture log-scale
#'  parameters
#'
#' @return A numeric vector of the mixture probabilities
#'
#' @noRd
#'
#' @author Matthew L. Fidler
.getMixJacFromLog <- function(val, idx) {
  do.call("rxode2::dmexpit", list(val[idx]))
}

#' @export
rxUiGet.thetaIniMix <- function(x, ...) {
  .ui <- x[[1]]
  .theta <- .ui$theta
  if (length(.ui$mixProbs) > 0) {
    .theta[.ui$mixProbs] <- rxode2::mlogit(.theta[.ui$mixProbs])
  }
  .theta
}
attr(rxUiGet.thetaIniMix, "rstudio") <- stats::setNames(1, "a")

#' @export
rxUiGet.thetaMixIndex <- function(x, ...) {
  .ui <- x[[1]]
  .theta <- .ui$theta
  if (length(.ui$mixProbs) > 0) {
    # Optimized: Use match() instead of which() + %in%
    # match(x, y) returns indices in y where x elements are found
    # This is equivalent to which(y %in% x) but faster
    .m <- match(.ui$mixProbs, names(.ui$theta))
    # Return only non-NA indices (where matches were found)
    .m[!is.na(.m)]
  } else {
    integer(0)
  }
}
attr(rxUiGet.thetaMixIndex, "rstudio") <- 1L
