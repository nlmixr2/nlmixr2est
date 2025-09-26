#' Get the mixture probabilities from the estimated log-scale parameters
#'
#'
#' @param val numeric vector of log-scale parameters
#' @return A numeric vector of the mixture probabilities
#' @noRd
#' @author Matthew L. Fidler
.getMixFromLog <- function(val) {
  v <- do.call("rxode2::mexpit", list(val))
  c(v, 1-sum(v))
}
