#' Free Poped memory (if any is allocated)
#'
#' This should not be called directly but is used in babelmixr2's
#' poped interface
#'
#' @return nothing, called for side effects
#'
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.popedFree <- function() {
  invisible(.Call(`_nlmixr2est_popedFree`))
}

#' Setup the PopED environment
#'
#' @param e environment with setup information for popEd
#' @return nothing, called for side effects
#' @export
#' @author Matthew L. Fidler
.popedSetup <- function(e) {
  invisible(.Call(`_nlmixr2est_popedSetup`, e))
}
#' After the solver has been setup in memory (with the appropriate
#' data), this solves for the appropriate sample times
#'
#' This really should not be called directly (if not setup correctly
#' can crash R)
#'
#' @param theta parameters (includes covariates)
#' @param id this is the design identifier
#' @param totn This is the total number of design points tested
#' @return a data frame with $f and $w corresponding to the function
#'   value and standard deviation at the sampling point
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.popedSolveIdN <- function(theta, id, totn) {
  .Call(`_nlmixr2est_popedSolveIdN`, theta, id, totn)
}
