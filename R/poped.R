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
#' This should not typically be called directly
#'
#' @param e environment with setup information for popEd
#' @param full setup the full model
#' @return nothing, called for side effects
#' @export
#' @keywords internal
#' @author Matthew L. Fidler
.popedSetup <- function(e, full=FALSE) {
  invisible(.Call(`_nlmixr2est_popedSetup`, e, full))
}
#' Solve poped problem for appropriate times (may already be setup)
#'
#' This really should not be called directly (if not setup correctly
#' can crash R)
#'
#' @param theta parameters (includes covariates)
#' @param xt original unsorted time (to match the f/w against)
#' @param id this is the design identifier
#' @param totn This is the total number of design points tested
#' @return a data frame with $f and $w corresponding to the function
#'   value and standard deviation at the sampling point
#' @export
#' @author Matthew L. Fidler
#' @keywords internal
.popedSolveIdN <- function(theta, xt, id, totn) {
  .Call(`_nlmixr2est_popedSolveIdN`, theta, xt, id, totn)
}
#' @rdname .popedSolveIdN
#' @export
.popedSolveIdN2 <- function(theta, xt, id, totn) {
  .Call(`_nlmixr2est_popedSolveIdN2`, theta, xt, id, totn)
}
