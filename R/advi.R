# advi.R -- orchestration for est="advi" (automatic differentiation variational
# inference, Kucukelbir et al. 2017).  Sets up the FOCEi inner problem (reused
# for the per-subject log-joint and eta-gradient) plus, when non-mu structural
# thetas are present, the impmap theta-sensitivity model (reused for the outer
# population gradient), then drives the ADVI optimization loop in C++.

#' Fit an ADVI model: set up the inner/outer problems and run the C++ loop.
#' @param env estimation environment (holds ui, data, adviControl)
#' @noRd
.adviFitModel <- function(env) {
  stop("est=\"advi\" is not yet implemented", call. = FALSE)
}
