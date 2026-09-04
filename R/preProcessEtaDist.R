## Declared non-Gaussian random effect (eta) distributions.
##
## `lotri` parses `dist(eta.cl) ~ dgamma(...)` and `rxode2` turns it into
## a model (`rxEtaDistExpand()`): a latent standard normal with a FIXED
## identity omega, plus `phiU()` + the family's inverse CDF, plus
## unconstrained `rxCor.*` thetas carrying the Gaussian copula's
## correlation.  Two things are left for nlmixr2est.
##
## 1. Run that expansion before estimation, which is all the support most
##    methods need: what they are handed afterwards is an ordinary model
##    with a fixed identity omega.
##
## 2. Refuse the methods for which it is NOT ordinary.  A declared
##    distribution that a method quietly ignored would fit a different
##    model than the one written, with nothing to say so -- the same
##    reasoning, and the same attribute-on-the-S3-method mechanism, as the
##    prior gate in R/priors.R.
##
##      attr(nlmixr2Est.myMethod, "etaDist") <- TRUE
##
## The methods that support it are the FOCEi family (whose inner problem
## needs only d(eta)/d(latent), which rxode2 differentiates exactly
## through `phiU()` and the inverse CDF), SAEM (a declared eta has no
## `theta + eta` form, so it lands in the already-exercised `nonMuEtas`
## path and is still Gibbs/Metropolis sampled with the same sample
## covariance update), and simulation.  Refused: `npag`/`npb`, which model
## the random effect distribution nonparametrically, so a declared one
## contradicts them outright; `nlme` and `nls`, which are Gaussian by
## construction; and `vae`/`emvi`/`fbvi`, whose ELBO hardcodes the normal
## family.

#' Does the dispatched estimation method support a declared eta distribution?
#'
#' Reads the `"etaDist"` attribute of the `nlmixr2Est.<method>` S3 method,
#' so a method registered by another package can declare support without
#' editing this file.  The attribute may be `TRUE`/`FALSE` or a
#' `function(control)` returning a logical.
#'
#' @param est estimation routine name
#' @param control control object
#' @return boolean
#' @noRd
#' @author Matthew L. Fidler
.isEtaDistMethod <- function(est, control=NULL) {
  if (!is.character(est) || length(est) != 1L) return(FALSE)
  .v <- as.character(utils::methods("nlmixr2Est"))
  if (!(paste0("nlmixr2Est.", est) %in% .v)) return(FALSE)
  .a <- attr(utils::getS3method("nlmixr2Est", est), "etaDist")
  if (is.null(.a)) return(FALSE)
  if (is.function(.a)) return(isTRUE(.a(control)))
  isTRUE(.a)
}

#' Refuse a declared eta distribution the dispatched method cannot use
#'
#' @param env nlmixr2 estimation environment
#' @return nothing, called for the error
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2AssertEtaDist <- function(env) {
  .ui <- try(get("ui", envir=env), silent=TRUE)
  if (inherits(.ui, "try-error") || is.null(.ui)) return(invisible())
  .d <- .rxUiEtaDists(.ui)
  if (nrow(.d) == 0L) return(invisible())
  .est <- class(env)[1]
  .control <- if (exists("control", envir=env)) get("control", envir=env) else NULL
  if (.isEtaDistMethod(.est, .control)) return(invisible())
  stop("est=\"", .est, "\" cannot use the declared non-normal random effect ",
       "distribution(s) on '", paste(.d$name, collapse="', '"), "'",
       call.=FALSE)
}

#' `rxode2::rxUiEtaDists()` when the installed rxode2 has it
#'
#' Looked up rather than called directly so that an older rxode2 -- which
#' cannot have produced a declaration in the first place -- degrades to
#' "no declarations" instead of erroring.
#'
#' @param ui rxode2 ui
#' @return the declaring random effects, zero rows when there are none
#' @noRd
#' @author Matthew L. Fidler
.rxUiEtaDists <- function(ui) {
  .ns <- asNamespace("rxode2")
  if (!exists("rxUiEtaDists", envir=.ns, inherits=FALSE)) {
    return(data.frame(name=character(0), etaDist=character(0),
                      stringsAsFactors=FALSE))
  }
  get("rxUiEtaDists", envir=.ns)(ui)
}

#' Pre-processing hook: expand declared eta distributions
#'
#' @param ui rxode2 ui object
#' @param est estimation routine name
#' @param data data
#' @param control control object
#' @return list with the expanded `ui`, or NULL when there is nothing to do
#' @noRd
#' @author Matthew L. Fidler
.preProcessEtaDist <- function(ui, est, data, control) {
  if (is.null(ui)) return(NULL)
  .d <- .rxUiEtaDists(ui)
  if (nrow(.d) == 0L) return(NULL)
  ## The gate in nlmixr2Est() has already refused a method that cannot use
  ## these; a hook also runs for `est="rxSolve"` (simulate/vpc/augPred),
  ## which is supported.
  list(ui=rxode2::rxEtaDistExpand(ui))
}

preProcessHooksAdd(".preProcessEtaDist", .preProcessEtaDist)

#' The correlation matrix a fit's `rxCor.*` thetas encode
#'
#' Inverts the row-normalized Cholesky parameterization
#' `rxEtaDistExpand()` writes:
#'
#'   L[i, j] = tanh(y[i, j]) * s[i, j - 1],  L[i, i] = s[i, i - 1]
#'
#' with `s[i, 0] = 1` and `s[i, j] = s[i, j-1]*sqrt(1 - tanh(y[i,j])^2)`,
#' then returns `R = L L'`.
#'
#' @param nms the block's random effect names, in block order
#' @param y named numeric vector of the `rxCor.<i>.<j>` estimates
#' @return the correlation matrix, with `nms` as dimnames
#' @noRd
#' @author Matthew L. Fidler
.etaDistCorFromY <- function(nms, y) {
  .k <- length(nms)
  .L <- diag(.k)
  for (.i in seq_len(.k)) {
    .s <- 1.0
    for (.j in seq_len(.i - 1L)) {
      .t <- tanh(y[[paste0("rxCor.", nms[.i], ".", nms[.j])]])
      .L[.i, .j] <- .t * .s
      .s <- .s * sqrt(1 - .t * .t)
    }
    .L[.i, .i] <- .s
  }
  .R <- .L %*% t(.L)
  dimnames(.R) <- list(nms, nms)
  .R
}

#' Post-estimation hook: report the declared distributions on the fit
#'
#' The expansion is not something the user wrote, so the fit says what it
#' did: `$etaDist` is the declarations as written, and `$etaDistCor` is
#' the copula correlation matrix of each declared block, rebuilt from the
#' `rxCor.*` estimates.
#'
#' The `rxCor.*` rows themselves already read as correlations in the
#' back-transformed column -- `rxEtaDistExpand()` gives them
#' `backTransform("tanh")`, and `tanh()` of one of them is the partial
#' correlation between its two random effects given the ones before them,
#' which for a 2x2 block is simply the correlation.
#'
#' @param ret the finalized fit
#' @return the fit, with `$etaDist` and `$etaDistCor` when the model
#'   declared a distribution
#' @noRd
#' @author Matthew L. Fidler
.postFinalEtaDist <- function(ret) {
  if (!inherits(ret, "nlmixr2FitCore") || is.null(ret$env)) return(NULL)
  .ui <- try(ret$ui, silent=TRUE)
  if (inherits(.ui, "try-error") || is.null(.ui)) return(NULL)
  .info <- try(get("etaDistInfo", envir=.ui), silent=TRUE)
  if (inherits(.info, "try-error") || is.null(.info)) return(NULL)
  .fix <- ret$env$fixef
  if (is.null(.fix)) return(NULL)
  .cor <- lapply(.info$blocks, function(.nms) {
    .need <- unlist(lapply(seq_along(.nms), function(.i) {
      if (.i == 1L) return(NULL)
      paste0("rxCor.", .nms[.i], ".", .nms[seq_len(.i - 1L)])
    }), use.names=FALSE)
    if (!all(.need %in% names(.fix))) return(NULL)
    .etaDistCorFromY(.nms, as.list(.fix[.need]))
  })
  names(.cor) <- vapply(.info$blocks, function(.n) .n[1], character(1),
                        USE.NAMES=FALSE)
  .cor <- .cor[!vapply(.cor, is.null, logical(1))]
  assign("etaDist", .info$etaDist, envir=ret$env)
  assign("etaDistCor", .cor, envir=ret$env)
  ret
}

postFinalObjectHooksAdd(".postFinalEtaDist", .postFinalEtaDist)
