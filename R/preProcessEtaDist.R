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

#' The `"etaDist"` attribute of the dispatched estimation method
#'
#' Read from the `nlmixr2Est.<method>` S3 method, so a method registered
#' by another package can declare support without editing this file.  The
#' attribute may be `TRUE`/`FALSE`, the string `"native"`, or a
#' `function(control)` returning one of those.
#'
#' `"native"` means the method translates the declaration ITSELF and must
#' see it unexpanded -- babelmixr2's `est="nonmem"` writes Bauer's own
#' `$ABBR FUNCTION GAMMACDFINV` control stream, which reads nothing like
#' the expansion and could not be recovered from it.
#'
#' @param est estimation routine name
#' @param control control object
#' @return `TRUE`, `FALSE` or `"native"`
#' @noRd
#' @author Matthew L. Fidler
.etaDistMethodAttr <- function(est, control=NULL) {
  if (!is.character(est) || length(est) != 1L) return(FALSE)
  .v <- as.character(utils::methods("nlmixr2Est"))
  if (!(paste0("nlmixr2Est.", est) %in% .v)) return(FALSE)
  .a <- attr(utils::getS3method("nlmixr2Est", est), "etaDist")
  if (is.null(.a)) return(FALSE)
  if (is.function(.a)) .a <- .a(control)
  if (identical(.a, "native")) return("native")
  isTRUE(.a)
}

#' Does the dispatched estimation method support a declared eta distribution?
#'
#' @param est estimation routine name
#' @param control control object
#' @return boolean
#' @noRd
#' @author Matthew L. Fidler
.isEtaDistMethod <- function(est, control=NULL) {
  !isFALSE(.etaDistMethodAttr(est, control))
}

#' The error a method that cannot use a declared distribution gets
#'
#' One message, raised from two places: the pre-processing hook (early,
#' before any work is done) and the gate in `nlmixr2Est()` (the backstop).
#'
#' @param d declaring random effects, as `rxUiEtaDists()` returns them
#' @param est estimation routine name
#' @param control control object
#' @return nothing, called for the error
#' @noRd
#' @author Matthew L. Fidler
.etaDistRefuse <- function(d, est, control) {
  if (nrow(d) == 0L) return(invisible())
  if (!is.character(est) || length(est) != 1L) return(invisible())
  if (.isEtaDistMethod(est, control)) return(invisible())
  stop("est=\"", est, "\" cannot use the declared non-normal random effect ",
       "distribution(s) on '", paste(d$name, collapse="', '"), "'",
       call.=FALSE)
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
  .etaDistRefuse(.d, .est, .control)
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
  ## Refuse HERE rather than leaving it to the gate in nlmixr2Est().  The
  ## hooks run first, so a method that cannot use a declared distribution
  ## would otherwise pay for the expansion and its own pre-processing
  ## before being told no -- for `est="npag"` that is the whole
  ## nonparametric mu-expansion, which is not a wait to impose on someone
  ## who is about to get an error.  The gate stays as the backstop for
  ## paths that reach nlmixr2Est() without running hooks.
  .etaDistRefuse(.d, est, control)
  ## a method that translates the declaration itself has to see it
  ## unexpanded (see `.etaDistMethodAttr()`)
  if (identical(.etaDistMethodAttr(est, control), "native")) return(NULL)
  list(ui=rxode2::rxEtaDistExpand(ui),
       control=.etaDistMceta(control, est))
}

#' Default the inner-problem multistart on for a declared distribution
#'
#' These are harder likelihood surfaces than the models nlmixr2 usually
#' sees.  Measured on Bauer's own gamma dataset (300 subjects), a cold
#' started FOCEi settles into a local optimum with the residual error
#' absorbing variability that belongs to the random effect -- the gamma
#' mean of CL comes back as 6.7 against a truth of 5.0, and its relative
#' variance as 0.22 against 0.085.  With `mceta = 10` the same cold start
#' gives 5.2 and 0.145.  That is exactly the knob Bauer's own recipe uses
#' (`MCETA=10`, and `100` for the harder cases).
#'
#' `mceta = n` evaluates the previous eta, eta = 0, and `n - 1` draws from
#' omega for each subject's inner problem and keeps the best.  It suits a
#' declared distribution particularly well: the latent omega is the
#' identity, so those draws are standard normal -- precisely the right
#' proposal on the scale the inner problem actually works on.
#'
#' Only the methods that HAVE the knob are touched -- the FOCEi family,
#' `imp`/`impmap` and `agq`.  `saemControl()` has no `mceta` (SAEM samples
#' the etas rather than optimizing them), so a saem control is returned
#' unchanged and nothing is said about it.
#'
#' The default is only changed when it is still the package default, and
#' the change is announced, so an explicit `mceta =` stays in charge.
#'
#' @param control control object, possibly NULL
#' @param est estimation routine name, used only to name the right control
#'   function in the message
#' @return the control with `mceta` defaulted on, or NULL when this control
#'   has no `mceta` or already carries a chosen one
#' @noRd
#' @author Matthew L. Fidler
.etaDistMceta <- function(control, est = NULL) {
  if (!is.list(control) || !any(names(control) == "mceta")) return(NULL)
  if (!identical(as.integer(control$mceta), .etaDistMcetaPkgDefault)) return(NULL)
  control$mceta <- .etaDistMcetaDefault
  .ctlFn <- if (is.character(est) && length(est) == 1L &&
                  exists(paste0(est, "Control"),
                         envir=asNamespace("nlmixr2est"), inherits=FALSE)) {
    paste0(est, "Control(mceta=)")
  } else {
    "mceta= in the control"
  }
  .minfo(paste0("declared random effect distribution: using mceta=",
                .etaDistMcetaDefault, " (set ", .ctlFn, " to override)"))
  control
}

## the foceiControl() default, and what a declared distribution uses instead
.etaDistMcetaPkgDefault <- -2L
.etaDistMcetaDefault <- 10L

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

#' Report the declared distributions on a fit
#'
#' Computed on demand through the `nmObjGet` accessors rather than stored
#' by a post-estimation hook, so they are there for every method that can
#' fit such a model -- the post-final hooks only run on the FOCEi path.
#'
#' `$etaDist` is the declarations as the model wrote them; `$etaDistCor`
#' is the copula correlation matrix of each declared block, rebuilt from
#' the `rxCor.*` estimates.
#'
#' The `rxCor.*` rows already read as correlations in the fit's
#' back-transformed column: `rxEtaDistExpand()` gives them
#' `backTransform("tanh")`, and `tanh()` of one is the partial correlation
#' between its two random effects given the ones before them -- which for
#' a 2x2 block, the usual case, is simply the correlation.
#'
#' @param x list of the fit environment and the exact flag, as
#'   `nmObjGet()` dispatches it
#' @param ... ignored
#' @return the declarations, or NULL when the model declared none
#' @export
#' @keywords internal
#' @author Matthew L. Fidler
nmObjGet.etaDist <- function(x, ...) {
  .info <- .etaDistInfo(x[[1]])
  if (is.null(.info)) return(NULL)
  .info$etaDist
}
attr(nmObjGet.etaDist, "desc") <-
  "The non-normal random effect distributions the model declared"

#' @rdname nmObjGet.etaDist
#' @export
nmObjGet.etaDistCor <- function(x, ...) {
  .env <- x[[1]]
  .info <- .etaDistInfo(.env)
  if (is.null(.info)) return(NULL)
  .fix <- try(get("fixef", envir=.env), silent=TRUE)
  if (inherits(.fix, "try-error") || is.null(.fix)) return(NULL)
  .cor <- lapply(.info$blocks, function(.nms) {
    .need <- unlist(lapply(seq_along(.nms), function(.i) {
      if (.i == 1L) return(NULL)
      paste0("rxCor.", .nms[.i], ".", .nms[seq_len(.i - 1L)])
    }), use.names=FALSE)
    if (length(.need) == 0L || !all(.need %in% names(.fix))) return(NULL)
    .etaDistCorFromY(.nms, as.list(.fix[.need]))
  })
  names(.cor) <- vapply(.info$blocks, function(.n) .n[1], character(1),
                        USE.NAMES=FALSE)
  .cor <- .cor[!vapply(.cor, is.null, logical(1))]
  if (length(.cor) == 0L) return(NULL)
  .cor
}
attr(nmObjGet.etaDistCor, "desc") <-
  "The Gaussian copula correlation of each declared random effect block"

#' What `rxEtaDistExpand()` recorded on the fit's model
#'
#' @param env fit environment
#' @return the `etaDistInfo` list, or NULL
#' @noRd
#' @author Matthew L. Fidler
.etaDistInfo <- function(env) {
  .ui <- try(get("ui", envir=env), silent=TRUE)
  if (inherits(.ui, "try-error") || is.null(.ui)) return(NULL)
  .ui <- try(rxode2::rxUiDecompress(.ui), silent=TRUE)
  if (inherits(.ui, "try-error")) return(NULL)
  .info <- try(get("etaDistInfo", envir=.ui), silent=TRUE)
  if (inherits(.info, "try-error")) return(NULL)
  .info
}
