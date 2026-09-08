## The peer model for the declared-distribution M-step.
##
## Design: R/etaDistMstep.R, section 4.  ONE rxode2 model, swapped in with
## odeSwap the way odeSlotPred / odeSlotThetaSens / odeSlotHess2 already are,
## whose lhs give the M-step, PER OBSERVATION RECORD,
##
##     rx_edll_<k>_                            log p( eta_k ; args_k )
##     rx__sens_rx_edll_<k>___BY_THETA_<j>___  d/d(THETA_j_) of it
##
## split by declared random effect, because each has its own family and its own
## thetas and the only coupling between them is the copula (section 5,
## estimated separately and in closed form).
##
## NOTHING HERE DIFFERENTIATES ANYTHING.  rxode2's symengine does the whole
## chain rule -- d(log p)/d(args) from the `.rxD` entries that already back
## llikGamma()/llikNorm()/..., and d(args)/d(theta) from the argument
## expressions as written.  Verified on all 22 declarable families:
## `D(llikGamma(x, sh, rt), sh)` comes back as `llikGammaDshape(x, sh, rt)`,
## and composes through `sh = 1/exp(THETA_2_)` by itself.  The families with no
## llik() counterpart (gumbel, frechet, pareto, ...) have elementary log
## densities that symengine differentiates natively; `lgammafn` differentiates
## to `digamma`.
##
## A covariate on a distribution parameter therefore needs no handling at all.
## The peer is evaluated per record against the data, so fixed or time-varying
## is simply a value the model has.
##
## HOW THE ETA GETS IN.  The M-step holds eta fixed at Q(phiU(w_i); theta_old)
## and moves theta only inside args -- if the peer recomputed eta from the
## candidate args the objective would collapse to the latent normal density and
## be constant in theta.  So the eta has to enter as a VALUE, and the peer
## reads it out of the declared eta's own latent slot, `ETA[k]`: the M-step
## writes the fixed eta there for the duration of the peer solve.  The peer
## never calls phiU() or the inverse CDF, so nothing else reads that slot while
## it is swapped in, and symengine sees `ETA[k]` as a free symbol -- which is
## exactly the fixed-eta derivative that is wanted.  This also keeps the peer
## on the shared `ind->par_ptr` layout, so it introduces no new parameter.

#' Log-density templates for the families that can be declared on an eta
#'
#' Each entry is `name|template`, where the template is an 'rxode2' expression
#' with `{x}` standing for the random effect and `{par}` for the distribution's
#' parameter of that name -- the same shape, and the same parameter names in
#' the same canonical order, as lotri's inverse-CDF templates, so
#' `lotriEtaDists()` stays the single source of truth for the parameterization.
#'
#' This is the OTHER half of what lotri supplies.  lotri gives the quantile
#' function, which builds the eta going forward; this gives the log density,
#' which scores it in the M-step.  Neither replaces the other.
#'
#' `llik*()` is used wherever rxode2 has one, because those carry exact,
#' numerically careful parameter derivatives.  The rest are written out; every
#' one of them is elementary.
#'
#' `stdNormal` has no parameters, so it has no M-step and no entry.
#'
#' @noRd
#' @author Matthew L. Fidler
.etaDistPeerDefs <- c(
  ## unbounded continuous
  "dnorm|llikNorm({x}, {mean}, {sd})",
  "studentT|llikT({x}, {nu}, {mu}, {sigma})",
  "dcauchy|llikCauchy({x}, {location}, {scale})",
  "doubleExponential|-log(2*({sigma})) - abs(({x}) - ({mu}))/({sigma})",
  "dlogis|-log({scale}) - (({x}) - ({location}))/({scale}) - 2*log(1 + exp(-(({x}) - ({location}))/({scale})))",
  "gumbel|-log({beta}) - (({x}) - ({mu}))/({beta}) - exp(-(({x}) - ({mu}))/({beta}))",
  ## positive continuous
  "dlnorm|llikNorm(log({x}), {meanlog}, {sdlog}) - log({x})",
  "dchisq|llikChisq({x}, {df})",
  "invChiSquare|-({nu})/2*log(2) - lgammafn(({nu})/2) - (({nu})/2 + 1)*log({x}) - 1/(2*({x}))",
  "scaledInvChiSquare|({nu})/2*log(({nu})*({sigma})*({sigma})/2) - lgammafn(({nu})/2) - (({nu})/2 + 1)*log({x}) - ({nu})*({sigma})*({sigma})/(2*({x}))",
  "dexp|llikExp({x}, {rate})",
  "dgamma|llikGamma({x}, {shape}, {rate})",
  "invGamma|({alpha})*log({beta}) - lgammafn({alpha}) - (({alpha}) + 1)*log({x}) - ({beta})/({x})",
  "dweibull|llikWeibull({x}, {shape}, {scale})",
  "frechet|log({alpha}) - log({sigma}) - (1 + ({alpha}))*log(({x})/({sigma})) - (({x})/({sigma}))^(-({alpha}))",
  "rayleigh|log({x}) - 2*log({sigma}) - ({x})*({x})/(2*({sigma})*({sigma}))",
  "pareto|log({alpha}) + ({alpha})*log({y_min}) - (({alpha}) + 1)*log({x})",
  "paretoType2|log({alpha}) - log({lambda}) - (({alpha}) + 1)*log(1 + (({x}) - ({mu}))/({lambda}))",
  ## bounded continuous
  "dbeta|llikBeta({x}, {shape1}, {shape2})",
  "betaProportion|llikBeta({x}, ({mu})*({kappa}), (1 - ({mu}))*({kappa}))",
  "dunif|llikUnif({x}, {min}, {max})")

.etaDistPeerCache <- new.env(parent=emptyenv())

#' The peer table: lotri's eta distributions plus a `logDensity` column
#'
#' Built on first use and cached.  A family lotri can declare but this has no
#' template for gets `NA`, which is what makes the peer decline for it -- with
#' a named reason, not silently.
#'
#' @return data frame with the columns of `lotri::lotriEtaDists()` plus
#'   `logDensity`
#' @noRd
#' @author Matthew L. Fidler
.etaDistPeerTable <- function() {
  if (!is.null(.etaDistPeerCache$tab)) return(.etaDistPeerCache$tab)
  .tab <- lotri::lotriEtaDists()
  .l <- strsplit(.etaDistPeerDefs, "|", fixed=TRUE)
  .nm <- vapply(.l, `[[`, character(1), 1L, USE.NAMES=FALSE)
  .ld <- vapply(.l, `[[`, character(1), 2L, USE.NAMES=FALSE)
  .w <- match(.nm, .tab$name)
  if (anyNA(.w)) {
    stop("eta log-density(s) not in the lotri eta distribution table: '", # nocov
         paste(.nm[is.na(.w)], collapse="', '"), "'", call.=FALSE) # nocov
  }
  .tab$logDensity <- NA_character_
  .tab$logDensity[.w] <- .ld
  rownames(.tab) <- NULL
  .etaDistPeerCache$tab <- .tab
  .tab
}

#' Substitute a declaration into its log-density template
#'
#' Mirrors rxode2's `.rxEtaDistQuantile()`: the arguments are positional in
#' lotri's canonical order, so the template's `{name}` placeholders line up by
#' position.
#'
#' @param txt the declaration as stored, eg `"dgamma(1/exp(lrv), 1/(exp(lrv)*exp(lm)))"`
#' @param x rxode2 expression for the random effect, eg `"ETA[1]"`
#' @param what the eta's name, for messages
#' @return an rxode2 expression, or NULL when the family has no log density
#'   here (the caller reports the decline)
#' @noRd
#' @author Matthew L. Fidler
.etaDistPeerLogDensity <- function(txt, x, what) {
  .call <- str2lang(txt)
  .nm <- as.character(.call[[1]])
  .tab <- .etaDistPeerTable()
  .w <- which(.tab$name == .nm)
  if (length(.w) != 1L || is.na(.tab$logDensity[.w])) return(NULL)
  .d <- .tab$logDensity[.w]
  .args <- as.list(.call)[-1]
  .parNames <- character(0)
  if (nzchar(.tab$parNames[.w])) {
    .parNames <- strsplit(.tab$parNames[.w], ",", fixed=TRUE)[[1]]
  }
  for (.i in seq_along(.args)) {
    .d <- gsub(paste0("{", .parNames[.i], "}"),
               paste0("(", deparse1(.args[[.i]]), ")"), .d, fixed=TRUE)
  }
  .d <- gsub("{x}", x, .d, fixed=TRUE)
  if (grepl("{", .d, fixed=TRUE)) {
    stop("'", what, "' does not supply every argument of '", .nm, "'", # nocov
         call.=FALSE) # nocov
  }
  .d
}

#' The lhs name carrying a declared family's log-likelihood
#' @param k 1-based index of the declared random effect
#' @return character
#' @noRd
#' @author Matthew L. Fidler
.etaDistPeerLlName <- function(k) paste0("rx_edll_", k, "_")
