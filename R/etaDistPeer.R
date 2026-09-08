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

#' The lhs name carrying d(log p_k)/d(THETA_j_)
#' @param k 1-based index of the declared random effect
#' @param j 1-based theta index (THETA_j_ / ntheta ordering)
#' @return character
#' @noRd
#' @author Matthew L. Fidler
.etaDistPeerSensName <- function(k, j) {
  paste0("rx__sens_rx_edll_", k, "___BY_THETA_", j, "___")
}

#' Peer model text for the declared-distribution M-step
#'
#' Built the way `rxUiGet.impmapThetaSens()` builds the sensitivity model:
#' load the model into symengine, define the new quantities there so the
#' argument expressions resolve against the model's own assignments, then
#' render each one back with `rxFromSE()`.
#'
#' Two differences from the sensitivity model, both because this peer is
#' ODE-free (section 4: "The peer has no states -- arithmetic on parameters and
#' covariates"):
#'
#'   * no `..ddt` / `..sens` lines are carried, so no state sensitivity ODEs and
#'     no integration when it is swapped in;
#'   * a declaration whose arguments reach a STATE is declined by name rather
#'     than quietly given a wrong derivative -- the peer would then need the
#'     state sensitivities it deliberately does not carry.
#'
#' WHICH THETAS BELONG TO WHICH FAMILY is decided symbolically, not by scanning
#' the declaration for names: `d(rx_edll_k_)/d(THETA_j_)` is emitted for every
#' estimated theta and the zero columns are dropped.  A theta that reaches the
#' family through a chain of intermediate model variables is found the same way
#' as one written into the declaration directly, and a theta that appears in the
#' text but cancels out is correctly left out.
#'
#' @param x rxode2 ui, in a list (rxUiGet convention)
#' @return a list with `peer` (the model text), `etaIdx` (the `ETA[k]` slot each
#'   declared random effect reads its fixed value from), `thetaIdx` (a list, per
#'   family, of the theta indices whose derivative is non-zero) and `declined`
#'   (a named character vector of families the peer will not score, with the
#'   reason), or `NULL` when the model declares nothing
#' @noRd
#' @author Matthew L. Fidler
#' @export
rxUiGet.etaDistPeer <- function(x, ...) {
  .ui <- x[[1]]
  .core <- .etaDistMstepCore(.ui)
  if (is.null(.core)) return(NULL)
  .st <- .etaDistDeclGet(.ui)
  if (is.null(.st)) return(NULL)
  .ini <- rxode2::rxUiDecompress(.ui)$iniDf
  .n <- length(.st$name)
  ## ETA[k]: the declared random effect's own latent slot, which the M-step
  ## overwrites with the fixed eta for the duration of the peer solve.
  .etaIdx <- vapply(.st$name, function(.nm) {
    .w <- which(.ini$name == .nm & .ini$neta1 == .ini$neta2)
    if (length(.w) == 1L) as.integer(.ini$neta1[.w]) else NA_integer_
  }, integer(1), USE.NAMES = FALSE)
  .idx <- .impmapEstTheta(.ui)$all
  .s <- rxUiGet.loadPruneSens(x, ...)
  if (!exists("..maxTheta", .s)) return(NULL)
  .stateVars <- .rxode2stateOdeNoOutput(.s)
  .declined <- character(0)
  .lines <- character(0)
  .thetaIdx <- vector("list", .n)
  for (.k in seq_len(.n)) {
    .nm <- .st$name[.k]
    if (is.na(.etaIdx[.k])) {
      .declined[.nm] <- "no latent random effect to read the fixed eta from"
      next
    }
    ## The symengine SYMBOL is ETA_k_ (it renders back as ETA[k]), the same
    ## spelling .impmapChainRule() uses for THETA_j_.
    .d <- .etaDistPeerLogDensity(.st$etaDist[.k], paste0("ETA_", .etaIdx[.k], "_"),
                                 .nm)
    if (is.null(.d)) {
      .declined[.nm] <- paste0("no log density for '",
                               as.character(str2lang(.st$etaDist[.k])[[1]]), "'")
      next
    }
    ## A distribution parameter is population-level plus covariates; it must not
    ## reach a state.  Checked on the DECLARATION as written -- after symengine
    ## substitutes, a state is indistinguishable from anything else it resolves.
    .reach <- all.vars(str2lang(.st$etaDist[.k]))
    if (length(intersect(.reach, .stateVars)) > 0L) {
      .declined[.nm] <- paste0("its arguments reach the state(s) ",
                               paste(intersect(.reach, .stateVars), collapse = ", "),
                               "; the peer is ODE-free and carries no state sensitivities")
      next
    }
    .ll <- .etaDistPeerLlName(.k)
    .ok <- tryCatch({
      .e <- eval(parse(text = paste0("with(.s, ", .d, ")")))
      assign(.ll, .e, envir = .s)
      TRUE
    }, error = function(e) {
      .declined[.nm] <<- paste0("could not be built: ", conditionMessage(e))
      FALSE
    })
    if (!.ok) next
    .lines <- c(.lines, paste0(.ll, "=", rxode2::rxFromSE(get(.ll, envir = .s))))
    ## The family's thetas, found by differentiating rather than by name.
    .keep <- integer(0)
    for (.j in .idx) {
      .g <- tryCatch(
        rxode2::rxFromSE(eval(parse(text = paste0("with(.s, D(", .ll, ", THETA_",
                                                  .j, "_))")))),
        error = function(e) "0")
      if (.g %in% c("0", "0.0", "-0")) next
      .keep <- c(.keep, .j)
      .lines <- c(.lines, paste0(.etaDistPeerSensName(.k, .j), "=", .g))
    }
    .thetaIdx[[.k]] <- .keep
    if (length(.keep) == 0L) {
      .declined[.nm] <- "no estimated theta reaches it, so there is nothing to maximize"
    }
  }
  if (length(.lines) == 0L) {
    return(list(peer = NULL, etaIdx = .etaIdx, thetaIdx = .thetaIdx,
                declined = .declined))
  }
  ## Lightweight return only -- never the symengine environment (see the note on
  ## rxUiGet.impmapThetaSens, where caching `.s` doubled memory per model).
  list(peer = paste(c(.lines, ""), collapse = "\n"), etaIdx = .etaIdx,
       thetaIdx = .thetaIdx, declined = .declined)
}
attr(rxUiGet.etaDistPeer, "rstudio") <- emptyenv()
