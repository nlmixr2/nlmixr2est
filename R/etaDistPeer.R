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

#' d(log p_k)/d(THETA_j_), taken by symengine
#'
#' Deliberately the same shape as `.impmapChainRule()`, including the parameter
#' name `s`: the symengine environment has to be reached with the idiom that
#' setup supports, and reads have to go through the `$` accessor via an
#' intermediate variable.  `with(s, <bare name>)` does NOT work here -- the `$`
#' NSE misbehaves when nested in a call and base `get`/`::` are shadowed, which
#' is the same caveat `rxUiGet.impmapThetaSens()` records for rx_pred_ / rx_r_.
#'
#' No state-sensitivity term, unlike `.impmapChainRule()`: the peer is ODE-free
#' and a declaration reaching a state is declined before it gets here.
#'
#' @param s symengine environment carrying the peer model
#' @param target lhs name, eg `"rx_edll_1_"`
#' @param j 1-based theta index
#' @return the derivative as an rxode2 expression
#' @noRd
#' @author Matthew L. Fidler
.etaDistPeerD <- function(s, target, j) {
  .l <- eval(parse(text = paste0("with(s, D(", target, ", THETA_", j, "_))")))
  rxode2::rxFromSE(.l)
}

#' Peer model text for the declared-distribution M-step
#'
#' Builds the log density of every declared random effect and its derivative
#' with respect to every estimated theta, symbolically, in ONE symengine load
#' over the model's own text with the peer lines appended.
#'
#' Two constraints force that shape, both learned the hard way:
#'
#'   * The peer lines have to be IN the text that is loaded.  Defining
#'     `rx_edll_k_` into an environment loaded from the fitted model instead
#'     leaves symengine with a model that never mentions `llik*()`, and
#'     `rxFromSE()` then does not recognize llikGamma as a model function --
#'     it falls through to treating the `.rxD` entry's own R closure as a
#'     user-defined function ("R user function 'paste0' has variable number of
#'     arguments").  Putting the call in the text is what makes the `.rxD`
#'     dispatch fire.
#'
#'   * It has to be ONE load.  A symengine load shadows base `get`, `parse`,
#'     `paste0` and `with` on the search path, and a second `rxS()` in the same
#'     session leaves them shadowed, so every subsequent `with(s, D(...))`
#'     fails.  This runs inside the package namespace, where those names
#'     resolve through the namespace imports and the shadow is not on the
#'     lookup path -- which is also why `.impmapChainRule()` can use the same
#'     idiom, and why none of this works from the global environment.
#'
#' The peer is ODE-free (design section 4: "no states -- arithmetic on
#' parameters and covariates"), so the endpoint lines are dropped and no state
#' sensitivities are carried.  A declaration whose arguments reach a STATE is
#' declined by name rather than quietly given a wrong derivative.
#'
#' WHICH THETAS BELONG TO WHICH FAMILY is decided symbolically, not by scanning
#' the declaration for names: `d(rx_edll_k_)/d(THETA_j_)` is taken for every
#' estimated theta and the zero columns are dropped.  A theta that reaches the
#' family through a chain of intermediate model variables is found the same way
#' as one written into the declaration directly, and a theta that appears in
#' the text but cancels out is correctly left out.  A derivative that ERRORS is
#' not one that is zero -- that declines the family and names itself, rather
#' than silently narrowing what the M-step maximizes over.
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
  .ui <- rxode2::rxUiDecompress(x[[1]])
  .st <- .etaDistDeclGet(.ui)
  if (is.null(.st)) {
    .d <- rxode2::rxUiEtaDists(.ui)
    if (nrow(.d) == 0L) return(NULL)
    .st <- .etaDistDeclStash(.ui, .d)
    if (is.null(.st)) return(NULL)
  }
  .n <- length(.st$name)
  if (.n == 0L) return(NULL)
  .ini <- .ui$iniDf
  ## ETA[k]: the declared random effect's own latent slot, which the M-step
  ## overwrites with the fixed eta for the duration of the peer solve.
  ##
  ## rxEtaDistExpand() renames that latent `rxz.<eta>` and turns the declared
  ## name into an ordinary lhs, so THAT is what carries the eta number here --
  ## the same convention .etaDistMstepInfo() uses.  The bare name is the
  ## fallback for an UNEXPANDED ui, which is what the tests and etaDistInit()
  ## hand in.
  ##
  ## For a correlated pair the expansion also writes `rxN.<eta>`, the copula
  ## combination of the latents.  The peer wants `rxz.<eta>`: it never maps a
  ## latent through the inverse CDF, it only needs a slot to read the
  ## already-computed eta out of, and rxz is the one nothing else reads while
  ## the peer is swapped in.
  .netaOf <- function(.nm) {
    .w <- which(.ini$name == .nm & .ini$neta1 == .ini$neta2)
    if (length(.w) == 1L) as.integer(.ini$neta1[.w]) else NA_integer_
  }
  .latName <- character(.n)
  .etaIdx <- integer(.n)
  for (.k in seq_len(.n)) {
    .z <- paste0("rxz.", .st$name[.k])
    .i <- .netaOf(.z)
    if (is.na(.i)) { .z <- .st$name[.k]; .i <- .netaOf(.z) }
    .latName[.k] <- .z
    .etaIdx[.k] <- .i
  }
  ## Every estimated theta, NOT .impmapEstTheta()'s non-mu subset.  That split
  ## exists for the sensitivity model, whose mu-referenced thetas are updated by
  ## the EM closed form; it says nothing about which thetas a DECLARATION
  ## reaches, and on Bauer's gamma model it excludes all four of them.  The
  ## zero-column drop is what narrows this, and it narrows it correctly:
  ## rxCor.* and the residual-error thetas appear in no rx_edll_ and fall out.
  .th <- .ini[!is.na(.ini$ntheta), ]
  .idx <- sort(as.integer(.th$ntheta[!(!is.na(.th$fix) & .th$fix)]))
  .states <- .ui$state
  if (is.null(.states)) .states <- character(0)
  .declined <- character(0)
  .peerLine <- rep(NA_character_, .n)
  for (.k in seq_len(.n)) {
    .nm <- .st$name[.k]
    if (is.na(.etaIdx[.k])) {
      .declined[.nm] <- "no latent random effect to read the fixed eta from"
      next
    }
    .dens <- .etaDistPeerLogDensity(.st$etaDist[.k], .latName[.k], .nm)
    if (is.null(.dens)) {
      .declined[.nm] <- paste0("no log density for '",
                               as.character(str2lang(.st$etaDist[.k])[[1]]), "'")
      next
    }
    ## A distribution parameter is population-level plus covariates; it must not
    ## reach a state.  Checked on the DECLARATION as written -- once symengine
    ## substitutes, a state is indistinguishable from anything else it resolves.
    .reach <- intersect(all.vars(str2lang(.st$etaDist[.k])), .states)
    if (length(.reach) > 0L) {
      .declined[.nm] <- paste0("its arguments reach the state(s) ",
                               paste(.reach, collapse = ", "),
                               "; the peer is ODE-free and carries no state sensitivities")
      next
    }
    .peerLine[.k] <- paste0(.etaDistPeerLlName(.k), " <- ", .dens)
  }
  if (all(is.na(.peerLine))) {
    return(list(peer = NULL, etaIdx = .etaIdx,
                thetaIdx = vector("list", .n), declined = .declined))
  }
  ## ONE load, through rxUiGet.loadPruneSens(), of a ui carrying the peer lines.
  ##
  ## Not a bare rxS() on assembled text: loadPruneSens() is what sets the
  ## environment up the way the rest of this file's idioms require (it is the
  ## same entry rxUiGet.impmapThetaSens() uses), and a bare load leaves reads
  ## and derivatives failing inside rxFromSE().  And not two loads: a second
  ## rxS() in one session leaves base `get`/`parse`/`paste0`/`with` shadowed and
  ## every later derivative fails.
  ##
  ## The peer lines have to be IN the loaded model.  Defining rx_edll_k_ into an
  ## environment loaded WITHOUT them leaves symengine with a model that never
  ## mentions llik*(), and rxFromSE() then does not recognize llikGamma as a
  ## model function -- it falls through to treating the `.rxD` entry's own R
  ## closure as user-defined ("R user function 'paste0' has variable number of
  ## arguments").  Appending them is what makes the `.rxD` dispatch fire.
  .uiP <- rxode2::rxUiDecompress(.ui)
  .uiP$lstExpr <- c(.uiP$lstExpr,
                    lapply(.peerLine[!is.na(.peerLine)], str2lang))
  .s <- tryCatch(rxUiGet.loadPruneSens(list(.uiP), ...), error = function(e) NULL)
  if (!is.null(.s) && !exists("..maxTheta", .s)) .s <- NULL
  if (is.null(.s)) {
    for (.k in which(!is.na(.peerLine))) {
      .declined[.st$name[.k]] <- "the peer model could not be loaded into symengine"
    }
    return(list(peer = NULL, etaIdx = .etaIdx,
                thetaIdx = vector("list", .n), declined = .declined))
  }
  .lines <- character(0)
  .thetaIdx <- vector("list", .n)
  for (.k in which(!is.na(.peerLine))) {
    .nm <- .st$name[.k]
    .ll <- .etaDistPeerLlName(.k)
    ## `$` through an intermediate variable, never `with(s, <bare name>)` --
    ## see the caveat on .etaDistPeerD().
    .llSym <- .s[[.ll]]
    .llTxt <- tryCatch(rxode2::rxFromSE(.llSym), error = function(e) NULL)
    if (is.null(.llTxt)) {
      .declined[.nm] <- "the log density did not survive the symengine load"
      next
    }
    .keep <- integer(0); .err <- NULL; .col <- character(0)
    for (.j in .idx) {
      .g <- tryCatch(
        .etaDistPeerD(.s, .ll, .j),
        error = function(e) {
          .err <<- paste0("d/d(THETA_", .j, "_) could not be taken: ",
                          conditionMessage(e))
          NULL
        })
      if (!is.null(.err)) break
      if (.g %in% c("0", "0.0", "-0")) next
      .keep <- c(.keep, .j)
      .col <- c(.col, paste0(.etaDistPeerSensName(.k, .j), "=", .g))
    }
    if (!is.null(.err)) { .declined[.nm] <- .err; next }
    if (length(.keep) == 0L) {
      .declined[.nm] <- "no estimated theta reaches it, so there is nothing to maximize"
      next
    }
    .thetaIdx[[.k]] <- .keep
    .lines <- c(.lines, paste0(.ll, "=", .llTxt), .col)
  }
  ## Lightweight return only -- never the symengine environment (see the note on
  ## rxUiGet.impmapThetaSens, where caching it doubled memory per model).
  list(peer = if (length(.lines) == 0L) NULL else paste(c(.lines, ""), collapse = "\n"),
       etaIdx = .etaIdx, thetaIdx = .thetaIdx, declined = .declined)
}
attr(rxUiGet.etaDistPeer, "rstudio") <- emptyenv()

#' Compile the peer model
#'
#' Mirrors `.impmapThetaSensModel()`: same FOCEi codegen parameter block
#' (`THETA[]`/`ETA[]`), same cmt/interpolation preamble, same role-tagged
#' artifact name so it cannot share a compiled shared object with another build
#' of the same text.
#'
#' `eventSens` is not offered.  The peer has no states, so there is nothing for
#' an event sensitivity to be taken with respect to.
#'
#' @param ui rxode2 ui, already expanded
#' @return the compiled model, or `NULL` when there is no peer to compile
#' @noRd
#' @author Matthew L. Fidler
.etaDistPeerModel <- function(ui) {
  ## Through `$`, by name, so it is cached and `ui$etaDistPeer` stays printable
  ## while debugging -- the same reason .impmapThetaSensModel() does.
  .p <- ui$etaDistPeer
  if (is.null(.p) || is.null(.p$peer)) return(NULL)
  .cmt <- ui$foceiCmtPreModel
  .interp <- ui$interpLinesStr
  if (.interp != "") .cmt <- paste0(.cmt, "\n", .interp)
  nlmixr2global$toRxParam <-
    paste0(.uiGetThetaEtaParams(ui, TRUE), "\n", .cmt, "\n")
  nlmixr2global$toRxDvidCmt <- .foceiToCmtLinesAndDvid(ui)
  .toRx(.p$peer, "compiling declared-distribution model...",
        role = "rxEtaDistLl")
}

#' The saem-side plan for the declared-distribution peer
#'
#' The compiled peer plus the lhs NAMES the M-step has to resolve once it is
#' registered.  Names rather than offsets: the peer emits one block per family,
#' so there is no single contiguous offset to hand over, and `odeSwapLhsIndex()`
#' resolves each by name at registration.
#'
#' Shaped like `rxUiGet.saemThetaSensPlan()` and attached the same way, so the
#' peer rides the machinery that already exists rather than a parallel one.
#'
#' @param x rxode2 ui, in a list
#' @return list with `ok`, the compiled `etaDistLl`, the per-family lhs name
#'   vectors and the ETA slot each family reads, or `NULL` when the model
#'   declares no distribution the peer can score
#' @noRd
#' @author Matthew L. Fidler
#' @export
rxUiGet.etaDistPeerPlan <- function(x, ...) {
  .ui <- x[[1]]
  .p <- .ui$etaDistPeer
  if (is.null(.p) || is.null(.p$peer)) return(NULL)
  .mod <- tryCatch(.etaDistPeerModel(.ui), error = function(e) NULL)
  if (is.null(.mod)) return(NULL)
  ## Only the families that actually produced columns.  A declined family has
  ## no lhs to resolve and must not occupy a slot in the M-step's index
  ## vectors, or the k-th entry would stop meaning the k-th declared eta.
  .keep <- which(!vapply(.p$thetaIdx, is.null, logical(1)))
  if (length(.keep) == 0L) return(NULL)
  list(ok = TRUE,
       etaDistLl = .mod,
       ## which declared random effect each retained block belongs to (1-based)
       etaDistLlFam = as.integer(.keep),
       ## ETA[k] slot each retained family reads its fixed eta out of
       etaDistLlEta = as.integer(.p$etaIdx[.keep]),
       ## lhs name of each retained family's log density
       etaDistLlName = vapply(.keep, .etaDistPeerLlName, character(1)),
       ## the theta indices behind each retained family's gradient columns,
       ## flattened with a per-family count so C++ can walk them without a
       ## ragged structure
       etaDistLlNth = as.integer(vapply(.p$thetaIdx[.keep], length, integer(1))),
       etaDistLlTheta = as.integer(unlist(.p$thetaIdx[.keep], use.names = FALSE)),
       etaDistLlGradName = unlist(lapply(.keep, function(.k) {
         vapply(.p$thetaIdx[[.k]], function(.j) .etaDistPeerSensName(.k, .j),
                character(1))
       }), use.names = FALSE),
       declined = .p$declined)
}
attr(rxUiGet.etaDistPeerPlan, "rstudio") <- emptyenv()
