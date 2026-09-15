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
## covariance update), and simulation.
## Refused: `npag`/`npb`, which model the random effect distribution
## nonparametrically, so a declared one contradicts them outright; `nlme`
## and `nls`, which are Gaussian by construction; `emvi`/`fbvi`, each of
## which needs its own audit before the same claim can be made; and `vae`,
## which was allowed on a structural argument and is refused again on a
## measured one -- see the block above `attr(nlmixr2Est.vae, ...)` in
## R/vae.R.  Its ELBO relocates the between-subject variability into the
## residual instead of estimating it, on every arm and at every dispersion.

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
#' Warn when a covariate coefficient on a declaration starts at exactly zero
#'
#' A covariate coefficient conventionally starts at 0, and on a declared
#' distribution that is the one starting value the outer search cannot leave.
#' Measured on a known allometric effect (true coefficient 0.75, gamma CL with
#' the covariate in its rate): focei returns 0.0017 from a start of 0 and
#' 0.7353 from a start of 0.1.  The objective is not flat there -- its central
#' difference at 0 is -53.4 -- but from 0 the search takes one step of about
#' 1e-3, the objective change falls under tolerance, and it stops.  Every start
#' from 0.1 to 0.5 converges to the same optimum, nine objective units better.
#'
#' Warning rather than adjusting the value: the start is the user's, and a
#' coefficient nudged off zero behind their back would change a reported
#' estimate with no record of why.
#'
#' @param d declaration table from `.rxUiEtaDists()`
#' @param ui rxode2 ui
#' @return nothing, called for the warning
#' @noRd
.etaDistWarnZeroSlope <- function(d, ui, est = "") {
  .iniDf <- try(ui$iniDf, silent = TRUE)
  if (inherits(.iniDf, "try-error") || is.null(.iniDf)) return(invisible())
  .th <- .iniDf[!is.na(.iniDf$ntheta), c("name", "est"), drop = FALSE]
  if (nrow(.th) == 0L) return(invisible())
  .hit <- character(0)
  for (.i in seq_len(nrow(d))) {
    .cl <- try(str2lang(d$etaDist[.i]), silent = TRUE)
    if (inherits(.cl, "try-error")) next
    .v <- all.vars(.cl)
    ## a symbol in the declaration that is not a theta is a covariate; with no
    ## covariate there is no slope to be trapped
    if (length(setdiff(.v, .th$name)) == 0L) next
    .own <- intersect(.v, .th$name)
    .z <- .own[.th$est[match(.own, .th$name)] == 0]
    if (length(.z)) {
      .hit <- c(.hit, paste0(.z, " (in dist(", d$name[.i], "))"))
    }
  }
  # Three estimators, three different behaviors, all measured on a known
  # allometric effect (true coefficient 0.75, covariate in a gamma rate) from
  # displaced starts.  Do not hand any of them another's advice.
  #
  #   focei  recovers it, given a non-zero slope start -- which
  #          foceiControl(zeroThetaRetry=) now arranges.  Falls through below.
  #   saem   ESTIMATES it, on both routes.  This said the opposite -- "estimates
  #          it with NOTHING", the thetas back BIT-EXACTLY at ini(), blamed on a
  #          coefficient having "no phi0 column to write back through".  The
  #          coefficient always had a column; rxode2's mu2 scan was claiming it
  #          out of the `rxEdA.*` role anchor and moving it into COV/MCOV, which
  #          is right on the cdf route (the decoder reads that anchor) and wrong
  #          on the direct route (the anchor feeds the prior alone).  rxode2's
  #          `.muRefAnchorIsDirect()` now separates the two.  Measured, true
  #          coefficient 0.75 throughout:
  #
  #            direct route, start 0.30              -> 0.7441
  #            cdf, subject-constant cov, start 0.35 -> 0.5526
  #            cdf, time-varying cov,     start 0.35 -> 0.4299
  #
  #          On the time-varying arm focei gets 0.4313 on the same data, so the
  #          two agree to 0.3% there.  The subject-constant cdf arm is the one
  #          that lands furthest from truth, which is what the warning is for.
  #   imp    moves them -- lclm 1.9 -> 1.71 and lclrv -2.0 -> -2.25, both
  #          toward truth -- but the coefficient stalls at 0.032 against 0.75.
  #          It gets neither the zeroTheta nudge nor the retry, both of which
  #          are FOCEi-family only.
  if (!identical(est, "") && !grepl("^(focei|foce|foi|posthoc)", est)) {
    .cov <- character(0)
    for (.i in seq_len(nrow(d))) {
      .cl <- try(str2lang(d$etaDist[.i]), silent = TRUE)
      if (inherits(.cl, "try-error")) next
      if (length(setdiff(all.vars(.cl), .th$name)) > 0L) {
        .cov <- c(.cov, paste0("dist(", d$name[.i], ")"))
      }
    }
    if (length(.cov) > 0L) {
      .what <- if (grepl("^saem", est)) {
        paste0("estimates them, but check the coefficient against its ",
               "standard error before trusting its magnitude.  Measured on a ",
               "120-subject arm with a true coefficient of 0.75: the direct ",
               "route (saemControl(etaDistParam=\"direct\")) returned 0.7441 ",
               "from a start of 0.30, while the cdf route returned 0.5526 ",
               "(subject-constant covariate) and 0.4299 (time-varying) from a ",
               "start of 0.35 -- focei gets 0.4313 on that last arm, so the ",
               "two agree there and differ most on the subject-constant one")
      } else if (grepl("^(imp|impmap|qrpem)", est)) {
        paste0("moves them but does not reach the coefficient: measured, the ",
               "declaration's own thetas travel toward truth while the ",
               "coefficient stalls near its starting value (0.032 against a ",
               "true 0.75).  The zeroTheta nudge and its retry, which fix ",
               "this for FOCEi, are FOCEi-family only")
      } else {
        # deliberately makes no measured claim: only focei, saem and imp have
        # been run on this arm, and quoting one of their numbers here would
        # attribute a measurement to an estimator that never produced it
        paste0("has not been verified to estimate a covariate on a ",
               "declaration.  The zeroTheta nudge and its retry, which make ",
               "this work for FOCEi, are FOCEi-family only, and saem is known ",
               "not to estimate such parameters at all.  Check the ",
               "coefficient against its starting value and standard error ",
               "before trusting it")
      }
      ## The FOCEi recommendation is not advice to give saem any more -- saem
      ## estimates these now -- so it is per-branch rather than appended to
      ## every message.
      .use <- if (grepl("^saem", est)) {
        ""
      } else {
        paste0("  Use est=\"focei\" for a covariate on a declaration, where ",
               "the coefficient is recovered.")
      }
      warning("est=\"", est, "\" ", .what, ": ",
              paste(.cov, collapse = ", "), ".", .use, call. = FALSE)
      return(invisible())
    }
  }
  if (length(.hit) == 0L) return(invisible())
  warning("a covariate-carrying declaration has a parameter starting at ",
          "exactly 0: ", paste(.hit, collapse = ", "),
          ".  A zero start carries no magnitude for the outer search to scale ",
          "by, so such a coefficient can come back sitting on its ",
          "foceiControl(zeroTheta=) nudge rather than estimated.  The FOCEi ",
          "family detects that and re-fits once from a larger nudge ",
          "(foceiControl(zeroThetaRetry=)), keeping the better objective, so ",
          "no action is usually needed -- but check the estimate against its ",
          "relative standard error, which is extreme for a coefficient that ",
          "is still on the nudge.  One thing the retry cannot fix: under ",
          "prop() alone the objective is not continuous in such a ",
          "coefficient, so bound the residual variance with a FIXED add() ",
          "alongside it.", call. = FALSE)
  invisible()
}

#' Warn when a covariate enters both a declaration and the structural model
#'
#' A covariate on a declaration's role anchor and the SAME covariate on the
#' structural expression that consumes that random effect are two ways of
#' saying the same thing about the same parameter.  When they enter with the
#' same functional form the two coefficients are not separately identifiable --
#' only their sum is -- and the fit will happily report both, trading them off
#' against each other run to run.
#'
#' Warns rather than refusing, deliberately, even though the plan called for
#' refusing the exactly-aliased case.  Deciding "exactly aliased" needs the two
#' functional forms compared, and doing that on expression text is guesswork:
#' `log(WT/70)` against `log(WT)-log(70)` is the same shape spelled two ways,
#' while a genuinely different shape can look similar.  A false warning costs a
#' line of output; a false refusal blocks a model the user cannot then fit at
#' all.  The message names both locations so the reader can judge.
#'
#' @param d declaration table from `.rxUiEtaDists()`
#' @param ui rxode2 ui (unexpanded -- this runs before the expansion)
#' @return nothing, called for the warning
#' @noRd
.etaDistWarnCovAliased <- function(d, ui) {
  .cov <- try(ui$allCovs, silent = TRUE)
  if (inherits(.cov, "try-error") || length(.cov) == 0L) return(invisible())
  .lst <- try(ui$lstExpr, silent = TRUE)
  if (inherits(.lst, "try-error") || length(.lst) == 0L) return(invisible())
  .hit <- character(0)
  for (.i in seq_len(nrow(d))) {
    .cl <- try(str2lang(d$etaDist[.i]), silent = TRUE)
    if (inherits(.cl, "try-error")) next
    .declCov <- intersect(all.vars(.cl), .cov)
    if (length(.declCov) == 0L) next
    .eta <- d$name[.i]
    for (.e in .lst) {
      .v <- all.vars(.e)
      if (!(.eta %in% .v)) next
      ## the line that ASSIGNS the eta is the declaration's own decoder, not a
      ## structural use of it
      if (is.call(.e) && length(.e) > 2L && identical(.e[[1]], quote(`<-`)) &&
            is.name(.e[[2]]) && identical(as.character(.e[[2]]), .eta)) {
        next
      }
      .both <- intersect(.declCov, .v)
      if (length(.both) > 0L) {
        .hit <- c(.hit, paste0(paste(.both, collapse = ", "),
                               " (on dist(", .eta, ") and on `",
                               deparse1(.e), "`)"))
      }
    }
  }
  if (length(.hit) == 0L) return(invisible())
  warning("a covariate enters both a declared distribution and the structural ",
          "model for the same parameter: ", paste(unique(.hit), collapse = "; "),
          ".  If the two enter with the same functional form their ",
          "coefficients are not separately identifiable -- only their sum is -- ",
          "and a fit will trade them off against each other.  Keep the ",
          "covariate in one place, or check that the two forms really are ",
          "different.", call. = FALSE)
  invisible()
}

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
  .etaDistWarnZeroSlope(.d, ui, est)
  .etaDistWarnCovAliased(.d, ui)
  ## Warm start, before the expansion and after the refusal.
  ##
  ## A declared family is very largely a STARTING VALUE problem: the E-step
  ## hands the M-step an eta sample drawn under the current family, the M-step
  ## fits the family to that sample, and from a poor start the pair simply walk
  ## off together to a self-consistent wrong answer.  etaDistInit() solves a
  ## log-normal surrogate -- an ordinary, well-conditioned saem fit -- and
  ## moment-matches its answer back onto the declared families, which stops the
  ## walk before it starts.  Measured on Bauer's four gamma arms, saem at
  ## nBurn=300/nEm=150, mean absolute relative error over CL, V1, both relative
  ## variances and the copula correlation:
  ##
  ##        g1     g2     g3     g4    mean
  ##  cold  10.6%   9.6%   5.1%  25.9%  12.8%
  ##  warm   6.9%   6.4%   5.9%   7.6%   6.7%
  ##
  ## g3 is the one arm it costs anything, and it is the arm that already fits
  ## best -- least to gain, most to disturb.
  ##
  ## No recursion: the surrogate has no dist() declarations left (they become
  ## exp(mu + eta)), so this hook returns NULL for the nested fit.
  if (!identical(control$etaDistWarmStart, FALSE)) {
    .warm <- try(etaDistInit(ui, data), silent=TRUE)
    if (inherits(.warm, "try-error")) {
      warning("the declared-distribution warm start failed; starting values ",
              "are unchanged\n  ",
              conditionMessage(attr(.warm, "condition")), call.=FALSE)
    } else {
      ui <- .warm
      ## the declarations themselves are untouched, but re-read so the stash
      ## below records the WARMED starting values rather than the original ones
      .d <- .rxUiEtaDists(ui)
    }
  }
  ## rxEtaDistExpand() clears the iniDf's `etaDist` column (rxode2
  ## R/etaDist.R), and the copula block it replaces with independent latents
  ## plus rxCor.* thetas -- so after expansion rxUiEtaDists() reports nothing
  ## and the declarations are unrecoverable.  Everything downstream that needs
  ## to know a random effect WAS declared (the ODE-free M-step in particular)
  ## reads this stash instead.
  ## Which representation the estimator asked for.  Read from the control
  ## rather than assumed, and defaulted to "cdf" so a method that has never
  ## heard of the argument (or an older control round-tripped through
  ## do.call) gets exactly the behaviour it had before.
  ## From the CONTROL ARGUMENT first.  The hook is handed `control` directly,
  ## and reading only `rxGetControl(ui, ...)` found nothing -- the value is not
  ## on the ui at this point -- so the route silently stayed "cdf" however it
  ## was set.  Measured: a fit with etaDistParam="direct" completed normally
  ## with no refusal and no expansion change.
  .param <- tryCatch({
    .p <- if (!is.null(control) && !is.null(control$etaDistParam)) {
      control$etaDistParam
    } else {
      rxode2::rxGetControl(ui, "etaDistParam", "cdf")
    }
    if (is.character(.p) && length(.p) >= 1L &&
          .p[1] %in% c("auto", "cdf", "direct")) {
      .p[1]
    } else {
      "auto"
    }
  }, error = function(e) "auto")
  ## "auto" -- the default -- means DIRECT wherever direct can express the
  ## model, and cdf where it cannot.
  ##
  ## The two routes are not equally good at a covariate on a declaration, and
  ## the difference is not close.  A covariate there is a parameter of the
  ## PRIOR, and only the direct route's etas are a genuine sample from that
  ## prior, so only there does an objective identify it.  Measured on a
  ## 120-subject arm, true coefficient 0.75, three seeds, mean |error|:
  ##
  ##   direct (prior/Q2)                     0.044
  ##   cdf, coefficient as a plain theta     0.382
  ##   cdf, coefficient through MCOV         0.237   (and erratic: 0.47-0.98)
  ##
  ## On the cdf route the declared eta is decoded from a latent USING THE
  ## CURRENT parameters, so the density of those decoded etas is stationary
  ## wherever the fit already is -- measured, the objective peaks at the
  ## current value with drift within one standard error at every value tried.
  ## No amount of work on that side recovers the coefficient.
  ##
  ## direct cannot express every model, and says so by name: a correlated block
  ## of more than two declared etas, and a declared eta correlated with an
  ## ordinary one.  Those fall back to cdf rather than erroring, which is why
  ## this is a trial rather than a flag.
  if (identical(.param, "auto")) {
    .estNm <- if (is.character(est) && length(est) == 1L) est else ""
    .param <- if (!grepl("^saem", .estNm)) {
      ## every other estimator refuses "direct" below; do not hand it one
      "cdf"
    } else if (inherits(tryCatch(rxode2::rxEtaDistExpand(ui, param = "direct"),
                                 error = function(e) e), "error")) {
      "cdf"
    } else {
      "direct"
    }
  }
  ## REFUSE "direct" for every estimator that does not consume it.
  ##
  ## The expansion does its half for all of them -- no latent, no decoder, the
  ## eta kept with a FIXED placeholder omega -- but an estimator that does not
  ## read the route reads that placeholder 1 as a Gaussian VARIANCE and
  ## completes normally, having fitted a standard normal random effect where a
  ## gamma was declared.  The wrong model, no error, no warning.  So the list is
  ## of estimators that have been TAUGHT the route, and everything else is
  ## refused by name.
  ##
  ## saem consumes it: `etaDistDirect` reaches saem.cpp, kernel 1 draws from the
  ## declared family (through the Gaussian copula for a correlated pair),
  ## kernels 2 and 3 walk on the bijected scale and pay the log-Jacobian, and
  ## the acceptance scores the declared columns with rxEtaDistLogD /
  ## rxEtaDistPairLogD instead of the quadratic.
  ##
  ## focei and imp do NOT, and there is a mathematical reason beyond "not
  ## written yet" -- see the message below.
  ##
  ## Add an estimator here only WITH a test that fails if the Gaussian prior is
  ## still in use; "it ran and the numbers look plausible" is exactly what the
  ## broken version does.
  .directOk <- c("saem")
  if (identical(.param, "direct")) {
    .est <- if (is.character(est) && length(est) == 1L) est else "this method"
    if (!(.est %in% .directOk)) {
      .why <- if (.est %in% c("focei", "foce", "posthoc", "imp", "impmap")) {
        paste0(
          "  this is not only unimplemented.  The inner MAP adds the prior's ",
          "curvature to the Hessian, and for the direct parameterization that ",
          "is d2 log p/d(eta)2, which is POSITIVE wherever the declared density ",
          "is convex -- every family with an interior-mode-free shape, gamma ",
          "with shape < 1 among them.  Measured on gamma(0.5, 0.5) with a ",
          "residual sd of 0.3, the inner Hessian 1/s^2 + (shape-1)/eta^2 is ",
          "negative for eta < 0.212, which is 35.5% of the prior mass: those ",
          "subjects have no interior mode and so no Laplace expansion at all\n",
          "  the cdf route does not have this problem, because the quantity it ",
          "expands around is the standard normal latent")
      } else {
        paste0(
          "  the declared random effect would be fitted as a standard normal ",
          "-- the wrong model, silently")
      }
      stop("etaDistParam=\"direct\" is not supported for est=\"", .est, "\"\n",
           .why, "\n",
           "  use etaDistParam=\"cdf\" (the default)",
           call. = FALSE)
    }
    ## The direct route USED to require etaDistMstep=TRUE here, and that
    ## refusal named the wrong cause.
    ##
    ## What is true: on this route the declared thetas reach the model only
    ## through the `rxEdA.*` anchors, which nothing reads, so the observation
    ## likelihood does not depend on them -- they parameterize the prior and
    ## nothing else.  What was WRONG was concluding that the family M-step is
    ## therefore mandatory.  The objective that identifies them is the eta
    ## DENSITY, sum_i log p(eta_i | theta), and saem has had that objective all
    ## along (src/saem.cpp, the etaDistLoglik block) -- it was simply trapped
    ## inside `etaDistMstep()`, whose loop is gated on the etaDistMstep control,
    ## so switching the family MLE off took the density objective with it.
    ##
    ## It is now owned by the Q1/Q2 partition instead (.etaDistThetaSplit),
    ## which is NoLimits.jl's rule: a parameter appearing in a random-effect
    ## distribution and in no observation-side expression is identified by the
    ## prior alone and gets its own M-step.  So there is nothing left to refuse.
  }
  .decl <- .etaDistDeclStash(ui, .d, param = .param)
  ## Decompress BEFORE stashing: rxUiDecompress() on a compressed ui returns a
  ## new object, so assigning into it would write to a temporary and the stash
  ## would never reach the ui that is returned.
  .ui <- rxode2::rxUiDecompress(rxode2::rxEtaDistExpand(ui, param = .param))
  ## In `meta`, which is the ONLY container that survives to the estimators.
  ## Measured, on a real saem fit, by planting a probe in each candidate and
  ## seeing which arrived: the ui environment, the control and an extra iniDf
  ## column were all gone by the time the M-step asked, because the ui is
  ## rebuilt and the est method installs its own freshly-built control after
  ## the hooks have run.  `meta` is rxode2's own metadata environment and is
  ## deliberately carried across model rewrites, so it is the one that holds.
  if (!is.null(.decl)) .etaDistDeclSet(.ui, .decl)
  list(ui=.ui)
}

#' Preserve the declarations across `rxEtaDistExpand()`
#'
#' The expansion is lossy by design: it rewrites the model into one with
#' ordinary standard-normal random effects, so the declaration it consumed is
#' no longer anywhere in the ui.  This records the parts that cannot be
#' reconstructed afterwards -- which etas were declared, with what family, and
#' how the copula paired them -- keyed so the expanded model's own thetas and
#' etas can be found from it.
#'
#' Copula pairing is read here rather than later because the block that carries
#' it (a unit-diagonal omega whose off-diagonal IS the correlation) is exactly
#' what the expansion removes.  Only a PAIR is recorded: the M-step drivers
#' reconstruct a single partner, so a larger block returns `NULL` and takes the
#' general path instead of being silently wrong.
#'
#' @param ui the UNEXPANDED rxode2 ui
#' @param d declared random effects, as `rxUiEtaDists()` returns them
#' @return a list, or `NULL` when the block is not one this can describe
#' @noRd
.etaDistDeclStash <- function(ui, d, param = "cdf") {
  .ini <- rxode2::rxUiDecompress(ui)$iniDf
  .n <- nrow(d)
  .netaOf <- function(.nm) {
    .w <- which(.ini$name == .nm & .ini$neta1 == .ini$neta2)
    if (length(.w) == 1L) as.integer(.ini$neta1[.w]) else NA_integer_
  }
  .id <- vapply(d$name, .netaOf, integer(1))
  .cw <- rep(-1L, .n)                       # 0-based partner, or -1
  .ct <- rep(NA_character_, .n)             # the rxCor.* theta carrying it
  .off <- .ini[!is.na(.ini$neta1) & !is.na(.ini$neta2) &
                 .ini$neta1 != .ini$neta2, , drop = FALSE]
  if (nrow(.off) > 0L) {
    for (.r in seq_len(nrow(.off))) {
      .a1 <- which(.id == .off$neta1[.r]); .a2 <- which(.id == .off$neta2[.r])
      if (length(.a1) != 1L || length(.a2) != 1L) next  # not a declared pair
      .hi <- max(.a1, .a2); .lo <- min(.a1, .a2)
      if (.cw[.hi] >= 0L) return(NULL)                  # >2 declared partners
      .cw[.hi] <- .lo - 1L
      ## rxEtaDistExpand() names the Cholesky theta rxCor.<later>.<earlier>
      .ct[.hi] <- paste0("rxCor.", d$name[.hi], ".", d$name[.lo])
    }
  }
  ## `param` rides along HERE, and not on `etaDistInfo` where the expansion
  ## also records it.
  ##
  ## `etaDistInfo` is an environment variable on the ui and it does NOT reach
  ## the estimator: measured on a direct fit, saem.R saw the expanded model
  ## (its eta was `rxd.eta.cl`) while `.etaDistIsDirect()` read FALSE off the
  ## same ui, so the eta-name lookup went looking for `rxz.eta.cl`, found
  ## nothing, and returned no metadata at all -- which the C++ side then read as
  ## "no declarations" and sampled a standard normal.  Marking it sticky is not
  ## enough; saem rebuilds the ui more than once.
  ##
  ## This stash lives in `ui$meta`, which is the one container already proven to
  ## survive that (see the comment at the end of .preProcessEtaDist()), so the
  ## route travels with the declarations it describes rather than separately.
  list(name = as.character(d$name), etaDist = as.character(d$etaDist),
       corWith = .cw, corTheta = .ct, param = as.character(param)[1])
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
  .fix <- try(get("fixef", envir=.env), silent=TRUE)
  if (inherits(.fix, "try-error") || is.null(.fix)) return(NULL)
  .info <- .etaDistInfo(.env)
  .blocks <- if (is.null(.info)) NULL else .info$blocks
  if (is.null(.blocks) || length(.blocks) == 0L) {
    ## `etaDistInfo` does not survive onto the fit's ui, so the blocks are
    ## recovered from the fit itself -- see .etaDistBlocksFromFit()
    .blocks <- .etaDistBlocksFromFit(x[[1]])
  }
  ## No blocks recovered does not mean no correlation: the recovery keys on
  ## `rxz.*` names and `rxCor.*` thetas, and the DIRECT route has neither.
  ## Falling through to NULL here is what made the fallback below unreachable.
  if (length(.blocks) == 0L) return(.etaDistCorFromSampler(.env))
  .cor <- lapply(.blocks, function(.nms) {
    .need <- unlist(lapply(seq_along(.nms), function(.i) {
      if (.i == 1L) return(NULL)
      paste0("rxCor.", .nms[.i], ".", .nms[seq_len(.i - 1L)])
    }), use.names=FALSE)
    if (length(.need) == 0L || !all(.need %in% names(.fix))) return(NULL)
    .etaDistCorFromY(.nms, as.list(.fix[.need]))
  })
  names(.cor) <- vapply(.blocks, function(.n) .n[1], character(1),
                        USE.NAMES=FALSE)
  .cor <- .cor[!vapply(.cor, is.null, logical(1))]
  if (length(.cor) == 0L) return(.etaDistCorFromSampler(.env))
  .cor
}

#' The declared copula correlations the SAMPLER used, for a route with no
#' `rxCor.*` theta
#'
#' `etaDistParam="direct"` keeps the correlation in the omega rather than in a
#' theta, so there is nothing for the `rxCor.*` rebuild above to read and the
#' value the fit USED was previously unavailable to the user at all -- it is
#' estimated (measured moving from a 0.30 start to 0.60 on a pair simulated at
#' 0.60) and then discarded at the end of the fit.
#'
#' saem now returns it, `.getSaemOmega()` stashes it, and this assembles the
#' same per-block correlation matrices the cdf route reports.
#'
#' The value is on the LATENT scale -- it is the Gaussian copula's parameter,
#' not the Pearson correlation of the random effects.  Those differ whenever the
#' marginals are not normal: measured, a latent 0.50 induces an eta-scale
#' correlation of 0.437 for gamma(shape 2) with gamma(shape 0.5), 0.454 with a
#' lognormal, 0.474 for two gammas, and 0.500 only for normal with normal.
#'
#' @param env fit environment
#' @return named list of correlation matrices, or NULL
#' @noRd
#' @author Matthew L. Fidler
.etaDistCorFromSampler <- function(env) {
  .rho <- try(get(".etaDistRhoFit", envir=env), silent=TRUE)
  if (inherits(.rho, "try-error") || is.null(.rho) || length(.rho) == 0L) {
    return(NULL)
  }
  .cw <- try(get(".etaDistCorWithFit", envir=env), silent=TRUE)
  if (inherits(.cw, "try-error") || length(.cw) != length(.rho)) return(NULL)
  .st <- try(.etaDistDeclGet(env$ui), silent=TRUE)
  if (inherits(.st, "try-error") || is.null(.st)) return(NULL)
  .nm <- .st$name
  if (length(.nm) != length(.rho)) return(NULL)
  .out <- list()
  for (.i in seq_along(.rho)) {
    .j <- .cw[.i]
    if (is.na(.j) || .j < 0L) next
    .p <- c(.nm[.j + 1L], .nm[.i])
    .R <- matrix(c(1, .rho[.i], .rho[.i], 1), 2, 2, dimnames=list(.p, .p))
    .out[[.p[1]]] <- .R
  }
  if (length(.out) == 0L) return(NULL)
  .out
}

#' Put the declared copula correlation in `parFixed`
#'
#' On the cdf route the correlation is already there: it is an `rxCor.*` theta
#' with `backTransform("tanh")`, so its back-transformed column IS the
#' correlation and it carries a standard error like any other row.
#'
#' On the direct route there is no such theta, so the correlation appeared
#' nowhere in the printed fit even though it was estimated and used.  This adds
#' it as a row named for the pair it joins, with no standard error -- there is
#' no theta to have one, and inventing a blank column entry is better than
#' implying the value is not an estimate at all.
#'
#' Labelled as the LATENT correlation, because for non-normal marginals it is
#' not the correlation of the random effects: a latent 0.50 induces 0.437 to
#' 0.474 across the family pairings measured, and 0.500 only for normal with
#' normal.
#'
#' @param ret fit object
#' @return `ret`, with the copula rows appended to `parFixed`
#' @noRd
#' @author Matthew L. Fidler
.postFinalEtaDistCorParFixed <- function(ret) {
  .env <- try(ret$env, silent=TRUE)
  if (inherits(.env, "try-error") || is.null(.env)) return(ret)
  .cor <- try(.etaDistCorFromSampler(.env), silent=TRUE)
  if (inherits(.cor, "try-error") || is.null(.cor) || length(.cor) == 0L) {
    return(ret)
  }
  .pfd <- try(get("parFixedDf", envir=.env), silent=TRUE)
  if (inherits(.pfd, "try-error") || is.null(.pfd)) return(ret)
  .rows <- NULL; .nms <- character(0)
  for (.R in .cor) {
    .p <- rownames(.R)
    if (length(.p) != 2L) next
    .nms <- c(.nms, paste0("cor(", .p[1], ",", .p[2], ")"))
    .r <- .pfd[1, , drop=FALSE]
    .r[1, ] <- NA
    if ("Estimate" %in% names(.r)) .r[1, "Estimate"] <- .R[1, 2]
    if ("Back-transformed" %in% names(.r)) {
      .r[1, "Back-transformed"] <- .R[1, 2]
    }
    .rows <- rbind(.rows, .r)
  }
  if (is.null(.rows)) return(ret)
  rownames(.rows) <- .nms
  assign("parFixedDf", rbind(.pfd, .rows), envir=.env)
  .pf <- try(get("parFixed", envir=.env), silent=TRUE)
  if (!inherits(.pf, "try-error") && !is.null(.pf) && ncol(.pf) > 0L) {
    .r2 <- .pf[rep(1, length(.nms)), , drop=FALSE]
    .r2[] <- ""
    .ec <- intersect(c("Estimate", "Back-transformed"), names(.pf))
    for (.k in seq_along(.nms)) {
      for (.c in .ec) .r2[.k, .c] <- format(signif(.rows[.k, "Estimate"], 4))
    }
    rownames(.r2) <- .nms
    assign("parFixed", rbind(.pf, .r2), envir=.env)
  }
  ret
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

#' Drop the latent random effects from the reported parameter table
#'
#' `rxEtaDistExpand()` leaves the latent standard normals (`rxz.<eta>`) in
#' `parFixed`, where they print as `NA` in every column.  Their variance is
#' fixed at one by construction -- that is what makes the copula a copula --
#' so they are not estimates and there is nothing to report for them.
#'
#' The copula correlations (`rxCor.<i>.<j>`) are NOT touched here: the
#' expansion already gives them `backTransform = "tanh"`, so their
#' back-transformed column is the correlation and has always been correct.
#' (An earlier version of this hook recomputed it, on the strength of my
#' having misread the raw Estimate column as a correlation.  It was the
#' reading that was wrong, not the table.)
#'
#' @param ret fit object
#' @return `ret`, without the latent rows
#' @noRd
#' @author Matthew L. Fidler
.postFinalEtaDistParFixed <- function(ret) {
  .env <- try(ret$env, silent=TRUE)
  if (inherits(.env, "try-error") || is.null(.env)) return(ret)
  .pfd <- try(get("parFixedDf", envir=.env), silent=TRUE)
  if (inherits(.pfd, "try-error") || is.null(.pfd)) return(ret)
  ## `rxd.` as well as `rxz.`: the DIRECT route's declared eta is not a latent,
  ## but it is not an estimate either.  Its omega entry is a placeholder the
  ## expansion fixes at 1 because the FAMILY carries the dispersion, so there is
  ## nothing about it to report -- and a printed variance for a gamma-distributed
  ## random effect invites exactly the wrong reading.
  .nm <- rownames(.pfd)
  if (is.null(.nm) || !any(grepl("^rx[zd][.]", .nm))) return(ret)
  assign("parFixedDf", .pfd[!grepl("^rx[zd][.]", .nm), , drop=FALSE], envir=.env)
  .pf <- try(get("parFixed", envir=.env), silent=TRUE)
  if (!inherits(.pf, "try-error") && !is.null(.pf)) {
    .nm2 <- rownames(.pf)
    if (!is.null(.nm2) && any(grepl("^rx[zd][.]", .nm2))) {
      assign("parFixed", .pf[!grepl("^rx[zd][.]", .nm2), , drop=FALSE], envir=.env)
    }
  }
  ret
}

#' Drop a directly-parameterized declared random effect from `$omega`
#'
#' On `etaDistParam="direct"` the declared eta IS the random effect and carries
#' its family as its prior, so it has no variance in the ordinary sense -- the
#' family holds the dispersion.  `rxEtaDistExpand()` fixes its omega entry at a
#' placeholder 1 for exactly that reason: it is machinery, not a parameter.
#'
#' What the fit printed instead was neither.  Measured on a correlated gamma
#' pair, every cell of the 2x2 came back 180.6704 -- so `$omegaR` reported a
#' correlation of 1.000 with an SD of 13.44 -- while the iniDf still carried
#' `fix = TRUE` on those rows.  Everything upstream is correct (the ui omega,
#' `saemModelOmegaFixed` and its values all hold 1.0/0.3), so a reporting path
#' writes into a row it knows is fixed.  That line has NOT been found.
#'
#' This does not paper over that.  Even with the placeholder reported correctly,
#' printing `omega = 1` for a gamma-distributed random effect is the wrong thing
#' to show: a reader takes it for the dispersion, and the dispersion is in the
#' family's own parameters, which are already in `parFixed`.
#'
#' NoLimits.jl has no omega concept at all for this reason -- a declared random
#' effect's distribution parameters are ordinary fixed effects, tagged with a
#' role -- and its governing rule where a family cannot supply a quantity is to
#' OMIT the row with a warning rather than emit a placeholder.  This is that
#' rule.
#'
#' The correlation is a separate loss and is NOT recovered here: on this route
#' it never reaches a theta, so the fit does not contain it (see
#' `.etaDistWarnCorFrozen()`).  Reporting it needs the value plumbed out of the
#' sampler first.
#'
#' @param ret fit object
#' @return `ret`, with the directly-parameterized rows removed from `$omega`
#' @noRd
#' @author Matthew L. Fidler
.postFinalEtaDistDirectOmega <- function(ret) {
  .env <- try(ret$env, silent=TRUE)
  if (inherits(.env, "try-error") || is.null(.env)) return(ret)
  .om <- try(get("omega", envir=.env), silent=TRUE)
  if (inherits(.om, "try-error") || is.null(.om) || !is.matrix(.om)) return(ret)
  .nm <- rownames(.om)
  if (is.null(.nm)) return(ret)
  .drop <- grepl("^rxd[.]", .nm)
  if (!any(.drop)) return(ret)
  assign("omega", .om[!.drop, !.drop, drop=FALSE], envir=.env)
  ret
}

#' Report the declared block in `$omega`, under the names the model used
#'
#' The user writes `eta.cl + eta.v1 ~ c(1, 0.5, 1)` -- a covariance block with
#' the correlation in it.  What comes back is the expansion's internals: a 2x2
#' identity named `rxz.eta.cl` / `rxz.eta.v1`, with the fitted correlation
#' living in a `rxCor.*` theta instead.  The block the user wrote is nowhere in
#' `$omega`.
#'
#' The latent random effects are standard normals, so their covariance matrix
#' IS the correlation matrix -- unit diagonal is what the declaration requires
#' -- and the fitted block goes back into `$omega` on the covariance scale,
#' named as the model named it.  `$omegaR` then derives the correlation view
#' through the machinery every other model uses.
#'
#' The `rxCor.*` rows stay in `parFixed`: that is where their standard error
#' is, on the estimated scale like every other row, and moving the value into
#' `$omega` must not take the uncertainty out of the output.
#'
#' @param ret fit object
#' @return `ret`, with `$omega` carrying the declared block
#' @noRd
#' @author Matthew L. Fidler
.postFinalEtaDistOmega <- function(ret) {
  .env <- try(ret$env, silent=TRUE)
  if (inherits(.env, "try-error") || is.null(.env)) return(ret)
  .om <- try(get("omega", envir=.env), silent=TRUE)
  if (inherits(.om, "try-error") || is.null(.om) || !is.matrix(.om)) return(ret)
  .blocks <- try(.etaDistBlocksFromFit(ret), silent=TRUE)
  if (inherits(.blocks, "try-error") || length(.blocks) == 0L) return(ret)
  .fix <- try(get("fixef", envir=.env), silent=TRUE)
  if (inherits(.fix, "try-error") || is.null(.fix)) return(ret)
  .nm <- rownames(.om)
  if (is.null(.nm)) return(ret)
  for (.b in .blocks) {
    .need <- unlist(lapply(seq_along(.b), function(.i) {
      if (.i == 1L) return(NULL)
      paste0("rxCor.", .b[.i], ".", .b[seq_len(.i - 1L)])
    }), use.names=FALSE)
    if (length(.need) == 0L || !all(.need %in% names(.fix))) next
    .R <- try(.etaDistCorFromY(.b, as.list(.fix[.need])), silent=TRUE)
    if (inherits(.R, "try-error")) next
    .row <- match(paste0("rxz.", .b), .nm)
    if (anyNA(.row)) next
    .om[.row, .row] <- .R
    .nm[.row] <- .b          # report them as the model named them
  }
  dimnames(.om) <- list(.nm, .nm)
  assign("omega", .om, envir=.env)
  ret
}

postFinalObjectHooksAdd(".postFinalEtaDistOmega", .postFinalEtaDistOmega)

postFinalObjectHooksAdd(".postFinalEtaDistParFixed", .postFinalEtaDistParFixed)

postFinalObjectHooksAdd(".postFinalEtaDistDirectOmega",
                        .postFinalEtaDistDirectOmega)

postFinalObjectHooksAdd(".postFinalEtaDistCorParFixed",
                        .postFinalEtaDistCorParFixed)

#' Recover the declared correlation blocks from the fit itself
#'
#' `rxEtaDistExpand()` records what it did in `etaDistInfo` on the ui it
#' returns, but that does not survive to the fit object (measured: present on
#' the expanded ui, absent on `fit$ui`), which is why `fit$etaDistCor` came
#' back NULL for a model that plainly has a declared block.
#'
#' Everything needed is still in the fit, so it is read from there instead of
#' carried: the latent random effects are `rxz.<declared eta>` and appear in
#' the omega in their block order, and the copula parameters are
#' `rxCor.<i>.<j>` thetas naming the pair they connect.  Two random effects
#' are in the same block exactly when such a theta joins them.
#'
#' @param ret fit object
#' @return list of character vectors, one per block, in omega order; empty
#'   when the model declares nothing
#' @noRd
#' @author Matthew L. Fidler
.etaDistBlocksFromFit <- function(ret) {
  .ui <- try(rxode2::rxUiDecompress(ret$ui), silent=TRUE)
  if (inherits(.ui, "try-error") || is.null(.ui)) return(list())
  .ini <- .ui$iniDf
  if (is.null(.ini) || !any(names(.ini) == "neta1")) return(list())
  .e <- .ini[!is.na(.ini$neta1) & .ini$neta1 == .ini$neta2, ]
  if (nrow(.e) == 0L) return(list())
  .e <- .e[order(.e$neta1), ]
  .lat <- .e$name[grepl("^rxz[.]", .e$name)]
  if (length(.lat) == 0L) return(list())
  .dec <- sub("^rxz[.]", "", .lat)
  .th <- .ini$name[!is.na(.ini$ntheta)]
  .cor <- .th[grepl("^rxCor[.]", .th)]
  ## adjacency from the copula thetas; a lone declared random effect is its
  ## own block and simply has no correlation to report
  .grp <- seq_along(.dec)
  for (.c in .cor) {
    .p <- sub("^rxCor[.]", "", .c)
    .i <- which(vapply(.dec, function(.d) startsWith(.p, paste0(.d, ".")),
                       logical(1)))
    for (.ii in .i) {
      .j <- which(.dec == sub(paste0("^", .dec[.ii], "[.]"), "", .p))
      if (length(.j) == 1L) {
        .keep <- min(.grp[.ii], .grp[.j])
        .drop <- max(.grp[.ii], .grp[.j])
        .grp[.grp == .drop] <- .keep
      }
    }
  }
  .out <- lapply(sort(unique(.grp)), function(.g) .dec[.grp == .g])
  .out[vapply(.out, length, integer(1)) > 1L]
}
