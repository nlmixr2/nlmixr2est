## Prior distributions declared in the `ini({})` block.
##
## A prior that an estimation method cannot use must not be silently
## ignored -- the fit would then do something other than what the model
## says, with nothing to tell the user.  Every method therefore either
## uses the priors it is given or refuses them, and that is enforced in
## one place (`nlmixr2Est()`) rather than method by method, so a method
## added later cannot forget to do it.
##
## The shared kernel this gate protects (`rxode2::rxPriorBuildSpec()` /
## `rxPriorLogDensity()`, nlmixr2/rxode2#1270) implements THREE genuinely
## different omega conventions, matching NONMEM/Monolix/textbook-Bayes
## respectively -- `"general"`, `"nwpri"`, `"tnpri"` -- see
## `.nlmixr2PriorMethod()` below for how a model's own `ini({})` syntax
## picks among them (never a user- or method-supplied default: the same
## `invWishart(nu)` syntax means a different number under "general" and
## "nwpri", so guessing wrong would silently fit the wrong penalty).
##
## A method says what it supports with an attribute on itself:
##
##   attr(nlmixr2Est.myMethod, "nlmixr2Priors") <- "general"
##
## The levels are
##
## - `"none"` (also the default, when the attribute is absent) -- the
##   method cannot use priors at all
## - `"theta"` -- priors on population parameters only (`dnorm()`,
##   `dcauchy()`, `stdNormal()`, the `multiNormal()` family the lotri
##   shorthand produces among thetas); anything that touches an omega
##   element is refused.  No FOCEi-family method uses this level (see
##   `"general"` below); it exists for a future method whose omega
##   optimization cannot yet honour a prior on it at all.
## - `"general"` -- everything the kernel's `"general"` method covers:
##   the above, plus a normal prior directly on an omega element and a
##   textbook inverse-Wishart on an omega block.  No further check here;
##   `rxPriorBuildSpec()` itself is the validator (nothing beyond what
##   `rxUiPriors()` reports can reach it).  This is what FOCEi's family
##   declares (#931): the natural-scale omega gradient is chain-ruled into
##   `op_focei.cholOmegaInv`'s own estimation-scale parameterization via
##   the SAME derivative data (`d.omegaInv`/`tr.28`, from the model's
##   `_rxInv` handle) FOCEi's own (non-prior) omega gradient already
##   relies on (`foceiPriorOmegaGradAdd()`, `src/inner.cpp`).
## - `"nwpri"` -- normal priors, and degrees of freedom on an omega
##   block (`invWishart(4)`), evaluated with NONMEM's own `$PRIOR NWPRI`
##   convention rather than the textbook one; a normal prior on the
##   omega values themselves (TNPRI) is refused
## - `"tnpri"` -- normal priors including directly on omega elements
##   (Monolix's/NONMEM's own-estimation joint-normal assumption);
##   `dcauchy()` and `invWishart()` are refused (that is `"nwpri"`'s
##   mechanism)
## - `"all"` -- everything, so nothing is checked here; for pseudo-methods
##   (`"output"`, `"posthoc"`) that evaluate an already-accepted model
##   rather than estimate one, so there is no prior they could ignore
##
## The list can grow as methods gain support.

#' Levels a method may declare for `nlmixr2Priors`
#'
#' @noRd
.nlmixr2PriorLevels <- c("none", "theta", "general", "nwpri", "tnpri", "all")

#' What priors does the method dispatched for this environment support?
#'
#' @param env nlmixr2 estimation environment
#' @return one of `.nlmixr2PriorLevels`
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2PriorSupport <- function(env) {
  for (.cls in class(env)) {
    .fn <- utils::getS3method("nlmixr2Est", .cls, optional=TRUE)
    if (is.null(.fn)) next
    .a <- attr(.fn, "nlmixr2Priors")
    if (is.null(.a)) return("none")
    if (length(.a) != 1L || !(.a %in% .nlmixr2PriorLevels)) {
      stop("the 'nlmixr2Priors' of the '", .cls, "' estimation method must be one of '",
           paste(.nlmixr2PriorLevels, collapse="', '"), "'", call.=FALSE)
    }
    return(.a)
  }
  "none"
}

#' Does this rxode2 supply the prior assertions?
#'
#' They arrived with the `prior` column; with an older rxode2 the
#' equivalent check is done here instead so that a prior is still never
#' quietly dropped.
#'
#' @param what name of the assertion
#' @return the function, or NULL
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2RxAssert <- function(what) {
  .ns <- asNamespace("rxode2")
  if (!exists(what, envir=.ns, inherits=FALSE)) return(NULL)
  get(what, envir=.ns)
}

#' Fallback used when rxode2 predates the prior assertions
#'
#' @param ui rxode2 ui
#' @param extra text appended to the error
#' @return nothing, called for the error
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2AssertNoPriorsFallback <- function(ui, extra="") {
  .iniDf <- ui$iniDf
  if (is.null(.iniDf) || !any(names(.iniDf) == "prior")) return(invisible(ui))
  .w <- which(!is.na(.iniDf$prior))
  if (length(.w) > 0L) {
    stop("the model specifies prior distribution(s) on ",
         paste0("'", .iniDf$name[.w], "'", collapse=", "),
         ", which this estimation method cannot use", extra,
         call.=FALSE)
  }
  invisible(ui)
}

#' Evaluate `expr` with the prior gate bypassed
#'
#' For internal re-entries that evaluate an already-accepted model at fixed
#' parameters: the zero-iteration focei runs behind `setOfv()`, `addCwres()`
#' and the impmap objective recompute all re-dispatch with `est="focei"`,
#' which declares no prior support.  By then the priors were already accepted
#' (or refused) by the estimation method that produced the fit, and the
#' fixed-parameter evaluation does not use them, so refusing again would only
#' break post-processing of a prior-carrying fit (#938).
#'
#' Scoped: the flag is restored on exit, so a user-initiated estimation is
#' never affected.  The field is deliberately absent from
#' `.nlmixr2globalReset()` -- the nested `nlmixr2()` call resets
#' `nlmixr2global`, and the flag has to survive it (reset assigns known
#' fields; it does not clear the environment).
#'
#' @param expr expression to evaluate
#' @return the value of `expr`
#' @noRd
.nlmixr2PriorGateBypass <- function(expr) {
  .saved <- nlmixr2global$nlmixr2PriorGateBypass
  nlmixr2global$nlmixr2PriorGateBypass <- TRUE
  on.exit(assign("nlmixr2PriorGateBypass", .saved, envir=nlmixr2global))
  force(expr)
}

#' Refuse a prior that references an omega element
#'
#' Used for the `"theta"` level, for a method that can honour a prior on a
#' population parameter but not yet one on an omega element -- FOCEi's
#' family declares `"general"` instead (#931), since it does wire that
#' natural-scale gradient through `op_focei.cholOmegaInv`
#' (`foceiPriorOmegaGradAdd()`, `src/inner.cpp`); no FOCEi-family method
#' uses this level.
#'
#' @param ui rxode2 ui
#' @param extra text appended to the error
#' @return the ui, invisibly; called for the error otherwise
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2AssertThetaOnlyPriors <- function(ui, extra="") {
  .p <- rxode2::rxUiPriors(ui)
  if (length(.p$name) == 0L) return(invisible(ui))
  .bad <- which(!is.na(.p$neta1))
  if (length(.bad) > 0L) {
    stop("the model puts a prior on the omega parameter(s) ",
         paste0("'", .p$name[.bad], "'", collapse=", "),
         ", which this estimation method cannot use yet", extra,
         call.=FALSE)
  }
  invisible(ui)
}

#' Refuse the priors the dispatched estimation method cannot use
#'
#' @param env nlmixr2 estimation environment
#' @return nothing, called for the error
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2AssertPriors <- function(env) {
  if (isTRUE(nlmixr2global$nlmixr2PriorGateBypass)) return(invisible())
  .support <- .nlmixr2PriorSupport(env)
  if (.support == "all") return(invisible())
  .ui <- get("ui", envir=env)
  .extra <- paste0(" with est=\"", class(env)[1], "\"")
  if (.support == "none") {
    .f <- .nlmixr2RxAssert("assertRxUiNoPriors")
    if (is.null(.f)) return(invisible(.nlmixr2AssertNoPriorsFallback(.ui, extra=.extra)))
    .f(.ui, extra=.extra)
  } else if (.support == "theta") {
    .nlmixr2AssertThetaOnlyPriors(.ui, extra=.extra)
  } else if (.support == "general") {
    ## rxPriorBuildSpec() itself is the validator; nothing beyond what
    ## rxUiPriors() reports can reach it.
    return(invisible())
  } else if (.support == "tnpri") {
    ## normal priors, including directly on omega, are fine; dcauchy()
    ## and invWishart() (that's "nwpri"'s mechanism) are not
    .f <- .nlmixr2RxAssert("assertRxUiNormalPriors")
    if (is.null(.f)) return(invisible())
    .f(.ui, extra=.extra)
  } else if (.support == "nwpri") {
    ## normal priors and omega degrees of freedom are fine; a normal
    ## prior on the omega values themselves is not
    .f <- .nlmixr2RxAssert("assertRxUiNoOmegaNormalPriors")
    if (is.null(.f)) return(invisible())
    .f(.ui, extra=.extra)
  }
  invisible()
}

#' Which kernel method does this model's own prior syntax mean?
#'
#' `rxPriorBuildSpec(ui, method=)` implements three genuinely different
#' omega conventions that are NOT interchangeable -- the same
#' `invWishart(nu)` syntax means a different number under `"general"`'s
#' textbook Wishart than under `"nwpri"`'s NONMEM degrees-of-freedom
#' convention.  So there is no safe default to fall back on; the method
#' has to be read off of what the model's own `ini({})` actually wrote:
#'
#' - a normal prior directly on an omega element (`om.eta ~ 0.01`) only
#'   makes sense as `"tnpri"` (Monolix's/NONMEM's own-estimation
#'   assumption) -- `"nwpri"` refuses it outright
#' - `invWishart()`/`wishart()` degrees of freedom on an omega block only
#'   makes sense as `"nwpri"` (NONMEM's own `$PRIOR NWPRI`) -- `"tnpri"`
#'   refuses it outright
#' - a model mixing both conventions, or using `dcauchy()` anywhere
#'   (neither `"nwpri"` nor `"tnpri"` has a Cauchy analogue), can only be
#'   `"general"` -- the one method that accepts everything `rxUiPriors()`
#'   can report, at the cost of not being either specialized convention
#' - a model with no omega-referencing prior and no Cauchy: the three
#'   methods agree exactly (they differ only in omega/Cauchy handling),
#'   so `"general"` is returned as a harmless default
#'
#' @param ui rxode2 ui
#' @return one of `"general"`, `"nwpri"`, `"tnpri"`
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2PriorMethod <- function(ui) {
  .hasWishart <- FALSE
  .hasOmegaNormal <- FALSE
  .hasCauchy <- FALSE
  .testOmegaDf <- .nlmixr2RxAssert("testRxUiOmegaDf")
  if (!is.null(.testOmegaDf)) .hasWishart <- isTRUE(tryCatch(.testOmegaDf(ui), error=function(e) FALSE))
  .testOmegaNormal <- .nlmixr2RxAssert("testRxUiOmegaNormalPriors")
  if (!is.null(.testOmegaNormal)) {
    .hasOmegaNormal <- isTRUE(tryCatch(.testOmegaNormal(ui), error=function(e) FALSE))
  }
  .stanName <- .nlmixr2RxAssert(".rxPriorStanName")
  if (!is.null(.stanName)) {
    .p <- rxode2::rxUiPriors(ui)
    if (length(.p$name) > 0L) {
      .hasCauchy <- any(vapply(.p$prior, function(p) {
        .fn <- try(str2lang(p)[[1]], silent=TRUE)
        if (inherits(.fn, "try-error")) return(FALSE)
        .fn <- as.character(.fn)
        if (length(.fn) != 1L) return(FALSE)
        .s <- try(.stanName(.fn), silent=TRUE)
        !inherits(.s, "try-error") && !is.na(.s) && identical(.s, "cauchy")
      }, logical(1), USE.NAMES=FALSE))
    }
  }
  if (.hasWishart && .hasOmegaNormal) return("general")
  if (.hasCauchy) return("general")
  if (.hasWishart) return("nwpri")
  if (.hasOmegaNormal) return("tnpri")
  "general"
}

#' Build this model's prior spec for the shared FOCEI-family C++ kernel
#'
#' A thin no-op (`NULL`) for a model with no priors, or when the installed
#' rxode2 predates `rxPriorBuildSpec()` (nlmixr2/rxode2#1270) -- so a
#' caller can invoke this unconditionally rather than needing its own
#' "does this model have a prior" branch first.  The gate
#' (`.nlmixr2AssertPriors()`) has already run before dispatch, so every
#' term reaching this build is one the calling method accepted -- but
#' accepting a level (e.g. `"general"`) does not mean every kernel method
#' can represent it; `rxPriorBuildSpec()` itself is what actually errors
#' (before any estimation starts) for a term the requested `method` cannot
#' express, e.g. `method="tnpri"` on a model with an `invWishart()` prior.
#'
#' @param ui rxode2 ui
#' @param method one of `"auto"` (default -- see `.nlmixr2PriorMethod()`
#'   for how the omega convention is read off the model), `"general"`,
#'   `"nwpri"`, `"tnpri"`.  Typically `foceiControl(priorMethod=)`, passed
#'   straight through by the caller.
#' @return an R external pointer (`rx_prior_spec_t*`) for
#'   `foceiControl(priorSpec=)` to carry into `op_focei`, or `NULL`
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2BuildPriorSpec <- function(ui, method=c("auto", "general", "nwpri", "tnpri")) {
  method <- match.arg(method)
  if (length(rxode2::rxUiPriors(ui)$name) == 0L) return(NULL)
  .build <- .nlmixr2RxAssert("rxPriorBuildSpec")
  if (is.null(.build)) return(NULL)
  if (identical(method, "auto")) method <- .nlmixr2PriorMethod(ui)
  .build(ui, method=method)
}
