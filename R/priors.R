## Prior distributions declared in the `ini({})` block.
##
## A prior that an estimation method cannot use must not be silently
## ignored -- the fit would then do something other than what the model
## says, with nothing to tell the user.  Every method therefore either
## uses the priors it is given or refuses them, and that is enforced in
## one place (`nlmixr2Est()`) rather than method by method, so a method
## added later cannot forget to do it.
##
## A method says what it supports with an attribute on itself:
##
##   attr(nlmixr2Est.myMethod, "nlmixr2Priors") <- "normal"
##
## The levels are
##
## - `"none"` (also the default, when the attribute is absent) -- the
##   method cannot use priors at all
## - `"normal"` -- normal priors only, ie `dnorm()`, `stdNormal()` and
##   the `multiNormal()` family the lotri shorthand produces
## - `"nwpri"` -- normal priors, and degrees of freedom on an omega
##   block (`invWishart(4)`), the way a NONMEM NWPRI model works; a
##   normal prior on the omega values themselves (TNPRI) is refused
## - `"all"` -- everything, so nothing is checked here
##
## The list can grow as methods gain support.

#' Levels a method may declare for `nlmixr2Priors`
#'
#' @noRd
.nlmixr2PriorLevels <- c("none", "normal", "nwpri", "all")

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

#' Refuse the priors the dispatched estimation method cannot use
#'
#' @param env nlmixr2 estimation environment
#' @return nothing, called for the error
#' @noRd
#' @author Matthew L. Fidler
.nlmixr2AssertPriors <- function(env) {
  .support <- .nlmixr2PriorSupport(env)
  if (.support == "all") return(invisible())
  .ui <- get("ui", envir=env)
  .extra <- paste0(" with est=\"", class(env)[1], "\"")
  if (.support == "none") {
    .f <- .nlmixr2RxAssert("assertRxUiNoPriors")
    if (is.null(.f)) return(invisible(.nlmixr2AssertNoPriorsFallback(.ui, extra=.extra)))
    .f(.ui, extra=.extra)
  } else if (.support == "normal") {
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
