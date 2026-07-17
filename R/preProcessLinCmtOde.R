#' Compartments that `linCmt()` owns
#'
#' @noRd
.linCmtOdeStates <- c("depot", "central", "peripheral1", "peripheral2")

#' Estimation methods that add sensitivity compartments to the solved model
#'
#' The FOCEi family (the callers of `.foceiFamilyControl()`) gains
#' `rx__sens_<state>_BY_ETA_#___` states; the nlm family (the callers of
#' `.nlmSetupEnv()`) gains `rx__sens_<state>_BY_THETA_#___` states when it uses
#' analytic derivatives.  Both push the `linCmt()` compartments off the numbers
#' the data was translated against.  The nlm family is listed whatever its
#' `solveType` is: only a derivative-free solve escapes the shift, and the
#' translation is equivalent (just solved numerically) in that case anyway.
#' Keep in sync with the callers of those two functions.
#'
#' @noRd
.linCmtOdeEstFamily <- c(
  # FOCEi family -- .foceiFamilyControl()
  "focei", "foce", "focep", "fo", "foi", "laplace", "agq", "posthoc",
  "impmap", "imp", "qrpem", "npb", "npag", "advi",
  "mfocei", "mfoce", "mfocep", "mlaplace", "magq", "mnpb", "mnpag",
  "ifocei", "ifoce", "ifocep", "ilaplace", "iagq", "inpb", "inpag",
  "foceif", "focef", "focepf", "mfoceif", "mfocef", "mfocepf",
  "ifoceif", "ifocef", "ifocepf",
  # nlm family -- .nlmSetupEnv()
  "nlm", "nlminb", "bobyqa", "newuoa", "uobyqa", "n1qn1", "lbfgsb3c",
  "optim", "nls"
)

#' Does this model mix `linCmt()` with ODE states?
#'
#' `ui$mvL` is only present for a `linCmt()` model, so a hand-written
#' depot/central ODE model is not caught here.
#'
#' @param ui rxode2 ui
#' @return logical, `TRUE` when `linCmt()` coexists with at least one ODE state
#' @noRd
#' @author Matthew L. Fidler
.uiIsMixedLinCmtOde <- function(ui) {
  if (is.null(ui$mvL)) return(FALSE)
  length(setdiff(ui$state, .linCmtOdeStates)) > 0L
}

#' Reorder the `d/dt()` lines of a linToOde() model to a target state order
#'
#' rxode2 numbers compartments by the order the `d/dt()` lines appear, and
#' `linToOde()` emits the linear compartments where `linCmt()` was called.  In
#' the original model the `linCmt()` compartments are appended last instead, so
#' without this the translation would silently renumber the compartments a
#' numeric `cmt` in the data refers to.
#'
#' @param ui the `linToOde()` translated ui
#' @param state the state order to restore (the original `linCmt()` model's)
#' @return ui with the `d/dt()` lines reordered, or `ui` when already in order
#' @noRd
#' @author Matthew L. Fidler
.linCmtOdeRestoreStateOrder <- function(ui, state) {
  if (identical(ui$state, state)) return(ui)
  ui <- rxode2::rxUiDecompress(ui)
  .lst <- ui$lstExpr
  .isDdt <- vapply(.lst, function(e) {
    is.call(e) && length(e) >= 2L && is.call(e[[2]]) &&
      identical(e[[2]][[1]], quote(`/`)) &&
      identical(e[[2]][[2]], quote(d))
  }, logical(1))
  if (!any(.isDdt)) return(ui)
  .ddtState <- vapply(.lst[.isDdt], function(e) as.character(e[[2]][[3]][[2]]),
                      character(1))
  if (!setequal(.ddtState, state)) return(ui)
  # Gather the d/dt() lines, in the target order, at the last d/dt() position.
  # They cannot simply be permuted among the slots they already occupy: the
  # `linCmt()` output assignment (e.g. C2 <- central/v) sits between them, and a
  # d/dt() moved ahead of it would read C2 before it is defined.  Every d/dt()
  # RHS only needs the assignments that already preceded the last d/dt().
  .idx <- which(.isDdt)
  .at <- max(.idx)
  .ddt <- .lst[.idx][order(match(.ddtState, state))]
  .lst <- append(.lst[-.idx], .ddt, after=sum(!.isDdt[seq_len(.at)]))
  .rebuildRxUiFromLstExpr(ui, .lst)
}

#' Rebuild an rxUi from a modified lstExpr
#'
#' Mirrors rxode2's `linToOde()`/`.rebuildRxUiFromExpr()` reconstruction.
#'
#' @param ui the (decompressed) ui to take `meta`/`iniFun` from
#' @param expr the replacement `lstExpr`
#' @return the rebuilt ui
#' @noRd
#' @author Matthew L. Fidler
.rebuildRxUiFromLstExpr <- function(ui, expr) {
  .ls <- ls(ui$meta, all.names=TRUE)
  .hasIni <- length(ui$iniDf$cond) > 0L
  .ret <- vector("list", length(.ls) + if (.hasIni) 3L else 2L)
  .ret[[1L]] <- quote(`{`)
  for (.i in seq_along(.ls)) {
    .ret[[.i + 1L]] <- rxode2::rxUiDeparse(ui$meta[[.ls[.i]]], .ls[.i])
  }
  .len <- length(.ls)
  if (.hasIni) {
    .ret[[.len + 2L]] <- ui$iniFun
    .ret[[.len + 3L]] <- bquote(model(.(as.call(c(quote(`{`), expr)))))
  } else {
    .ret[[.len + 2L]] <- bquote(model(.(as.call(c(quote(`{`), expr)))))
  }
  .fun <- function() {}
  body(.fun) <- as.call(.ret)
  if (is.function(ui$model)) environment(.fun) <- environment(ui$model)
  suppressMessages(rxode2::as.rxUi(.fun))
}

#' Translate a mixed `linCmt()`/ODE model to all-ODEs
#'
#' rxode2 requires the `linCmt()` compartments to be the last states of the
#' solve (`op$linOffset = neq - numLin - numLinSens`), so their compartment
#' number is one past the ODE states.  The FOCEi inner model adds an ODE state
#' per eta, and the nlm family one per theta, which pushes depot/central past
#' the compartment numbers the data was translated against
#' (`.foceiPreProcessData()` uses the plain model).  A dose then silently lands
#' in a sensitivity state -- every prediction comes back 0.  Solving the linear
#' part as ODEs removes the `linCmt()` block, so the numbering agrees again.
#' SAEM and NLME build no sensitivity states (NLME solves a pred-only
#' finite-difference model) and are left alone.
#'
#' @param ui rxode2 ui
#' @inheritParams nlmixr2
#' @return list with the translated ui, or `NULL` when nothing to do
#' @export
#' @author Matthew L. Fidler
.preProcessLinCmtOde <- function(ui, est, data, control) {
  if (!(est %in% .linCmtOdeEstFamily)) return(NULL)
  if (!.uiIsMixedLinCmtOde(ui)) return(NULL)
  .state <- ui$state
  .ui <- try(rxode2::linToOde(ui), silent = TRUE)
  if (inherits(.ui, "try-error")) return(NULL)
  # The model no longer mixes solved and ODE compartments, which is what was
  # asked for; say so rather than quietly changing how the model is solved.
  warning("'", est, "' cannot use the analytic 'linCmt()' in a model that also has ODEs ",
          "(the sensitivity compartments it adds renumber the linear compartments); ",
          "the linear compartments are solved as ODEs instead, which is slower. ",
          "'est=\"saem\"' uses the analytic 'linCmt()'.",
          call.=FALSE)
  .ui <- .linCmtOdeRestoreStateOrder(.ui, .state)
  if (!identical(.ui$state, .state)) {
    # never renumber the data's compartments quietly
    warning("mixed 'linCmt()'/ODE model: compartments were renumbered from '",
            paste(.state, collapse="', '"), "' to '",
            paste(.ui$state, collapse="', '"),
            "'; refer to compartments by name in 'cmt' rather than by number",
            call.=FALSE)
  }
  list(ui = .ui)
}

preProcessHooksAdd(".preProcessLinCmtOde", .preProcessLinCmtOde)
