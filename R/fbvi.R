# est="fbvi" -- full-Bayes variational inference.  Shares the whole engine with
# est="emvi" (src/inner.cpp adviOptimize_); the only difference is that the
# variational posterior also covers the unconstrained population vector under
# flat priors instead of point-estimating it by an M-step.  Selected by
# pointEstimate = FALSE, the same way est="imp" is impmap with mapIter = 0.

#' Control for the fbvi (full-Bayes variational inference) method
#'
#' A convenience wrapper around [emviControl()] with `pointEstimate = FALSE`,
#' i.e. the variational posterior covers the unconstrained population vector
#' (under flat priors) as well as the per-subject etas, rather than
#' point-estimating the population parameters by an M-step.  See [emviControl()]
#' for the full parameter list.
#'
#' @inheritParams emviControl
#' @param ... Parameters passed to [emviControl()].
#' @return An `emviControl` object with `pointEstimate = FALSE`.
#' @export
#' @author Matthew L. Fidler
#' @examples
#'
#' fbviControl()
fbviControl <- function(..., pointEstimate = FALSE) {
  ## Construct then mutate, exactly as impControl() does over impmapControl().
  ## Passing `pointEstimate = ` INTO emviControl() instead would bind formal #5
  ## by name and drop it out of positional matching, so a caller supplying five
  ## or more positional arguments would have every later one shifted a slot
  ## (fbviControl(42, 300, 1, "fullRank", FALSE, "adam") landing "adam" on
  ## adaptEta).  Mutating afterwards leaves emviControl()'s own matching intact.
  checkmate::assertLogical(pointEstimate, len = 1, any.missing = FALSE)
  .ctl <- emviControl(...)
  .ctl$pointEstimate <- pointEstimate
  .ctl
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.fbvi <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) return(fbviControl())
  # everything else -- a bare list, an emviControl(), an fbviControl() -- goes
  # through the shared validation, which normalizes it and then lets `est` force
  # pointEstimate = FALSE, announcing the override if it contradicts
  .viValidCtl(control, FALSE, "fbvi")
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.fbvi <- function(x, ...) .viGetControl(x)

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.fbvi <- function(env, ...) .viEst(env, FALSE, "fbvi")
attr(nlmixr2Est.fbvi, "covPresent") <- TRUE
## optimization runs in the unconstrained real coordinate space
attr(nlmixr2Est.fbvi, "unbounded") <- TRUE
attr(nlmixr2Est.fbvi, "iov") <- TRUE
