#'@rdname nlmixr2Est
#'@export
nlmixr2Est.rxSolve <- function(env, ...) {
  .ui <- get("ui", envir=env)
  .rxControl <- get("control", envir=env)
  if (is.null(.rxControl)) {
    .rxControl <- rxode2::rxControl()
  }
  .events <- get("data", envir=env)
  do.call(rxode2::rxSolve, c(list(object = .ui, params = NULL,
                                  events = .events, inits = NULL), .rxControl,
                             list(theta = NULL, eta = NULL)))
}
