# Expand independent ETA distributions before the SAEM model and mu-reference
# hooks run. The sampler operates on fixed unit-normal latents, while the
# observation model decodes them through each declared marginal's quantile.
.preProcessSaemEtaDist <- function(ui, est, data, control) {
  if (is.null(ui)) return(NULL)
  .rx <- asNamespace("rxode2")
  .decl <- get0("rxUiEtaDists", envir = .rx, inherits = FALSE)
  if (is.null(.decl)) return(NULL)
  .distributions <- .decl(ui)
  if (nrow(.distributions) == 0L) return(NULL)

  if (!identical(est, "saem")) {
    stop("Declared ETA distributions currently require est = 'saem' on this branch.",
         call. = FALSE)
  }
  .omega <- ui$omega
  .declared <- match(.distributions$name, rownames(.omega))
  if (anyNA(.declared)) {
    stop("Could not map the declared ETA distributions to omega.", call. = FALSE)
  }
  diag(.omega) <- 0
  if (any(.omega[.declared, , drop = FALSE] != 0)) {
    stop("Declared ETA distributions must be independent.", call. = FALSE)
  }
  .expand <- get0("rxEtaDistExpand", envir = .rx, inherits = FALSE)
  if (is.null(.expand)) {
    stop("This rxode2 build cannot expand declared ETA distributions.", call. = FALSE)
  }

  # Keep the default SAEM controls for ordinary models. Declared marginal
  # parameters enter the observation model through the quantile transform.
  if (is.null(control)) control <- saemControl()
  control$nonMuThetaOpt <- "bhhh"
  list(ui = .expand(ui, param = "cdf"), control = control)
}

preProcessHooksAdd(".preProcessSaemEtaDist", .preProcessSaemEtaDist)
