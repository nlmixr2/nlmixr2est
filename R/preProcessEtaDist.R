# Expand ETA distributions before the model and mu-reference hooks run.
# rxode2 uses fixed unit-normal latents and a Gaussian copula when a declared
# block has multiple correlated ETAs.
.preProcessEtaDist <- function(ui, est, data, control) {
  if (is.null(ui)) return(NULL)
  .rx <- asNamespace("rxode2")
  .decl <- get0("rxUiEtaDists", envir = .rx, inherits = FALSE)
  if (is.null(.decl)) return(NULL)
  .distributions <- .decl(ui)
  if (nrow(.distributions) == 0L) return(NULL)

  if (!isTRUE(est %in% c("saem", "focei", "imp", "impmap", "vae"))) {
    stop("Declared ETA distributions require SAEM, FOCEI, IMP, or VAE.",
         call. = FALSE)
  }
  .omega <- ui$omega
  .declared <- match(.distributions$name, rownames(.omega))
  if (anyNA(.declared)) {
    stop("Could not map the declared ETA distributions to omega.", call. = FALSE)
  }
  .undeclared <- setdiff(seq_len(nrow(.omega)), .declared)
  if (any(.omega[.declared, .undeclared, drop = FALSE] != 0)) {
    stop("Declare a marginal for every ETA in a correlated non-normal block.",
         call. = FALSE)
  }
  .expand <- get0("rxEtaDistExpand", envir = .rx, inherits = FALSE)
  if (is.null(.expand)) {
    stop("This rxode2 build cannot expand declared ETA distributions.", call. = FALSE)
  }

  expanded <- .expand(ui, param = "cdf")
  if (!identical(est, "saem")) return(list(ui = expanded))

  # Declared marginal parameters enter SAEM's observation model through the
  # quantile transform. Use BHHH to update these non-mu parameters.
  if (is.null(control)) control <- saemControl()
  control$nonMuThetaOpt <- "bhhh"
  list(ui = expanded, control = control)
}

preProcessHooksAdd(".preProcessEtaDist", .preProcessEtaDist)
