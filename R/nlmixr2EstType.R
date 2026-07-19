# Tagged estimation-method list used by nlmixr2Est.default() to print a helpful,
# category-grouped list of the available est= methods when an unsupported method
# is requested (issue #750).  Each user-facing method carries a `type` (the
# mathematical category from the NLME overview) and a `description`.  Third-party
# packages can join the printed list by setting attr(nlmixr2Est.<method>, "type")
# and attr(nlmixr2Est.<method>, "description") on their own S3 method, matching
# the existing "iov"/"covPresent"/"mu" attribute convention.

# Canonical display order of the categories; unknown categories sort after these.
.nlmixr2EstTypeOrder <- c(
  "Linearized",
  "Integral approximation",
  "Stochastic EM",
  "Nonparametric",
  "Machine learning",
  "Model Based Meta Analysis",
  "Optimal Design",
  "Optimizer (NLM family)",
  "External"
)

# est -> list(type=, description=) for the built-in, user-facing methods.
.nlmixr2EstTypeInfo <- list(
  fo      = list(type="Linearized", description="First-Order"),
  foi     = list(type="Linearized", description="First-Order with Interaction"),
  foce    = list(type="Linearized", description="First-Order Conditional Estimation"),
  focei   = list(type="Linearized", description="FOCE with Interaction"),
  focep   = list(type="Linearized", description="FOCE+ (residual at conditional eta)"),
  nlme    = list(type="Linearized", description="Lindstrom-Bates alternating (nlme)"),
  laplace = list(type="Integral approximation", description="Laplace approximation"),
  agq     = list(type="Integral approximation", description="Adaptive Gaussian Quadrature"),
  imp     = list(type="Integral approximation", description="Importance sampling (no MAP search)"),
  impmap  = list(type="Integral approximation", description="Importance sampling (MAP)"),
  saem    = list(type="Stochastic EM", description="Stochastic Approximation EM"),
  qrpem   = list(type="Stochastic EM", description="Quasi-Random Parametric EM"),
  npag    = list(type="Nonparametric", description="NonParametric Adaptive Grid"),
  npb     = list(type="Nonparametric", description="Nonparametric Bayes"),
  advi    = list(type="Machine learning", description="Automatic Differentiation Variational Inference"),
  vae     = list(type="Machine learning", description="Variational autoencoder NLME"),
  nlm      = list(type="Optimizer (NLM family)", description="nlm quasi-Newton"),
  nlminb   = list(type="Optimizer (NLM family)", description="PORT nlminb"),
  bobyqa   = list(type="Optimizer (NLM family)", description="BOBYQA (derivative-free)"),
  newuoa   = list(type="Optimizer (NLM family)", description="NEWUOA (derivative-free)"),
  uobyqa   = list(type="Optimizer (NLM family)", description="UOBYQA (derivative-free)"),
  n1qn1    = list(type="Optimizer (NLM family)", description="n1qn1 (BFGS)"),
  lbfgsb3c = list(type="Optimizer (NLM family)", description="L-BFGS-B"),
  optim    = list(type="Optimizer (NLM family)", description="Nelder-Mead / BFGS (optim)"),
  nls      = list(type="Optimizer (NLM family)", description="nonlinear least squares")
)

#' Stamp the type/description attributes onto the built-in nlmixr2Est methods
#'
#' Runs from `.onLoad()` so that `attr(nlmixr2Est.focei, "type")` and
#' `attr(nlmixr2Est.focei, "description")` are populated for the built-in
#' methods; called for side effects.
#' @return nothing
#' @noRd
.nlmixr2EstTypeApply <- function() {
  .ns <- asNamespace("nlmixr2est")
  for (.est in names(.nlmixr2EstTypeInfo)) {
    .nm <- paste0("nlmixr2Est.", .est)
    if (!exists(.nm, envir=.ns, inherits=FALSE)) next
    .fn <- get(.nm, envir=.ns, inherits=FALSE)
    if (!is.function(.fn)) next
    .info <- .nlmixr2EstTypeInfo[[.est]]
    attr(.fn, "type") <- .info$type
    attr(.fn, "description") <- .info$description
    .locked <- bindingIsLocked(.nm, .ns)
    if (.locked) unlockBinding(.nm, .ns)
    assign(.nm, .fn, envir=.ns)
    if (.locked) lockBinding(.nm, .ns)
  }
  invisible()
}

#' Collect the tagged estimation methods grouped by category
#'
#' Reads the central registry first, then falls back to the `type`/`description`
#' attributes on each S3 method so third-party methods can join the list.
#' @return list of `list(est=, type=, description=)`, one per tagged method
#' @noRd
.nlmixr2EstTypeTagged <- function() {
  .all <- nlmixr2AllEst()
  .rows <- lapply(.all, function(.est) {
    .info <- .nlmixr2EstTypeInfo[[.est]]
    if (is.null(.info)) {
      .fn <- try(utils::getS3method("nlmixr2Est", .est), silent=TRUE)
      if (!inherits(.fn, "try-error")) {
        .type <- attr(.fn, "type")
        if (!is.null(.type)) {
          .desc <- attr(.fn, "description")
          .info <- list(type=.type, description=if (is.null(.desc)) "" else .desc)
        }
      }
    }
    if (is.null(.info)) return(NULL)
    list(est=.est, type=.info$type, description=.info$description)
  })
  .rows[!vapply(.rows, is.null, logical(1))]
}

#' Format the tagged estimation methods as colored, grouped display lines
#'
#' @param current optional est= string typed by the user, highlighted if found
#' @return character vector of display lines
#' @noRd
.nlmixr2EstTypeLines <- function(current=NULL) {
  .rows <- .nlmixr2EstTypeTagged()
  if (length(.rows) == 0L) return(character(0))
  .uTypes <- unique(vapply(.rows, `[[`, character(1), "type"))
  .ord <- c(.nlmixr2EstTypeOrder[.nlmixr2EstTypeOrder %in% .uTypes],
            sort(setdiff(.uTypes, .nlmixr2EstTypeOrder)))
  unlist(lapply(.ord, function(.ty) {
    .sub <- Filter(function(.r) .r$type == .ty, .rows)
    c(paste0(cli::symbol$bullet, " ", crayon::bold(.ty)),
      vapply(.sub, function(.r) {
        .name <- if (!is.null(current) && identical(.r$est, current)) {
          crayon::yellow(.r$est)
        } else {
          crayon::blue(.r$est)
        }
        paste0("   ", cli::symbol$line, " ", .name, " -- ", .r$description)
      }, character(1)))
  }), use.names=FALSE)
}

#' Print the tagged, category-grouped estimation methods to the console
#'
#' Used when `nlmixr2()` is called with no arguments, and shares its formatting
#' with the unsupported-`est=` error.
#' @return the `nlmixr2AllEstType()` data frame, invisibly
#' @noRd
.nlmixr2EstTypePrint <- function() {
  .lines <- .nlmixr2EstTypeLines()
  if (length(.lines) > 0L) {
    message("nlmixr2 estimation methods (specify with `est=`):\n",
            paste(.lines, collapse="\n"))
  }
  invisible(nlmixr2AllEstType())
}

#' Tagged list of the available nlmixr2 estimation methods
#'
#' Returns the built-in (and any attribute-tagged third-party) `est=` methods
#' grouped by their estimation category, as used when an unsupported method is
#' requested.
#'
#' @return data.frame with columns `est`, `type` and `description`
#' @examples
#' nlmixr2AllEstType()
#' @export
nlmixr2AllEstType <- function() {
  .rows <- .nlmixr2EstTypeTagged()
  .uTypes <- unique(vapply(.rows, `[[`, character(1), "type"))
  .ord <- c(.nlmixr2EstTypeOrder[.nlmixr2EstTypeOrder %in% .uTypes],
            sort(setdiff(.uTypes, .nlmixr2EstTypeOrder)))
  .rows <- .rows[order(match(vapply(.rows, `[[`, character(1), "type"), .ord))]
  data.frame(
    est=vapply(.rows, `[[`, character(1), "est"),
    type=vapply(.rows, `[[`, character(1), "type"),
    description=vapply(.rows, `[[`, character(1), "description"),
    stringsAsFactors=FALSE
  )
}
