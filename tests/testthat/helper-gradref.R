# Cached numeric references for the analytic-derivative tests (gradient AND covariance).
#
# The gradient tests compare the analytic outer gradient against a central difference of
# the objective.  Building that reference is by far the dominant cost: each parameter
# needs two POSTHOC FITS at perturbed thetas, and the inner problem has to be fully
# converged at each one (see below), so a single test can spend minutes producing numbers
# that never change.
#
# They never change because the reference is a property of the MODEL, the DATA and the
# THETA POINT -- not of the gradient implementation being tested.  That is the whole point
# of a reference.  So unlike helper-zzz-fits.R, whose key hashes R/ + src/ because a fit
# IS implementation-dependent, this cache is deliberately NOT invalidated by source
# changes.  Refactoring the gradient must not silently regenerate the thing that is
# supposed to catch the refactor going wrong.
#
# Regenerate deliberately, when the OBJECTIVE legitimately changes (a solver default, a
# likelihood fix -- not a gradient refactor):
#
#     NLMIXR2EST_REGEN_GRADREF=TRUE Rscript -e 'testthat::test_file("tests/testthat/test-focei-fast-grad.R")'
#
# and commit the updated .rds with a note saying what moved and why.
#
# Why the REFERENCE cannot simply be made cheaper: it differences the MARGINAL (Laplace)
# objective, whose log|H(eta*)| term depends explicitly on eta*, and d log|H|/d eta does
# not vanish at the inner optimum -- that is precisely why the analytic gradient carries
# the Almquist Eq-46 etaP term.  Freezing the etas, warm-starting them, or capping the
# inner iterations all perturb eta*(theta +/- h), and the shi step selection is
# state-dependent so the difference lands on a different h as well.  Measured on the
# one-compartment model: warm starting ofvAt() from the base fit's etaMat cut the time
# 10.2s -> 6.5s but moved the reference up to 1.7%, degrading analytic-vs-reference
# agreement from 0.0003 to 0.0173 -- past the 1% the tests assert.  Cache it; do not
# approximate it.
#
# NOTE that measurement says nothing about warm-starting the GRADIENT side.  It varied the
# warm start only inside ofvAt(); the analytic gradient was computed once, cold, in both
# arms, so what moved was the reference.  With the reference frozen here, warm-starting
# the remaining (base) fits is an open and promising way to cut these tests further -- the
# analytic gradient has no equivalent step-size sensitivity.  It still needs eta* actually
# converged, since Almquist assumes the inner optimum, so verify against the cached
# reference rather than assuming.

## ABSOLUTE, resolved once at source time.  testthat::test_path() is relative to the
## CURRENT directory, and helpers are sourced from the package root while the tests
## themselves run with the working directory set to tests/testthat -- so a relative path
## captured here silently points somewhere that does not exist by the time it is used.
## baselines/, NOT fixtures/: .nlmixr2estbuild() unlinks every fixtures/*.rds on
## devtools::document(), and fixtures/*.rds is gitignored.  Either one alone would defeat
## this cache -- an untracked reference is rebuilt by whatever code is checked out, which
## is precisely the regression it exists to catch.  baselines/ is the repo's documented
## home for checked-in reference data that must survive cache clearing.
.gradRefDir <- normalizePath(file.path(testthat::test_path(), "baselines"),
                             winslash = "/", mustWork = FALSE)
if (!dir.exists(.gradRefDir)) dir.create(.gradRefDir, recursive = TRUE)

.gradRefRegen <- isTRUE(as.logical(Sys.getenv("NLMIXR2EST_REGEN_GRADREF", "false")))

#' Load a cached numeric reference, or compute and cache it.
#'
#' Used for BOTH kinds of reference in these tests:
#'   * central differences of the objective, for the analytic outer GRADIENT;
#'   * the covMethod="r" finite-difference standard errors, for the analytic
#'     COVARIANCE (i.e. the analytic Hessian that builds it).
#' Both are properties of the model/data/theta and independent of the analytic
#' implementation under test, so both are safe to freeze and expensive to rebuild.
#'
#' @param name cache key; becomes fixtures/gradref-<name>.rds
#' @param fn zero-argument function returning the reference (numeric vector or matrix)
#' @return the reference
#' @noRd
.numRef <- function(name, fn) {
  .path <- file.path(.gradRefDir, paste0("gradref-", name, ".rds"))
  if (!.gradRefRegen && file.exists(.path)) {
    .r <- tryCatch(readRDS(.path), error = function(e) NULL)
    if (!is.null(.r)) return(.r)
  }
  .r <- fn()
  ## Deliberately NOT silent: a reference that cannot be written means every run pays to
  ## rebuild it, which is the entire cost this cache exists to remove.  Warn loudly.
  tryCatch(saveRDS(.r, .path),
           error = function(e)
             warning("could not cache reference '", name, "' to ", .path, ": ",
                     conditionMessage(e), call. = FALSE))
  .r
}

## Back-compat alias: the gradient tests read better with the specific name.
.gradRef <- .numRef
