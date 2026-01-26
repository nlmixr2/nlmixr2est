# Common model fits for testing
# This file is loaded before tests run (helper- prefix ensures early loading)
# Fits are cached to disk to avoid recomputation across test runs
#
# This dramatically reduces test execution time by:
# 1. Computing fits once during test suite initialization
# 2. Reusing the same fit objects across multiple test files
# 3. Saving fits to disk (.rds files) to avoid recomputation even on first run
#
# WHEN TO USE CENTRALIZED FITS:
# - Use these fits when testing fit properties, not the fitting process itself
# - Only for standard model/data/estimator combinations
# - When you don't modify the model or data within the test
# - When you don't require specific control parameters (beyond standard saemControlFast)
#
# WHEN TO KEEP LOCAL FITS:
# - Testing the estimation process itself
# - Using custom control parameters
# - Modifying models (e.g., with model() or ini() pipes)
# - Testing error conditions or edge cases
# - When you need to modify the fit object (create a copy first: fitLocal <- fit)

# Version identifier for cache invalidation
# Increment this when models, data, or nlmixr2est version changes significantly
.fitCacheVersion <- paste0("v1.0.2-nlmixr2est-", utils::packageVersion("nlmixr2est"))

# Cache directory for fit objects
.fitCacheDir <- file.path(testthat::test_path(), "fixtures")

# Ensure cache directory exists
if (!dir.exists(.fitCacheDir)) {
  dir.create(.fitCacheDir, recursive = TRUE)
}

#' Helper function to load or create cached fit
#'
#' @param name Human-readable name for the fit (for logging)
#' @param fitFn Function that computes the fit (called if cache miss)
#' @param cacheFile Name of the cache file (e.g., "fit-one-compartment-saem.rds")
#' @return The fit object (either from cache or freshly computed)
.getCachedFit <- function(name, fitFn, cacheFile) {
  cachePath <- file.path(.fitCacheDir, cacheFile)

  # Try to load from cache
  if (file.exists(cachePath)) {
    tryCatch({
      cached <- readRDS(cachePath)
      if (!is.null(cached$version) && cached$version == .fitCacheVersion) {
        message(sprintf("✓ Loading cached fit: %s", name))
        return(cached$fit)
      } else {
        message(sprintf("⚠ Cache version mismatch for %s (expected %s, got %s)",
                       name, .fitCacheVersion, cached$version))
      }
    }, error = function(e) {
      warning(sprintf("Failed to load cache for %s: %s", name, e$message))
    })
  }

  # Compute fit if cache miss or invalid
  message(sprintf("⚙ Computing fit: %s (this may take a while...)", name))
  fit <- fitFn()
  AIC(fit)  # Force evaluation of AIC

  # Save to cache
  tryCatch({
    saveRDS(list(version = .fitCacheVersion, fit = fit), cachePath)
    message(sprintf("✓ Cached fit saved: %s", name))
  }, error = function(e) {
    warning(sprintf("Failed to save cache for %s: %s", name, e$message))
  })

  return(fit)
}

# =============================================================================
# ONE COMPARTMENT MODEL FITS (most common)
# Uses: one.compartment model from helper-models.R
# Data: theo_sd (from nlmixr2data)
# =============================================================================

#' One compartment SAEM fit (most commonly used in tests)
#' Used in: test-ini-ui.R, test-coerce.R, test-broom.R, test-print.R,
#'          test-saem-aic.R, test-00-reload-ll.R, and many more
one.compartment.fit.saem <- .getCachedFit(
  name = "one.compartment.fit.saem",
  fitFn = function() {
    .nlmixr(one.compartment, theo_sd, est = "saem", control = saemControlFast)
  },
  cacheFile = "fit-one-compartment-saem.rds"
)

#' One compartment FOCEI fit (second most common)
#' Used in: test-broom.R, test-nmobj.R, test-predict.R, test-keep.R, etc.
one.compartment.fit.focei <- .getCachedFit(
  name = "one.compartment.fit.focei",
  fitFn = function() {
    .nlmixr(one.compartment, theo_sd, est = "focei",
            control = foceiControl(print = 0, maxOuterIterations = 0L))
  },
  cacheFile = "fit-one-compartment-focei.rds"
)

# =============================================================================
# ONE COMPARTMENT WITH LAG MODEL FITS
# Uses: one.compartment.with.lag model from helper-models.R
# Data: warfarin (filtered to dvid=="cp")
# =============================================================================

#' One compartment with lag FOCEI fit
#' Used in: test-lag.R, test-model.R
one.compartment.with.lag.fit.focei <- .getCachedFit(
  name = "one.compartment.with.lag.fit.focei",
  fitFn = function() {
    d <- nlmixr2data::warfarin |>
      dplyr::filter(dvid == "cp")
    .nlmixr(one.compartment.with.lag, d, est = "focei",
            control = foceiControl(print = 0))
  },
  cacheFile = "fit-one-compartment-lag-focei.rds"
)

#' One compartment with lag SAEM fit
#' Used in: test-model.R
one.compartment.with.lag.fit.saem <- .getCachedFit(
  name = "one.compartment.with.lag.fit.saem",
  fitFn = function() {
    d <- nlmixr2data::warfarin |>
      dplyr::filter(dvid == "cp")
    .nlmixr(one.compartment.with.lag, d, est = "saem", control = saemControlFast)
  },
  cacheFile = "fit-one-compartment-lag-saem.rds"
)

message("✓ Helper-fits.R loaded successfully")
message(sprintf("  Cache version: %s", .fitCacheVersion))
message(sprintf("  Cache directory: %s", .fitCacheDir))
