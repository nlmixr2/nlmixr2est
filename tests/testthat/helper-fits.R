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
# - When you need to modify the fit object (create a copy first: fit_local <- fit)

# Version identifier for cache invalidation
# Increment this when models, data, or nlmixr2est version changes significantly
.fit_cache_version <- paste0("v1.0.0-nlmixr2est-", packageVersion("nlmixr2est"))

# Cache directory for fit objects
.fit_cache_dir <- file.path(testthat::test_path(), "fixtures")

# Ensure cache directory exists
if (!dir.exists(.fit_cache_dir)) {
  dir.create(.fit_cache_dir, recursive = TRUE)
}

#' Helper function to load or create cached fit
#'
#' @param name Human-readable name for the fit (for logging)
#' @param fit_fn Function that computes the fit (called if cache miss)
#' @param cache_file Name of the cache file (e.g., "fit-one-compartment-saem.rds")
#' @return The fit object (either from cache or freshly computed)
.get_cached_fit <- function(name, fit_fn, cache_file) {
  cache_path <- file.path(.fit_cache_dir, cache_file)
  
  # Try to load from cache
  if (file.exists(cache_path)) {
    tryCatch({
      cached <- readRDS(cache_path)
      if (!is.null(cached$version) && cached$version == .fit_cache_version) {
        message(sprintf("✓ Loading cached fit: %s", name))
        return(cached$fit)
      } else {
        message(sprintf("⚠ Cache version mismatch for %s (expected %s, got %s)", 
                       name, .fit_cache_version, cached$version))
      }
    }, error = function(e) {
      warning(sprintf("Failed to load cache for %s: %s", name, e$message))
    })
  }
  
  # Compute fit if cache miss or invalid
  message(sprintf("⚙ Computing fit: %s (this may take a while...)", name))
  fit <- fit_fn()
  
  # Save to cache
  tryCatch({
    saveRDS(list(version = .fit_cache_version, fit = fit), cache_path)
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
one.compartment.fit.saem <- .get_cached_fit(
  name = "one.compartment.fit.saem",
  fit_fn = function() {
    .nlmixr(one.compartment, theo_sd, est = "saem", control = saemControlFast)
  },
  cache_file = "fit-one-compartment-saem.rds"
)

#' One compartment FOCEI fit (second most common)
#' Used in: test-broom.R, test-nmobj.R, test-predict.R, test-keep.R, etc.
one.compartment.fit.focei <- .get_cached_fit(
  name = "one.compartment.fit.focei",
  fit_fn = function() {
    .nlmixr(one.compartment, theo_sd, est = "focei", 
            control = foceiControl(print = 0, maxOuterIterations = 0L))
  },
  cache_file = "fit-one-compartment-focei.rds"
)

# =============================================================================
# ONE COMPARTMENT WITH LAG MODEL FITS
# Uses: one.compartment.with.lag model from helper-models.R
# Data: warfarin (filtered to dvid=="cp")
# =============================================================================

#' One compartment with lag FOCEI fit
#' Used in: test-lag.R, test-model.R
one.compartment.with.lag.fit.focei <- .get_cached_fit(
  name = "one.compartment.with.lag.fit.focei",
  fit_fn = function() {
    d <- nlmixr2data::warfarin |>
      dplyr::filter(dvid == "cp")
    .nlmixr(one.compartment.with.lag, d, est = "focei", 
            control = foceiControl(print = 0))
  },
  cache_file = "fit-one-compartment-lag-focei.rds"
)

#' One compartment with lag SAEM fit
#' Used in: test-model.R
one.compartment.with.lag.fit.saem <- .get_cached_fit(
  name = "one.compartment.with.lag.fit.saem",
  fit_fn = function() {
    d <- nlmixr2data::warfarin |>
      dplyr::filter(dvid == "cp")
    .nlmixr(one.compartment.with.lag, d, est = "saem", control = saemControlFast)
  },
  cache_file = "fit-one-compartment-lag-saem.rds"
)

message("✓ Helper-fits.R loaded successfully")
message(sprintf("  Cache version: %s", .fit_cache_version))
message(sprintf("  Cache directory: %s", .fit_cache_dir))
