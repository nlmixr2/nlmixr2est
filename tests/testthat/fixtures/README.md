# Test Fixtures - Cached Model Fits

This directory contains cached `.rds` files for commonly used model fits in the test suite.

## Purpose

The test suite repeatedly fits the same models across multiple test files, which is time-consuming. These cached fits:
- Speed up test execution by computing fits once and reusing them
- Save developer time during test runs
- Improve CI/CD pipeline efficiency

## Cache Strategy

Fits are cached with a version identifier that includes:
- Cache format version (currently v1.0.0)
- nlmixr2est package version

Caches are automatically invalidated when:
- The package version changes
- The cache format changes
- Cache files are manually deleted

## What's Cached

See `helper-fits.R` for the complete list of cached fits. Common ones include:
- `one.compartment.fit.saem` - One compartment model with SAEM estimation
- `one.compartment.fit.focei` - One compartment model with FOCEI estimation
- `one.compartment.with.lag.fit.focei` - Lag model with FOCEI
- `one.compartment.with.lag.fit.saem` - Lag model with SAEM

## Git Ignore

This directory is excluded from version control via `.gitignore` because:
- Fit objects can be large (MB range)
- They are automatically regenerated on first test run
- Different R versions/platforms may produce slightly different results

## Manual Cache Management

To clear the cache and force recomputation:
```bash
rm -rf tests/testthat/fixtures/*.rds
```

Caches will be automatically recreated on the next test run.
