# CRAN Comments for nlmixr2est 3.0.1

## New features

- Now when optimizing only a single parameter with `focei`-family,
  will change to use `stats::optimize()` for the outer problem (#481)

- When estimating with all fixed population parameters, do a posthoc
  estimation.

- Internally removed `assignInMyNamespace()` replacing with
  `nlmixr2global`, which fixes some edge case bugs where the nlmixr2
  environment was not reset properly.

- Treated edge case where all initial parameters are zero and change
  scaling from scaled to unscaled (#486)

- Added `mu`4 referencing that will change string expressions to
  `rxode2` numeric values.  This allows derived strings to also be
  treated as `mu` expressions (#484)

## Bug Fixes

- Fix `focei` covariance step when many `omega` values are fixed #482
