
- Fixed examples that called `qs` since it isn't available on CRAN
  anymore (as requested)

- Add backward compatibility fixes for babelmixr2

Additionally

# nlmixr2est 6.0.0

- `focei`, `foce`, `fo`, `laplace`, and `agq` have all been
   successfully made thread safe and parallelized (for a single
   CPU). The default tolerance relaxation for difficult to solve ODEs
   has been changed to per individual instead of for the entire
   population (which is a breaking change, so major release).  This
   should allow more precision for a majority of the subjects in the
   optimization process.

- Add `predict(fit, level="ipred")`, `predict(fit,
  level="individual")` or `predict(fit, level=1)` to predict
  individual fits (with possibly a new dataset).

- Change test files to `.rds` files

- Drop magrittr `%>%` in favor of `|>`.

- **Breaking change:** Minimum R version increased from 4.0 to 4.1.0.
  This change is required to support the native pipe operator `|>`.
  Users on R < 4.1.0 will need to upgrade R to install this version
  of nlmixr2est.

- Bug fixes for deparsing nlmixr2 control objects

- `nlm` and related pooled methods now run in parallel (based on ID)

- Tests are optimized to reduce redundant fits and run in parallel.

- `nlm` (and related pooled optimizers: `bobyqa`, `newuoa`, `uobyqa`,
  `n1qn1`, `lbfgsb3c`, `optim`, `nlminb`) now support the same
  censoring behavior (M2/M3/M4) as FOCEI and SAEM.  The
  `$censInformation` field is populated for these fits in the same way
  as FOCEI/SAEM.

- `agqControl()` and `laplaceControl()` now have `rxUiDeparse()`
  methods so they can be saved better in packages like `nlmixr2save`
  and `shinyMixR`.

- Added new `outerOpt`; methods to `focei` and related methods (`agq`,
  `laplace`, `foce`, `fo`, `foi`): "uobyqa" and "newuoa".

- `saem` and other methods now respect bounds by default by internally
  adding the appropriate transform and then applying the
  back-transformation just before returning.

  For parameters that are mu-referenced, this breaks
  mu-referencing. When it breaks mu-referencing there is a warning
  issued.  The best practice is still to have unbounded parameters
  with mu-referencing.

  If you want to ignore this behavior you may
  use `control=list(boundedTransform=FALSE)` or for saem
  `control=saemControl(boundedTransform=FALSE)`

- The mu referencing covariate procedure was made less fragile to
  support mu referencing in conjunction with iov and bounded parameter
  transformations.
