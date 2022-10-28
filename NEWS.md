# nlmixr2est (development version)

- Fixes `$etaH` memory corruption so the standard errors of etas are now correct

- Removed the memory requirements for focei by `neta*neta*nsub`

- Fixed character based covariates so the work correctly (again) with
  focei.  Added a test for this as well.

# nlmixr2est 2.1.1

- Fixes `$dataMergeInner` so that observation-based log-likelihoods
  work with infusions.  Should fix tests with `ggPMX`

- Fixes `$etaSE` and `$etaRSE` to work correctly when there is only 1
  eta.

- Fixes npde valgrind observed on CRAN machines

# nlmixr2est 2.1.0

## Breaking changes

### FOCEi

 - Gill forward differences will not repeat now (by default), You can
   change back to prior behavior with `foceiControl(repeatGillMax=3)`
 
 - Number of sticky recalculation is reduced to 4; to have the old
   behavior use `foceiControl(stickyRecalcN=5)`
   
 - `n2ll` has been changed to `ll` to specify individual
   log-likelihoods.  This was only used in simulation and was not well
   documented. 
   
 - Generalized log-likelihood is only supported with `rxode2` `2.0.8` or later.
 
### FOCEi covariance calculation

 - The `S` matrix calculation was made a bit more robust to errors in
   individual gradients.  When there are errors in the individual
   gradient calculation, assume the gradient is the same as the
   overall gradient.  In the tests cases, were reasonable using this
   adjusted S matrix.  This means if some individuals do not have very
   much data to support a specific parameter, a `S` matrix calculation
   for the population will still be generated. When there is some
   patients/subject combinations that do not have sufficient data, we
   will add the following to the run information: `S matrix had
   problems solving for some subject and parameters`. The `S` matrix
   calculation will still fail if the percentage of parameters that
   are being reset is lower than `foceiControl(smatPer=0.6)` or
   whatever you specify.
 
 - The `r,s` covariance matrix will now also check for unreasonably
   small values (controlled by `foceiControl(covSmall=...)`) and
   select a different covariance estimate method even when the "r" and
   "s" matrices are calculated "correctly".

## New features

- What type(s) censoring (if any) is now stored in `fit$censInformation`

- Standard errors of `$etas` can now be obtained with `fit$phiSE`,
  also available are `fit$phiRSE` (relative standard error),
  `fit$phiH`, (individual hessian), `fit$phiC` (individual
  covariances), `fit$phiR` (individual correlation matrices)
  
- Can also use Shi 2021 differences in addition to Gill differences.
  In our tests (using the same datasets as CPT) these produced worse
  estimates than the Gill 1983, though it is unclear why since it
  should be a faster more accurate method.  A modified version is used
  in calculating the individual Hessians of numerically for the
  generalized likelihood approach.
  
- Generalized likelihood estimation is now present in `nlmixr2est` for
  `focei`, `foce` and `posthoc`
  
- `nmNearPD()` is a function you may use for nearest positive definite
  matrix.  This is derived from `Matrix::nearPD()` but is implemented
  in C/C++ to be used in (possibly threaded) optimization.
  
- Individual Hessians can be accessed by `$phiH`, covariance by
  `$phiC`, eta standard errors by `$phiSE` and eta RSEs can be
  accessed by `$phiRSE`.  There are `eta` aliases for these as well
  (`$etaH`, `$etaC`, `$etaSE`, and `$etaRSE`).
  
- Can now access the individual point's contribution to the overall
  likelihood when merging to the original dataset. These merges can be
  accessed with `$dataMergeFull`, `$dataMergeLeft`, `$dataMergeRight`,
  and `$dataMergeInner`.  The columns with the individual data column
  is `nlmixrLlikObs`.
  
  To calculate the total `focei`/`foce` objective function, the sum of the
  likelihoods still need to be adjusted by the omega/eta contribution,
  and the individual Hessians, and possibly the NONMEM objective
  function offset constant.
  
## Censoring fixes

 - Fixed bug where datasets with censoring that are not lower case `cens` and `limit` do not
   produce the correct table output (#180)

## FOCEi updates

- Resets now scale properly when a value is simulated outside the limit
- Models with zero gradients on the first step now switch to `bobyqa`
  by default.  With this, it is more important to examine the model
  parameters and fits for plausibility.

# nlmixr2est 2.0.8

## New features

- Add `pd`/`npd` as an output as well as `npd`/`npde`

## SAEM bug fix

- When loading a `nlmixr2` "saem" fit from another R session,
  `nlmixr2` will no longer crash with `fit$objf`

## NPDE/NPD fixes

- `NPDE` was identical to `NPD` even with correlated models, this was
  fixed (prior output was actually `NPDE`).

## Censoring fixes

- FOCEi censoring fixes:
  - M4 method equation bug fix
  - M4 method derivative change based on equation fix
  - M2 method added missing derivative 
  - Censoring already dTBS

- SAEM Censoring fixes:
  - SAEM method M4 method equation bug fix
  - Censoring limit changed to dTBS

- Censoring handling was unified

## Internal changes

- Added `ui$getSplitMuModel` which is used in `babelmixr2` and will be
  used in the refined stepwise covariate selection of `nlmixr2extra`

- Added work-around to remove `_nlmixr2est_RcppExport_registerCCallable`
  since the registering of C callable are handled manually at the moment.

# nlmixr2est 2.0.7

- Use `.zeros()` for the matrices in armadillo in addition to relying
  on `calloc` to give zero matrices.

- Fixed one uninitialized object

- Fix for `augPred` so it works on population only models

- `nlme` no longer sets options to treat all covariates as non
  mu-referenced covariates, but directly calls a function that can
  turn on or off the mu-reference covariate selection.

- `vpcSim` now tries to simulate IDs that didn't simulate correctly (with a warning)

- Export nmObjHandleControlObject

# nlmixr2est 2.0.6 -- new package

`nlmixr2est` contains the estimation functions within `nlmixr2`.

## FOCEI family changes

- Remove lower level `foceiFit` function.  Focei, foce, fo, foi, and
  posthoc now directly takes rxode2 ui objects

- New error types are supported in focei including mixing theta and
  etas in residual errors and different types of proportional errors

- Different types of additive and proportional errors can be used for
  each endpoint using ` + combined1()` or `+ combined2()` otherwise it
  takes the supplied `addProp` option to figure out which type of
  combined model is run (by default `combined2()`)

- Focei model cache is now named `focei-md5Digest.qs` and uses `qs`
  compression/saving/loading.

- `foceiControl()` aligned between other methods.

- `foceiControl(adjLik=TRUE)` uses the NONMEM-style objective function
  throughout.  `foceiControl(adjLik=FALSE)` uses the adjusted
  objective function throughout, and adjusts it back to the NONMEM
  objective function.

- Lag time and other between subject variability differences no longer
  calculate an ideal relative step size, but an absolute step size
  when using Gill differences (default)

- Objective function checks for infinite/NaN/NA values for the entire
  solving space and ensures no overflow occurs when calculating the
  inner hessian

## SAEM changes

- mu referencing is no longer required for `saem`; Internally non
  mu-referenced values are converted to mu referenced values and the
  converted back when calculating the nlmixr2 object.

- `nlmixr2` forced the parameter ordering to (1) population effects,
  (2) non mu-referenced between subject effects (3) omega estimates
  and (4) residual effects. This changes the order that `nlmixr2` sees
  the parameters. Since this is based on a random number generator,
  the optimization trajectory will be different and have different
  results than `nlmixr`

- Components of `omega` can now be fixed.

- Residual error components can also be fixed.

- When optimizing only one residual value, nlmixr2's saem uses `nlm`
  from R, which is more efficient than the nealder-meade method.

- Lower level `saem` functions (like `configsaem()`) are not exported
  because they are increasingly difficult to use and convert to
  something standard; a few methods (like `print`, `summary` etc) are
  maintained to view the lower level object and for debugging it.

- Parameter history and print-out no longer includes fixed parameters.

- The model to calculate the residuals more closely matches the model
  used for estimation to remove small rounding differences that may
  occur in the models.

- Different types of additive and proportional errors can be used for
  each endpoint using ` + combined1()` or `+ combined2()` otherwise it
  takes the supplied `addProp` option to figure out which type of
  combined model is run (by default `combined2()`)

- Parameter history and printout now uses standard deviation for
  additive only components, matching the estimation of the components.

- `rxode2` solving options are now saved in the `rxControl` part of
  the `saemControl()`.  That is
  `saemControl(rxControl=rxControl(...))`; This fixes any conflicting
  option names as well as allowing alignment between the control
  structures in `focei`, `nlme` and `saem`

- `saemControl()` aligned between other methods.


## nlme changes

- `nlme` has been completely rewritten to directly run from the
  `rxode2` UI

- `nlme` always tries to use mu-referencing (when available)

- Internally `nlme` now uses parallel processing for solving so it
  should be faster.

- `nlmixr2NlmeControl()` (which will overwrite `nlmeControl()`)
  documents and adds more options to `nlme`. Also aligned with other
  methods.

- `weights`, `fixed`, `random` can be specified in
  `nlmixr2NlmeControl()`.  If so, then the `nlme` object will be
  returned.

- `returnNlme` is a new option that will return the `nlme` object
  instead of the traditional `nlme` object.

- `nlme_ode` and `lme_lin_cmpt` are both removed.

- `rxode2` solving options are now saved in the `rxControl` part of
  the `saemControl()`.  That is
  `nlmeControl(rxControl=rxControl(...))`; This fixes any conflicting
  option names as well as allowing alignment between the control
  structures in `focei`, `nlme` and `saem`

## nlmixr2 object change

- With `saem`, the nlmixr2 function now saves/compresses the `phiM`
  information.  This means the gaussian and Laplacians likelihoods can
  be calculated when you save the nlmixr object and then restore it
  later.

- The nlmixr2 object compresses infrequently used and removes many
  unneeded objects. Even with compression, the `saem` objects are
  often a bit bigger since they include the large `phiM` object.

- `nlmixr2` now supports non-mu referenced ETAs in the `fit$parFixed`
  and `fit$parFixedDf`

## nlmixr2 interface change

- `nlmixr2` interface changed to use `rxode2` UI

- `keep` and `drop` are added to `tableControl` to influence the end data-frame

- `$simInfo` uses a quoted expression for `$rx` instead of a string

- `$simInfo$sigma` is a diagonal matrix since now the normal
  simulation is controlled by the variability modeled as a population
  value.

- `nlmixr2` now allows etas that have initial omega estimates of zero
  to be dropped from the model (instead of issuing an error about a
  non-positive definite `$omega` matrix)

## NPDE changes

- Fixed a bug where the number of simulations for a NPDE calculation
  are correctly passed by `addNpde(fit, table=tableControl(nsim=500))`

## VPC changes

- `vpc` function rewritten and split out to `vpcSim()` and
  `vpcPlot()` (which is a replacement for `vpc()`).

- There were too many mismatches between `vpc::vpc` and `nlmixr::vpc`
  which caused inconsistencies in code based on load order of `vpc`
  and `nlmixr`.  This way both coexist, and you can use the `vpc`
  simulation for other packages more easily (like `ggPMX`) without
  creating or summarizing data since `ggPMX` has its own methods for
  summarizing and creating plots.

- VPC now directly uses `rxode2::rxSolve`

## augPred() changes

- `augPred()` has been written to use the new fit object.

- `nlmixr2AugPred` was changed to `nlmixr2AugPredSolve()`

- `augPred` uses the new interface and supports multiple endpoints.
  The endpoint name is now always on the `plot(augPred(fit))`.

## getFitMethod() change

- Internally, fit estimation method is saved in `fit$est`, and now
  `getFitMethod(fit)` simply returns `fit$est`

## Delete methods

- Many methods lower level utility functions have been deleted.

- `nmDocx`, `nmLst` and `nmSave` have been removed.

## Bug fixes

- Now will reset the cache when items cannot be loaded. In the past
  error messages like `function
  'rx_0ba247452048de33b1ffb8af516714fc__calc_lhs' not provided by
  package 'rx_0ba247452048de33b1ffb8af516714fc_'` would cause the
  estimation to stop.  Now `rxode2::rxClean()` is run when this occurs.

