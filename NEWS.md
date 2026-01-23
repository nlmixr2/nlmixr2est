# nlmixr2est (development version)

- Optimize `%in%` operations throughout the codebase for better
  performance by replacing with `match()` for O(n) instead of O(n*m)
  complexity. This significantly improves performance in iterative
  algorithms like SAEM and FOCEI.

- Add `predict(fit, level="ipred")`, `predict(fit,
  level="individual")` or `predict(fit, level=1)` to predict
  individual fits (with possibly a new dataset).

- Change test files to `.rds` files

- Drop magrittr `%>%` in favor of `|>`.

- **Breaking change:** Minimum R version increased from 4.0 to 4.1.0.
  This change is required to support the native pipe operator `|>`.
  Users on R < 4.1.0 will need to upgrade R to install this version
  of nlmixr2est.


- Turn back on object compress and use roxde2's default for
  compression.

- Bug fixes for deparsing nlmixr2 control objects

# nlmixr2est 5.0.0

- Remove `qs` and change to `qs2`.  This breaks backward
  compatibility.

- Default to non-compressed nlmixr2 objects

# nlmixr2est 4.1.1

- Request nlmixr2est's pre-processing hooks for `augPred()`, `vpcSim()` and
  `$simInfo`, which fixes augPred in cases where `etas=0` are used in
  `nlmixr2` (#587)

- Fix scale.h so that `scaleType="none"` does not also require
  `scaleTo=0`

- Request Armadillo 15 with the special flag in the new `RcppArmadillo`

- Fix `focei` without etas (and without log-likelihood normal) to run
  `ELS` (See #590).

- Change the IOV implementation (#596):
   - Now shows estimates as `CV%` or `sd` without shrinkage calculation.
   - Allow different forms of `iov` estimation, controlled by
     `iovXform`.
   - Retains the `iov` parameter(s) in the output `data.frame`.
   - With `iov`, the `$omega` shows a list of variability by the
     conditioning variable(s).
   - `fit$iov` will show the IOV deviations by the conditioning
     variables(s) with the exception of `id`
   - IOV models can be used in other estimation methods and inherits
     the ETA values.

 - Added `$etaMat` method for `nlmixr2` fits to give the value that
   needs to be passed between each estimation method (related to iov #596)


# nlmixr2est 4.1.0

- Updated inferring the estimation method from the control
  object. Requires the control object to have a class of length one
  and match the estimation method.  For example `foceiControl()` would
  assume that the estimation method is related to `focei`.

- Changed Rstudio completion to not evaluate (in case it gets turned
  on for data.frames) (See #568)

- Turned on data completion for items like `$fitMergeInner`

- **Breaking change:** Changed the estimation method `posthoc` to add
  tables and calculate the covariance by default.  It is now a method
  with it's own control, `posthocControl()`.  As previously the
  default is not to include the interaction term (but you can turn it
  on with `posthocControl(interaction=TRUE)`).

- Added `foceControl()`, `foControl()` and `foiControl()` for the
  `foce`, `fo` and `foi` methods, respectively.  They try to convert
  the related control structures to the correct control structure for
  the estimation method.

- Added iov support for `focei`,  `foce`, and `saem` (#614)

- Added new estimation method `agq` which uses adaptive Gauss-Hermite
  Quadrature to fit a nonlinear-mixed effect model. In this method,
  you can choose the number of quadrature points to estimate the
  likelihood, with higher numbers giving more accurate likelihoods.
  The AGQ implementation in nlmixr2est allows you to specify the
  number of quadrature points via the `agqControl()` function, and
  supports both single and multiple subject models. This method is
  particularly useful for models where accurate likelihood estimation
  is critical.

- Also added a `laplace` method which is the same as
  `agq` with 1 node (and is numerically the same as `focei`, `foce` or
  log-likelihood `focei`/`laplace`, etc), but uses the `agq` routine.

- Fixed saem mu-reference display by not compressing the internal item
  `saem0`.

# nlmixr2est 4.0.2

- The loading and unloading of DLLs has been minimized in this version
  of nlmixr2est. This avoids loading/reloading the same DLLs and causing the
  CRAN mac m1 ASAN/USBAN false positive issue observed in CRAN.

- Additionally a new function `nlmixr2fix(fit)` has been added to
 `nlmixr2est`.  It attempts to make the fit loaded from a different
 version of nlmixr2 compatible with nlmixr2 4.0.  It also prints out
 the versions of `nlmixr2` that were used when creating this fit.
 With this information you are more likely to find a way to use the
 fit in your current session (or in an old session). (Issue #562)

# nlmixr2est 4.0.1

- Initialize lbfgsb3 error message to an empty string to address
  valgrind finding (as requested by CRAN).

# nlmixr2est 4.0.0

- When using a model to start a new focei model, the ETAs from the
  last fit are used as the starting point.  Now you can use
  `foceiControl(etaMat=NA)` to skip this and use `eta=0` for all
  items.

- When using `foceiControl(etaMat=fit)`, this will extract the ETAs
  from a fit for use in the next optimization.

- When using a `foceiControl(etaMat=)` option nlmixr2 no longer only
  evaluates the inner problem with the `etaMat` value.

- Add `mceta` option to `"focei"`.

  - `mceta=-1` is the default; the eta restarts at the best eta from
    the last step to start the inner optimization.
  - `mceta=0` the eta starts at `0` to start the inner optimization.
  - `mceta=1` the eta starts at either `0` or the best `eta`, which
    ever gives the lowest objective function to start the inner
    optimization.
  - `mceta=n` under the assumption of `omega` sample `n-1` `eta`
     values and use the lowest objective function of eta sampled, last
     best eta and eta=0 to start the inner optimization.

- Fix Rstudio print (issue #536)

- Support rxode2's new `+var()` definition in `saem`

- Support literal fixing of residuals (#524).  All methods that
  support a literal fix of residuals have an option `literalFixRes`
  which defaults to `TRUE`.  To get the behavior from older models you can use
  `literalFixRes=FALSE`

# nlmixr2est 3.0.4

- More robust covariance calculation in `focei`.

- Allow hook mechanism to handle piped arguments.

- Fix for when output message from optimizing doesn't print well
  (#325)

# nlmixr2est 3.0.3

- Moved data check for covariates and required data items to a
  pre-processing step. This fixes #499.  Each method that needs to
  have a covariate check needs to have a property `covPresent`. For
  example to apply the covariate data check to the `focei` method you
  need `attr(nlmixr2Est.focei, "covPresent") <- TRUE`.

- Bug fix for non-mu referenced etas when combined with mu referenced
  covariate values. (See #498)

- Changed option for `"saem"` to have `literalFix=FALSE`. This makes
  mu-referencing work better when fixing a population value.

# nlmixr2est 3.0.2

- Fix bug where models where omega boundary warnings caused problems
  in estimation (#490)

- Created a new api for pre-processing ui, allowing adding arbitrary
  hooks.  As written now, this includes literal fix and zero omega as
  well as added the new rxode2 ui processing.

- Fixed compilation to only use -I in most systems for maximum
  compatibility

# nlmixr2est 3.0.1

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

# nlmixr2est 3.0.0

- No binary linking to `rxode2`, `lbfgsb3c` and `n1q1`, which means
  that updating these will not make `nlmixr2est` crash without
  recompiling.

- New `mu`3 referencing will take context from the model to see if the
  algebraic expression can be completed from defined model variables;
  These variable would have to be unique.

# nlmixr2est 2.2.2

## Breaking changes

- Saem non-mu reference input parameters/covariates were fixed so they
  work correctly with fixed parameters (Issue #445)

- Focei changed back to having a lower bound for standard deviations
  when not specified. This means that best model fits may change.  You
  can revert to the old settings by using
  `foceiControl(sdLowerFact=0.0)`.  You can also change the factors to
  other values than the default value, that is
  `foceiControl(sdLowerFact=0.000001)` for instance which would
  multiply the initial value by `0.000001` when either the lower bound
  isn't specified or the lower bound is specified as zero for the
  error estimates related to error-based standard deviations.

- In `nlmixr2`, expressions are optimized.  Because of that
  optimization, numerical rounding differences can cause different
  directions in optimization when fixing parameters in the model
  vs. fixing the parameters manually.

  This means that the fixed parameters in a model vs hard-coded fixed
  parameters could give different values in the final model.

  A new option `literalFix` was introduced which change the fixed
  population parameters to constants in the model while running the
  optimization.  This makes the output of fixing within the model and
  fixing manually the same (which is what is likely expected). The
  default is for this to be turned on (ie. `literalFix=TRUE`).  You
  can get back the old behavior by using the option
  `literalFix=FALSE`.

- In `saem`, the monte-carlo sampling occurs for all parameters
  including non-informative ETAs.  A fix ensure that non-informative
  etas in `saem` are fixed to zero while sampling the `phi` values.
  This may change results for models with uninformative etas. To
  ignore the uninformative etas with `saem` you ca use use the prior
  `saem` handling with `saemControl(handleUninformativeEtas=FALSE)`.

## New features

- Gracefully degrade when $cov is not in the right form (see #423)

- Add support for PopED in place solving (used in babelmixr2)

- If `est=foceiControl()` or other nlmixr2 control with the class
  `foceiControl` infer the estimation method is `focei`

- Add back the warnings when estimation methods ignore the boundaries

- When using `rxSolve`, now respects the values from `tableControl()`
  (#465 and #297)

## Bug fixes

- Will emit warnings when the return object is not a nlmixr2 fit
  (#453)

## Other things

- Moved actual code of some matrix libraries to `lotri` and import
  them via function pointers

# nlmixr2est 2.2.1

- Align with the possibility that linCmt sensitivities may not be
  present (like intel c++)


## Bug fix
- `focei` cache needs to be based on the parameter order as well as
  the model information (#415)

# nlmixr2est 2.2.0

## New Features

- Algebraic mu referencing has been implemented in `nlme` and `saem`.

- New estimation method "nlm" has been added to estimate population
  only likelihoods using `stats::nlm` and possibly return a
  standardized `nlmixr2` fit.

- New estimation method "nls" has been added to estimate population
  only problems.  This uses `minpack.lm::nlsNM` by default if
  present, or the `stats::nls`

- New estimation method "optim" has been added to estimate population
  only likelihoods.  This uses `stats::optim` and returns a
  standardized `nlmixr2` fit.

- New estimation method "nlminb" has been added to estimate population
  only likelihoods.  This uses `stats::nlminb` and returns a
  standardized `nlmixr2` fit.

- New estimation methods from the `minqa` package: "bobyqa", "uobyqa"
  and "newuoa" have been added to estimate population only
  likelihoods.  These methods returns a standardized `nlmixr2` fit.

- New estimation method "lbfgsb3c" to estimate population only
  likelihoods.  This returns a standardized `nlmixr2` fit.

- New estimation method "n1qn1" to estimate population only
  likelihoods.  This returns a standardized `nlmixr2` fit.

- Added new feature for `vpcSim()` where a minimum number of subjects
  are simulated from the model when trying to fill in ODEs that were
  not solved successfully.  By default this is `10`.  This also
  works-around a bug when there is only one subject simulated and the
  `data.frame` has a slightly different output.

## Breaking changes

- Removed `fit$saemTransformedData` since it isn't actually used in
  `saem` anymore (but will break anyone's code who is using it)

- Now the internal function `.foceiPreProcessData()` requires the
  rxode2 control `rxControl()` because some of the new steady state
  lag features need to translate the data differently based on
  `rxControl()` options.


## Bug fixes

- Printing models with correlated omega values and omega values fixed
  to zero no longer fails (#359)

- Add back values for $parHistData (#368)

- This requires a new `rxode2` which will fix multiple endpoint issues observed (#394)

- Manual back-transformed values in `$parFixed` are now displaying
  correctly and are calculated based on the confidence interval in the
  control instead of 95% confidence no matter what (#397)

## Other changes

- An `as.rxUi()` method was added for fit models (#377)

# nlmixr2est 2.1.8

- Version bump and a minor documentation update (same as nlmixr2est
  2.1.7).  This version bump is to simply allow correct binary linkage
  to rxode2 2.0.14. Otherwise `nlmixr2` models will crash R.

# nlmixr2est 2.1.7

- As requested by CRAN, remove `Rvmmin`

- Values in `$parFixed` for BSV without exponential transformation are now
  correctly shown (#366)


# nlmixr2est 2.1.6

## Breaking changes

- Since `rxode2` now allows simulation with `omega` having diagonal
  zero elements, `$omega` and `$omegaR` now reflects this information
  including the zero omega elements in the output. On the other hand,
  the other eta-information and standard error information for zero
  etas are still excluded in `$phiR`, `$phiSE`, `$eta` etc.

## Bug fixes

- `vpcSim()` works when an eta value is fixed to 0 (#341)

- `augPred()` now consistently uses the simulation model (instead of
  the inner model used for `CWRES` calculation).

## Other changes

- Dropped dependence on orphaned package `ucminf`

# nlmixr2est 2.1.5

- Add `$fitMergeFull`, `$fitMergInner`, `$fitMergeLeft`,
  `$fitMergeRight` as a complement to `$dataMergeFull`,
  `$dataMergInner`, `$dataMergeLeft`, `$dataMergeRight`.  The fit
  variants prefer columns in the fit dataset instead of the original
  dataset.  This is useful for goodness of fit plots with censoring
  since the `DV` in the fit simulates values under the ipred/residual
  assumption and will give more appropriate goodness of fits,
  otherwise these values are the limit of whatever censoring is
  applied

- Moved the mu reference fix for the split mu referenced model here
  (from babelmixr2)


# nlmixr2est 2.1.4

- Breaking change, now calculate condition number based on covariance
  and correlation, the names have changed to be more explicit.
  `conditionNumber` changed to `conditionNumberCov` and a new metric
  `conditionNumberCor` has been added.

- A bug in boundary value detection prevented automatic covariance calculation
  with FOCEi estimation (#318)

- Fix `vpcSim` so that it will be a bit more robust when it is
  difficult to simulate.

- A bug in model piping which did not allow models to be appended to was fixed
  (rxode2#364)

- An internal change was made in `nlmixr2.rxUi()` to better support the
  babelmixr2 PKNCA estimation method (babelmixr2#75)

- Fixed bug where `$iniUi` did not return the initial ui when running
  non `focei` related methods.  Also added alias of `$uiIni` to the
  same function.

- Dropped Stan headers for this package, also updated to C++17

# nlmixr2est 2.1.3

- Allows `$etaH` and related family to be integrated into a `saem` fit
  if `cwres` is calculated.

- Fixed a bug where `nlmixrLlikObs` in the merged dataset is sometimes
  named `llikObs`, now it is always named `nlmixrLlikObs`

- Fixed a bug where `nlmixrLlikObs` shows up in merged dataset when
  `cwres` is not calculated (it was always `0`), also allow `cwres`
  calculation to pick up `nlmixrLlikObs` in merged dataset.

- Dropped `dparser` dependency

# nlmixr2est 2.1.2

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
