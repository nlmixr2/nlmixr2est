# Submission to allow log-likelihood estimation

- Note checked with win-builder and ubuntu development R as well as
  windows, mac and linux current release.
  
- CRAN issues with `anonymous non-C-compatible type given name for
  linkage purposes by typedef declaration; add a tag name here` where
  changed to simple `C++` structures to suppress the "significant" warning.
  
- Also linkage from 'rxode2' to be strict removing the clang error
  messages

## Breaking changes

### FOCEi

 - Gill forward differences will not repeat now (by default), You can
   change back to prior behavior with `foceiControl(repeatGillMax=3)`
 
 - Number of sticky recalculation is reduced to 4; to have the old
   behavior use `foceiControl(stickyRecalcN=5)`
   
 - `n2ll` has been changed to `ll` to specify individual
   log-likelihoods.  This was only used in simulation and was not well
   documented.
   
 - log-likelihood requires a more recent version of `rxode2`.  Since
   it takes a while to compile, this version was made compatible with
   the old version of rxode2 while compilation issues will be
   addressed later.
 
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

