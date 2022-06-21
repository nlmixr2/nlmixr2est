# Bugfix submission

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

- Added work-around to remove
  `_nlmixr2est_RcppExport_registerCCallable` for the most recent
  `Rcpp` since the registering of C callable are handled manually at
  the moment.

