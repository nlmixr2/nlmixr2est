# RPEM -- Fit Object, Residuals, Predictions

Full residual/prediction parity in the first increment (D10): a proper
`nlmixr2FitData` with population parameters, objective/log-likelihood, EBEs,
predictions, and residuals -- indistinguishable in shape from a SAEM/FOCEI fit.

## Assembly

- `R/rpem_fit.R` / `rpem_fit_aux.R`: build the `nlmixr2FitData` from the C++
  engine output, templated on `saem_fit.R`. Post-fit accessors go through
  existing `nmObjGet` / `nmObjHandle` -- add RPEM branches, do not fork the
  accessor system.
- Populate the standard slots: `theta` (typical values), `omega`
  (`Sigma^(k)`), residual-error params, `objDf`/objective (`-2*lnL` or the
  nlmixr2 objective convention -- confirm sign/constant vs SAEM), iteration
  trace, and SEs (`08-standard-errors.md`).

## EBEs / individual parameters

- RPEM's individual estimate is the posterior mean `<theta_i>_opt` (paper Eq 53)
  computed from the stored converged samples with the normalized posterior
  weights `g_ik`. Store per-subject EBEs and etas.

## Predictions

- **Population prediction** `<y_pred>_pop` (Eq 52): mixes over components using
  `w^(k) p(theta|mu^(k),Sigma^(k))`. For the population case, guard unrealistic
  `theta` (paper sets likelihood to zero for non-physical draws / truncation
  factor `N^(k)`); with transforms this is moot (see `03-parameter-space.md`).
- **Individual prediction** `<y_pred>_ind` (Eq 54/58): plug `<theta_i>_opt` into
  the model, or the posterior-expectation form. Provide PRED and IPRED columns.

## Residuals and diagnostics

- Reuse the lazy residual machinery: `addCwres()` (`src/cwres.cpp`),
  `addNpde()` (`src/npde.cpp`), `ires`/`res`. RPEM supplies the fitted
  predictions and the population/individual parameter draws these functions
  need; it does not reimplement CWRES/NPDE.
- VPC: works through the standard nlmixr2 simulation path once the fit carries
  the population model + parameter uncertainty; the paper's VPC (Eq VPC / Fig 4)
  is then just `vpcResultsUI`-compatible output.

## Store kept for post-fit

The converged-sample store (theta samples, per-sample likelihoods, labels) from
the terminal iterations is retained on the fit (or recomputable) to serve SEs,
EBEs, and predictions from one source of truth.

## Open items

- OI-1: Objective definition and constant: align `-2 lnL` reporting with how
  nlmixr2 reports SAEM objective so cross-method comparison in `11-validation.md`
  is apples-to-apples.
- OI-2: Confirm CWRES/NPDE inputs RPEM must provide match what `addCwres`/
  `addNpde` expect from a SAEM fit (population predictions, variance pieces).
