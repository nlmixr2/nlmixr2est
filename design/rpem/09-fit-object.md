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

## Concrete integration plan (mirror SAEM's finalize)

RPEM has random effects, so the full `nlmixr2FitData` is assembled the SAEM way,
not the nlm way: estimate the population parameters, set them on the UI, then run
FOCEI in eval-only mode to compute EBEs/residuals/tables at those fixed values.

Steps (see `R/saem.R` `nlmixr2Est.saem` lines ~1085-1122 as the template):
1. `.rpemFit` already returns `mu`/`omega`/`add.sd` and per-subject EBEs
   (`$ebe`, posterior-mean etas via Eq 53 -- DONE).
2. Build `env$fullTheta` (full theta vector in UI order: RPEM `mu` for the
   mu-referenced thetas, `add.sd`, and the held-fixed structural thetas at ini),
   `env$omega` (etaNames x etaNames from `omega`), `env$.etaMat`/`.etaMatBase`
   (the N-row EBE matrix), and `env$etaObf` (ID + eta cols + OBJI). Mirror
   `.getSaemTheta`/`.getSaemOmega`.
3. `.nlmixr2FitUpdateParams(env)` to push estimates into `ui$iniDf`.
4. `rpemControlToFoceiControl`: `foceiControl(maxOuterIterations=0,
   maxInnerIterations=0, covMethod=0, etaMat=<EBEs>, calcTables=..., est="rpem")`
   so FOCEI only evaluates at the RPEM estimates + supplied EBEs.
5. `nlmixr2CreateOutputFromUi(ui, data, control, table, env, est="rpem")`.
This yields CWRES/NPDE/IPRED/tables for free via the FOCEI eval path.

STATUS (2026-07-09): `.rpemBuildFit` implements steps 1-5 and the FOCEI-eval
assembly SUCCEEDS to a correct fit CORE -- `parFixedDf` (RPEM mu/add.sd, held
tcl/tv), `objDf` (OBJF/AIC/BIC/logLik), `omega`, per-subject EBEs (`ranef`),
shrinkage. Two fixes got this far: (a) do NOT call `nmObjHandleControlObject`
with the focei control (it nulls `env$control` so `.foceiFamilyControl` hits the
buggy `do.call(type)`); (b) the eval-only `foceiControl(maxOuter=0, maxInner=0,
etaMat=EBEs)`.

REMAINING BLOCKER: the per-observation residual table fails --
`addTable` -> `.calcTables` -> `.foceiPredIpredList` ->
`.foceiSolvePars(fit, .ipredModel, thetaEtaParameters$ipred, ...)` errors with
"argument is of length zero" because `fit$foceiThetaEtaParameters$ipred`
(`nmObjGet.foceiThetaEtaParameters`, built from `fit$ranef`/`fit$fixef` via
`_nlmixr2est_nlmixr2Parameters`) comes back empty. So `nlmixr2CreateOutputFromUi`
returns the bare `nlmixr2FitCore` instead of the wrapped `nlmixr2FitData`.
Next: figure out why `_nlmixr2Parameters(fixef, ranef)` yields an empty `$ipred`
for the RPEM fit (likely a missing `dataSav`/ipred-model or ranef-format field
the eval-only FOCEI pass didn't populate). Until then, `.rpemBuildFit` requires a
full `nlmixr2FitData` and `est="rpem"` falls back to the lightweight estimates
object (which stays fully working and tested).

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
