# RPEM -- Standard Errors

In scope for the first increment (D9).

STATUS (2026-07-09): SEs are DONE via the FOCEI-covariance path -- the eval-only
finalize (`.rpemBuildFit`) now uses `foceiControl(covMethod="r,s")`, so FOCEI
computes the R/S sandwich covariance at the fixed RPEM estimates and `parFixedDf`
reports SE / %RSE / CI. This reuses existing machinery and gives immediate
uncertainty. Two caveats: (1) these are FOCEI-covariance SEs, not the paper's
Fisher-score SEs (below) -- a later refinement; (2) M1 holds non-mu-ref
structural params fixed, yet FOCEI reports SEs for them, so those SEs are not
meaningful until the numeric fixed-effect update (D20) estimates them.

The paper's own method (Fisher-score on the converged Gaussian samples) remains
the target refinement:

## Method

1. After convergence, collect the accepted/converged Gaussian samples over the
   last consecutive iterations (paper uses iterations that are already
   stabilized, e.g. last ~5-6 collecting ~6000 samples/subject). Expose the
   collection window as a control.
2. Form the Fisher information via the score (gradient of the individual
   log-likelihood contributions) accumulated over subjects, per ref [4] section
   5. SE is proportional to `1/sqrt(n)` (paper confirms empirically).
3. Report SEs for: mixture means `mu^(k)` (and their `exp` typical values for
   log-transformed params), covariance/`Sigma^(k)` diagonal (the `sigma_*`
   in the paper's tables), weights `w^(k)`, and fixed effects.

## Reuse and consistency

- Reuse the converged-sample store already kept for predictions
  (`09-fit-object.md`) -- SEs and predictions draw from the same pool, so no
  extra ODE solves.
- Report SEs on the same scale nlmixr2 reports other methods (transformed-scale
  omega + back-transformed typical values), so tables line up with SAEM/FOCEI.
- For log-normal params, provide the typical-value SE via delta-method /
  error propagation from the `mu` SE, matching the paper's "standard error of
  the typical values" note.

## Validation

- Compare RPEM SEs against the paper's reported `+/-` values on the analytic and
  voriconazole models (the `11-validation.md` reproduction bar).
- Sanity: SEs shrink ~`1/sqrt(n)` when `n` grows (paper Table 1 shows ~10x
  smaller SE from n=100 to n=10000).

## Open items

- OI-1: Exact estimator form from ref [4] section 5 (Fisher score) -- pull the
  formula before coding; confirm whether it needs per-sample scores or can reuse
  the stored per-observation likelihood pieces.
- OI-2: Collection-window defaults (how many terminal iterations, how many
  samples/subject) as controls.
