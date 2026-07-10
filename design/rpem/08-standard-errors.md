# RPEM -- Standard Errors

In scope for the first increment (D9).

STATUS (2026-07-09): Fisher-score SEs DONE for the single-eta additive/proportional
case (the paper's native method). At the converged estimates `rpemFisherReg` (C++)
forms each subject's marginal-likelihood score from the stored samples via the
Fisher identity s_i = E_{eta|Y_i}[grad complete-data loglik], with the posterior
expectation taken over the samples using the self-normalized importance weights
w_ij = softmax_j logp_ij (samples were drawn from the prior, so these ARE the
posterior). Score blocks: coef_k = design_i,k * eta_ij/omega (typical value +
covariate coefs); residual sd = -nobs_i/sd + acc_ij/sd^3 (acc = SS additive / WSS
proportional); omega = -1/(2 omega) + eta_ij^2/(2 omega^2). The empirical Fisher
information I = t(S) %*% S is inverted (`.rpemFisherCov`); the theta block (typical
value + covariate coefs + residual sd) is installed as the fit's `cov` and the
parFixedDf SE/%RSE/CI overwritten (`.rpemInstallFisherCov`, mirrors SAEM's
`.saemInstallFullCov`) with covMethod "fisher". Matches FOCEI SEs closely (tka SE
0.097 vs 0.097; add.sd SE 0.0028 vs 0.0030) and shrinks ~1/sqrt(n). KEY gotcha: the
eval-only FOCEI reports the ui iniDf theta (maxOuter=0 holds the starting value), so
`.rpemBuildFit` now SEEDS iniDf$est with the RPEM estimates before the eval -- this
also fixed a latent non-mixture estimate-display bug (parFixedDf showed ini values).
Multiple diagonal etas are supported via `rpemFisherDiag` (one typical-value and one
`om.<eta>` variance score per random effect); nEta==1 with covariates uses
`rpemFisherReg` (regression).  The omega-variance SEs are reported in the fit `$cov`
(rows/cols `om.<eta>`, on the variance scale -- the score is already in those units so
no delta method), matching SAEM.  The empirical Fisher (score outer product) is
asymptotically exact (E[EBE^2]/om^2 == I_1), so SEs converge to FOCEI's with enough
subjects (verified at n=100; small n like theo's 12 gives noisier SEs).  Other
structures (combined/tbs/power residual, non-mu-ref structural betas, mixtures) keep
the FOCEI-covariance (`covMethod="r,s"`) SEs.  Follow-ups: Fisher scores for
combined/power residual and structural betas; mixtures.

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
