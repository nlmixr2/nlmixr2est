# RPEM -- M-step

The M-step is where RPEM is novel: a single randomized Metropolis-Hastings chain
samples the discrete labels `(i, k)` and the continuous `theta_i` together, then
the population parameters are updated as sample averages of the target
distribution `pi(s) = g_ik(theta_i)/N` (Eq 30). No new ODE solves -- everything
reuses the E-step sample pool and stored likelihoods.

Implemented in C++ (`src/rpem.cpp`, D17) inside the iteration loop, with no R
round-trip. The MH chain, conjugate updates, and the numeric fixed-effect
re-scoring all operate on in-memory C++ buffers from the E-step.

## Joint Metropolis-Hastings (Eq 31-33)

State `s = {i, k, theta_i}`. Proposal `s' = {i', k', theta'_{i'}}`:

- draw `i'` uniform on subjects, `k'` uniform on components,
- draw `theta'_{i'}` from `N(mu^(k'), Sigma^(k'))` (reuse an E-step sample).

Acceptance ratio (Eq 32), with the symmetric proposal terms cancelling:

```
A = min(1, [p(Y_{i'}|theta'_{i'}) / p(Y_i|theta_i)] * [N_i / N_{i'}] * [w^(k') / w^(k)])
```

All three factors are already cached from the E-step. Run the chain to collect
`m` samples of `f_ik(theta_i)` (Eq 33); default trial count ~`100*n` per
component (paper's autocorrelation-informed default).

## Conjugate updates from the accepted samples (Eq 15-21)

Cast every population update as the generalized average Eq 28,
`<f> = (sum accepted f) / (count)`:

- `mu^(k)` (Eq 15/18): sample mean of `theta` (resp. `alpha`) over accepted
  samples labelled `k`.
- `Sigma^(k)` (Eq 16/19): sample covariance about the new `mu^(k)`.
- `w^(k)` (Eq 27): fraction of accepted samples labelled `k` (equivalently
  `(1/n) sum_i n_ik/N_i`).

These operate on the transformed-scale random params (class 1 in
`03-parameter-space.md`).

## Fixed effects -- hybrid update (D6)

Fixed effects (class 3: residual error, non-mu-ref structural/covariate) are the
paper's `beta`. Two paths, chosen per parameter:

- **Sigma / residual-error params -- reuse SAEM's error-model updates (D20).**
  SAEM already derives the sigma updates for every supported normal error type:
  additive, proportional, combined1/combined2 (add+prop), and the power/lambda
  (`propT`/Box-Cox-Yeo-Johnson) variants -- see `src/saem.cpp` (`rmProp`,
  `rmAddProp`, `rmAddPropLam`, `handleF`, the `ares`/`bres` machinery) and
  `saem_fit_aux.R`. RPEM's M-step reuses these SAME formulas, applied to the
  accepted MH samples (the sums over accepted `(i,j)` replace SAEM's stochastic-
  approximation sums). For pure additive this is the paper's Eq-17 closed form
  (`new addSd^2 = sum accepted SS / sum accepted nobs`), already implemented in
  `rpemMstepK1` by backing `SS_ij` out of the stored log p (exact per C1.1).
  Proportional/combined/power reuse SAEM's `ares`/`bres`/`yptp` update
  expressions over the accepted samples.
- **Numeric** for the residuals SAEM does NOT close-form, and for non-mu-ref
  structural coefficients: maximize the Monte-Carlo Q-function
  `Q(beta) = sum over accepted samples of log p(Y_i | theta_i, beta)`
  over `beta`, using the nlm-family optimizer. Because the accepted `theta`
  samples and their cached structural predictions (`rx_pred_f_`) are stored,
  re-scoring a candidate `beta` re-evaluates only the residual-likelihood
  (`rx_r_`) branch of the RPEM likelihood model (`13-likelihood-model.md`) and
  does NOT re-solve ODEs. This is what makes the numeric path cheap per
  iteration.

The classifier decides the path; a model can have some closed-form and some
numeric fixed effects simultaneously.

## Correctness ordering within an iteration

Update sequence: (1) `mu^(k)`, (2) `Sigma^(k)`, (3) `w^(k)`, (4) fixed effects.
Fixed effects last so their Q uses the freshly-updated population draws'
labels/weights. Confirm against the paper's ordering; the EM monotonicity
property (`L(phi^(r+1)) >= L(phi^(r))`) is the invariant to watch in tests.

## Reuse points

- Cached E-step likelihoods, `N_i`, `w^(k)` -> acceptance ratio.
- nlm-family optimizer -> numeric `beta` maximization.
- Residual-likelihood re-scoring -> RPEM likelihood model `rx_r_` branch, cached
  `rx_pred_f_` (no ODE re-solve; see `13-likelihood-model.md`).

## Open items

- OI-1: Confirm which residual shapes qualify for the Eq 17 conjugate path in
  nlmixr2's error-model taxonomy; default everything else to numeric.
- OI-2: Chain warm-start: seed the MH state from the previous iteration's
  accepted mean to cut burn-in (speed, not correctness). Verify it does not bias
  estimates before enabling by default.
- OI-3: Decide burn-in / thinning defaults for the MH chain; expose as controls.
