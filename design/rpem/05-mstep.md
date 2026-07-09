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
  approximation sums). DONE for additive AND proportional (both closed form) in
  `rpemMstepK1Reg` via an `errType` arg: the E-step (`rpemSolveSubject`)
  accumulates per sample the additive `SS = sum (DV-cp)^2` and the proportional
  `WSS = sum ((DV-cp)/cp)^2` from the cached `cp = rx_pred_f_`, and the update is
  `new sd = sqrt(sum accepted / sum accepted nobs)` (-> add.sd or prop.sd). Both
  verified vs FOCEI (`test-rpem-prop.R`); the additive path is algebraically the
  same as the earlier back-out-from-log-p (C1.1), so no-cov tests are unchanged.
  DONE: combined (add+prop) too. Its variance V = add^2 + prop^2 cp^2 has no
  closed form, so the E-step caches per-obs cp^2 (`rpemOp.cp2v`) and residual^2
  (`rpemOp.r2v`), and `rpemMstepK1Comb` maximizes the Gaussian log-likelihood
  sum_visited count * sum_obs [-0.5 log V - r^2/(2V)] over (a=add^2, b=prop^2) by a
  GUARDED optimizer: a Newton direction only when the local Hessian is
  negative-definite (V is not globally concave in (a,b) -- raw Newton diverges to
  1e39), else gradient ascent, with a backtracking line search that accepts a step
  only if it raises the objective and keeps a,b>0. Verified vs FOCEI
  (`test-rpem-comb.R`).
  TBS (transform-both-sides, Box-Cox / Yeo-Johnson) -- IN PROGRESS. Two hard
  blockers resolved and committed:
    1. Model generation: rxCombineErrorLines emitted each eta as an omega
       DECLARATION (`eta ~ omega`) rather than an `eta <- ETA[k]` assignment; with
       a transform present that undeclared-then-declared symbol crashed the
       symengine load (rxTBS parse error), which is why only untransformed models
       worked. rpemModel0 now strips the eta `~` declarations and injects
       `eta <- ETA[k]` assignments before first use.
    2. Jacobian: the generated llikNorm is the normal loglik of the TRANSFORMED
       data and omits the change-of-variables Jacobian log|dt/dDV| (FOCEI adds it
       in C). For a dynamic/estimated lambda the Jacobian depends on lambda, so
       rpemModel0 now adds `log|rxTBSd(DV,...)|` (rxTBSd = d t/d DV) before negating
       rx_pred_ when a boxCox/yeoJohnson transform is present. Verified the full
       Box-Cox model -sum(rx_pred_) matches the hand transformed-normal loglik with
       the Jacobian.
  DONE: TBS lambda AND power, both via a 1-D PROFILE optimize (profile out the
  scale in closed form, golden-section on the exponent/lambda) over the accepted
  MH samples. The shared `rpemMHReg` runs the joint MH + regression; per-obs raw
  cp/DV are cached (`rpemOp.cpv/dvv`, replacing the combined-only cp^2/r^2) so any
  candidate residual param can be re-scored.
    - TBS (errType 3): additive on the transformed scale + dynamic Box-Cox/
      Yeo-Johnson lambda. `rpemMstepK1TBS` maximizes
      -0.5 N log(SS(lambda)/N) + sum log|dt/dDV| with `_powerD`/`_powerDD` using
      the model's own yj/low/hi (`.rpemExtractTBS`), so the C++ transform matches
      the model exactly. Matches FOCEI on tka/add.sd/lambda (`test-rpem-tbs.R`).
    - Power (errType 4): variance (prop.sd*cp^power)^2, exponent estimated.
      `rpemMstepK1Pow` maximizes -0.5[N log(SSc/N) + 2c sum log|cp|]. Matches
      FOCEI on tka/prop.sd/power (`test-rpem-pow.R`; power location needs more MC
      precision -- late-time low-cp points dominate).
  REMAINING: combined on the transformed scale (add+prop+TBS); power/lambda (propT)
  closed forms could still reuse SAEM's `ares`/`bres`/`yptp` if faster.
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
