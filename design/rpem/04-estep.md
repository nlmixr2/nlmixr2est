# RPEM -- E-step

Purpose: for each subject `i` and mixture component `k`, estimate
`n_ik = integral p(Y_i|theta_i) p(theta_i|mu^(k),Sigma^(k)) dtheta_i` (Eq 23) by
Monte Carlo, and cache the per-sample likelihoods for the M-step.

Implemented entirely in C++ (`src/rpem.cpp`, D17) -- no per-iteration R
round-trip. The steps below run inside the C++ iteration loop; the batched solve
uses rxode2's C-level `par_solve` in-process and reads `rx_pred_` straight from
the solve buffer.

## Procedure (per iteration)

1. For each `k`, factor `Sigma^(k) = L_k L_k^T` (Cholesky). Draw `m_Gauss`
   standard-normal vectors per subject; form `theta = mu^(k) + L_k z` on the
   transformed scale (Eq 24 sampling). Redraw every iteration (D7).
2. Back-transform `theta` to natural scale; assemble the population solve input
   for all (subject x sample x component) rows.
3. **Single batched `par_solve`** over the whole matrix (see
   `06-parallelization.md`), called from C++ via rxode2's C API. One solve call
   per iteration, not per sample, and no R round-trip.
4. Evaluate `p(Y_i | theta)` per sample from the **dedicated RPEM likelihood
   model** (`13-likelihood-model.md`): supply `THETA` (population/fixed) and
   `ETA` (the drawn eta) as inputs, no optimization. Read `rx_pred_` from the
   solve buffer; `log p(Y_i | theta) = -sum_obs(rx_pred_)` per (i, k, sample).
   Exponentiate/stabilize (log-sum-exp) to get the contribution to `n_ik`.
5. `n_ik = mean over samples of p(Y_i | theta)` (Eq 24). Use a log-sum-exp with
   a per-subject max for numerical stability, then divide by `m_Gauss`.
6. `N_i = sum_k w^(k) n_ik` (Eq 22); `tau_i(k) = w^(k) n_ik / N_i` (Eq 25).
7. `lnL = sum_i ln(N_i)` (Eq 26).

For K=1 (M1), `w^(1)=1`, `N_i = n_i1`, and the per-component loop collapses.

## Storage for M-step reuse

The M-step (Eq 33 MH) draws proposals from the *same* per-iteration sample pool,
so it must reuse E-step results without new ODE solves. Cache per iteration:

- `theta` samples: `n * K * m_Gauss * nRandom` doubles (transformed scale).
- per-sample log-likelihood `log p(Y_i | theta)`: `n * K * m_Gauss` doubles.
- per-sample residual-likelihood pieces needed to re-score `beta` proposals in
  the numeric M-step without re-solving (see `05-mstep.md`).

Memory guard: with `n=50, K=1, m_Gauss=1000, nRandom=7` this is ~3M doubles
(~24 MB) -- fine. Expose `m_Gauss` so users can trade memory/variance. For very
large `n*K*m_Gauss`, allow a chunked mode (solve+score in blocks, keep only what
the M-step needs). Flag if the store must be disk-backed (not expected M1).

## Controls (see `10-stopping-control.md`)

- `nGauss` (`m_Gauss`): samples per subject per component. Default 1000
  (paper uses 200-3000; 600-1000 typical).
- Reuse policy: fresh draws each iteration (fixed by D7).

## Reuse points

- rxode2 `par_solve` for the ODE batch.
- The dedicated RPEM likelihood model (`13-likelihood-model.md`) for
  `p(Y_i | theta)` -- a hybrid of `nlmModel0` (no-optimization eval) and
  `foceiModel0ll` (THETA+ETA). This is what makes multiple endpoints / error
  structures / autocorrelation / censoring reusable.

## Open items

- OI-1: RESOLVED (D2) -- do NOT reuse the FOCEI inner path or
  `foceiControl(maxInnerOpt=0)`. Build the dedicated no-optimization phi/eta
  likelihood model; see `13-likelihood-model.md`. `ETA` is confirmed suppliable
  in the solve path because this model emits a pred/r output (`rx_pred_`/`rx_r_`)
  rather than a DV prediction.
- OI-2: Whether censoring likelihood (M2/M3/M4) flows through the same model
  (it should, via the `rx_pred_f_`/`rx_r_` lines and `censResid`/`censEst`),
  which would pull censoring earlier than the roadmap assumes. Keep deferred
  unless free.
