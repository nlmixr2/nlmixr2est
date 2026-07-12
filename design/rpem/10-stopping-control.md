# RPEM -- Stopping Rule and rpemControl()

## Convergence (paper "Stopping criterion")

- Track the log-likelihood `lnL = sum_i ln(N_i)` per iteration.
- Fit a least-squares line to the last `slopeWindow` (default 30) `lnL` values.
- Before convergence the slope is positive. Stop at the first iteration where
  the slope is no longer positive (`slope <= 0`), i.e. the trace has flattened.
- Guard with `maxIter` (hard cap) and a minimum iteration count so the slope
  test has enough points (>= `slopeWindow`).
- Because MC integrals fluctuate, expect the slope to hover near zero at
  convergence -- the "first non-positive slope" rule is the paper's; keep it but
  make the window and a small tolerance configurable.

## rpemControl() options (draft)

Inherits `sharedControl()`. RPEM-specific (names to finalize against SAEM/FOCEI
control naming conventions):

| Option | Meaning | Default |
|---|---|---|
| `nGauss` | E-step samples per subject per component (`m_Gauss`) | 1000 |
| `nMH` | M-step MH trials per component | `100 * n` |
| `mhBurn`, `mhThin` | MH burn-in / thinning | tbd |
| `slopeWindow` | iterations in the slope test | 30 |
| `slopeTol` | slope tolerance for "non-positive" | 0 |
| `maxIter` | hard iteration cap | tbd |
| `nMix` | number of mixture components K | 1 |
| `mixBackend` | `"rpem"` or `"saem"` mixture representation (see `07`) | tbd via benchmark |
| `seed` | RNG seed for reproducibility | tbd |
| `cores` | threads; defers to `setRxThreads()` | rxode2 setting |
| `seWindow` | terminal iterations collected for SEs | tbd |
| `stiff` | non-stiff vs stiff ODE solver selection | rxode2 default |
| `atol`, `rtol` | ODE tolerances | rxode2 defaults (paper uses 1e-6) |

## Iteration reporting

- Reuse `R/iterPrintControl.R` conventions so the RPEM progress line matches
  SAEM/FOCEI (iteration, lnL, slope, key params).

## Open items

- OI-1: Finalize option names to match existing control vocabulary (grep
  `saemControl`/`foceiControl` for the analogous knobs before naming).
- OI-2: Decide sensible `maxIter`, `nMH`, burn-in/thin defaults from the M1
  benchmark runs rather than guessing now.
