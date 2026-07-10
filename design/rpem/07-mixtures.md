# RPEM -- Mixtures (K > 1)

Deferred past the M1 breadth-first core, but the engine and control surface must
not hard-code K=1. This spec covers the K>1 design and the benchmark-and-choose
plan (D8).

STATUS (2026-07-10): K=2 single-`mix()` split-ETA mixtures work end-to-end via the
nlmixr2 `mix()` interface (E-step `rpemEstepMixDraw`, joint (i,j,k) MH M-step
`rpemMstepMix`; full `nlmixr2FitData` with mixProbabilities/mixNum/mixList).  M2 guards
DONE: (1) label switching -- components are exchangeable, so `.rpemFitMix` enforces an
ascending-mu canonical order every iteration (`order(muK)`), keeping the collected trace
coherent (verified: a reversed init still returns components in truth order); (2) empty /
collapsing components -- a weight of exactly 0 makes `logw = -Inf` next iteration and
permanently kills the component, so `rpemMstepMix` floors each weight at 1e-3 and
renormalizes so every component stays reachable (verified finite on a single population
fit as K=2).  Recovers unequal weights (0.75/0.25) and mu/omega/add.sd.

PER-COMPONENT Sigma^(k) DONE: the split-ETA-per-component form
`mix(exp(tka1 + eta.ka1), p1, exp(tka2 + eta.ka2))` -- each component with its own eta
and its own omega -- is now supported.  `.rpemMixInfo` parses the per-component etas and
sets `perComp`/`etaForComp` (component k -> its eta index); the classifier allows
nEta == number of distinct mixed etas; `.rpemFitMix` draws each component's eta from its
own omega (diag draw), and `rpemMstepMix` (C++) takes `etaForComp` and returns a per-eta
omega vector (components sharing an eta pool into one Sigma, so shared-eta stays a
special case).  Label-switching reorder is applied only for the shared-eta (exchangeable)
form; split-ETA components are tied to distinct symbols.  On split-ETA data SAEM collapses
the component means (tka1/tka2 ~ 0.45/0.67 for a true 0.0/1.4) and reports a single merged
omega, whereas RPEM separates the components (~0.11/1.71) and reports one omega per
component -- so RPEM meets/exceeds SAEM's mixture values (the target bar).  The
high-variance component's omega is still under-estimated (the known multi-eta E-step
under-coverage), but SAEM does not separate per-component BSV at all.  test-rpem-mix-percomp.R.
Still open: >1 mix() call, and the SAEM-representation backend benchmark.  Mixture
component-probability SEs are intentionally not reported (not reported elsewhere in the
ecosystem; mixtures keep covMethod "r,s").

## Two backends to build and compare

RPEM's identity is the joint `(i, k, theta)` MH sampling of discrete labels and
continuous params. There are two ways to declare/store the mixture; build both,
benchmark for correctness, and make the more accurate/robust one the default.

1. **RPEM-specific backend** -- mixture declaration and bookkeeping tuned to
   RPEM's M-step: components carry their own `mu^(k)`, `Sigma^(k)`, `w^(k)`;
   the MH proposes `k'` uniformly and accepts by Eq 32; weights update by Eq 27.
   Maximum fidelity to the paper; independent of SAEM's constraints.
2. **SAEM/NONMEM-representation backend** -- reuse SAEM's existing mixture
   plumbing (`nMix`, subpopulation omega sharing via `omegaShareSubpop`, `mu2`
   referencing in `R/mu2.R`) for how the user declares mixtures and how
   components are stored. RPEM swaps in its MH kernel but consumes the SAEM
   representation. Consistent UX with SAEM; less new UI surface.

## Benchmark-and-choose protocol

- Correctness targets: the paper's K=2 analytic model (Table 1) and the K=2
  parameterization referenced in the paper's ref [4]. Recover true `mu^(k)`,
  `Sigma^(k)`, `w^(k)` within the paper's tolerances.
- Compare the two backends on: parameter bias/coverage, label-switching
  stability, convergence iterations, and wall-clock.
- Pick the winner as the documented default; keep the other reachable via
  `rpemControl(mixBackend=)` if it is materially better on some model class.

## Mixture-specific concerns

- **Label switching**: enforce an identifiability convention (e.g. order
  components by a chosen `mu` coordinate) when reporting, as SAEM does.
- **Empty/collapsing components**: guard `w^(k) -> 0` and singular `Sigma^(k)`;
  reuse `nearPD` for covariance repair.
- **Weights**: Eq 27 update; ensure `sum_k w^(k) = 1` each iteration.
- **Predictions**: population prediction mixes over components (Eq 52); the fit
  object must store all components (`09-fit-object.md`).

## Open items

- OI-1: How mixtures interact with non-mu-ref fixed effects -- fixed effects are
  shared across components (not per-`k`) unless the model says otherwise.
  Confirm the sharing rule matches SAEM.
- OI-2: Decide the `rpemControl` grammar for specifying K and per-component
  starting values; prefer reusing SAEM's grammar for continuity.
