# RPEM -- Mixtures (K > 1)

Deferred past the M1 breadth-first core, but the engine and control surface must
not hard-code K=1. This spec covers the K>1 design and the benchmark-and-choose
plan (D8).

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
