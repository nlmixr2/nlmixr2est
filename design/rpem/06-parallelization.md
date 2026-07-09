# RPEM -- Parallelization

Strategy (D3): OpenMP within-process, plus rxode2's threaded population solve.
No MPI. Must stay CRAN-portable and honor nlmixr2's thread policy
(`tests/testthat.R`, `aaaCranNlmixrThreads.R`).

## Where the time goes

Per the paper, the E-step ODE solves dominate. So the primary parallel win is
already delivered by handing one big population solve to rxode2 `par_solve`,
which threads over (subject x sample x component). RPEM's own OpenMP is for the
non-ODE reductions and the M-step chain.

## Parallel regions

1. **E-step ODE batch**: single `par_solve` over the full sample matrix. RPEM
   does not add its own OpenMP here -- it lets rxode2 thread it. Build the input
   matrix contiguously so rxode2 balances the load.
2. **E-step likelihood reduction**: `n_ik`, `N_i`, `lnL` via log-sum-exp.
   Parallelize over subjects/components with an OpenMP reduction; each thread
   owns disjoint (i,k) so no races.
3. **M-step MH chains**: run K (or per-component) independent chains in parallel;
   within a chain the Markov steps are serial. Parallelize across independent
   chains / across subjects' contributions to the conjugate sums.
4. **Numeric fixed-effect Q**: the Q evaluation sums over accepted samples ->
   OpenMP reduction; the optimizer itself stays serial.

Respect `ARMA_DONT_USE_OPENMP` / `EIGEN_DONT_PARALLELIZE` (already set in
Makevars) -- do not nest Armadillo/Eigen threading inside RPEM's OpenMP regions.

## RNG and determinism (A2)

- Use rxode2's thread-safe parallel RNG (threefry-style) seeded from a
  `seed`/`sim` control, so results are reproducible for a fixed thread count.
- Each (thread, iteration, subject, component, sample) must map to a
  deterministic RNG substream so the same seed gives the same draws regardless
  of scheduling. This is the standard rxode2 pattern -- reuse it, do not roll a
  new generator.
- Document that bitwise reproducibility holds for a fixed number of threads;
  across thread counts, expect statistically-equivalent but not identical draws
  (same policy rxode2 already documents).

## Thread count

- Honor `rxode2::setRxThreads()` and the CRAN/CI caps in `tests/testthat.R`
  (cap within-solve threads to 2 on CRAN). Add an `rpemControl(cores=)` that
  defers to the rxode2 thread setting by default.

## Open items

- OI-1: Confirm the E-step batched-solve input layout that maximizes `par_solve`
  load balance (order rows by subject or by component?).
- OI-2: Measure whether the M-step is ever a parallel bottleneck; if it stays
  cheap vs the E-step, keep its parallelism minimal to reduce complexity.
