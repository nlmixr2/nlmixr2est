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

- Use rxode2's thread-safe **threefry** RNG. NEVER use R's RNG
  (`norm_rand`/`unif_rand`/`GetRNGstate`) in the engine -- it is not thread-safe
  and would corrupt state under the OpenMP solve loop.
- **E-step draws (implemented):** draw all etas for the iteration UP FRONT via
  the already-registered `_rxode2_rxRmvnSEXP_` (threefry, thread-safe, used in
  `censResid.h`), then let the parallel solve loop read pre-drawn etas. Result:
  no RNG call inside any OpenMP region -- thread-safe by construction, and
  CRAN-compliant because rxRmvn is exported through the `rxode2ptr.h`
  function-pointer struct.
- **M-step RNG (M1): pre-drawn pools, no rxode2 change.** Within one M-step the
  population `mu`/`Omega` are fixed, so the MH proposals (`theta' ~ N(mu,Omega)`)
  and acceptance uniforms are iid, and for K=1 the chain is a single sequential
  walk. Draw all of them UP FRONT with the already-exported threefry `rxRmvn`
  (uniforms via `pnorm` of standard normals -- a deterministic, thread-safe math
  transform), then consume the pool in the chain. Fully thread-safe and
  CRAN-compliant with no scalar-RNG export.
- **M-step scalar RNG (deferred to M2):** parallel mixture chains (K>1) want
  per-chain scalar draws (`rxunif`/`rxnorm`/`rinorm`) that are NOT in
  `rxode2ptr.h`. Since CRAN forbids linking rxode2's C symbols directly, at M2
  export them through the registered function-pointer struct (rxode2 worktree:
  register `rxnorm`/`rinorm`/`rxunif`/`riunif`/`_setThreadInd`, append at the end
  of the `_rxode2_rxode2Ptr` list and `iniRxodePtrs0` to preserve indices). Note
  the cross-package build implication: nlmixr2est must compile against the
  modified header, so the dev rxode2 must be reachable at build time.
- Seed from a `seed`/`sim` control via rxode2's seed mechanism. Bitwise
  reproducibility holds for a fixed thread count; across thread counts, expect
  statistically-equivalent draws (rxode2's documented policy).

## Thread count

- Honor `rxode2::setRxThreads()` and the CRAN/CI caps in `tests/testthat.R`
  (cap within-solve threads to 2 on CRAN). Add an `rpemControl(cores=)` that
  defers to the rxode2 thread setting by default.

## Open items

- OI-1: Confirm the E-step batched-solve input layout that maximizes `par_solve`
  load balance (order rows by subject or by component?).
- OI-2: Measure whether the M-step is ever a parallel bottleneck; if it stays
  cheap vs the E-step, keep its parallelism minimal to reduce complexity.
