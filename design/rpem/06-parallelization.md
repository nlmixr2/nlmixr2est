# RPEM -- Parallelization

Strategy (D3): OpenMP within-process, plus rxode2's threaded population solve.
No MPI. Must stay CRAN-portable and honor nlmixr2's thread policy
(`tests/testthat.R`, `aaaCranNlmixrThreads.R`).

## Where the time goes

Per the paper, the E-step ODE solves dominate. So the primary parallel win is
handing the solve to rxode2 `par_solve`, which owns the per-thread ODE workspace.

STATUS (DONE, commit 3349a2b5): the E-step is parallel. With `cores>1`, for each
Monte Carlo sample the engine sets every subject's parameters (population THETA +
that sample's eta draw) and solves all subjects in ONE `par_solve` call -- the
thread-safe route (a manual OpenMP loop over `ind_solve` segfaults because that
per-thread ODE scratch is not set up). The cheap prediction read (`calc_lhs`, no
ODE work) stays serial; `cores==1` keeps the serial `ind_solve` reference path.
~2x speedup on 4 cores. The eta sampling was moved from C++ into R (`.drawEtas`
via `rxode2::rxRmvn`), drawn before the solve, so the sampling RNG is decoupled
from the solve core count: a fixed core count is exactly reproducible; across core
counts fits agree to Monte-Carlo precision (the MH is a Markov chain, so one
accept/reject flip from the ~1e-8 `ind_solve`-vs-`par_solve` difference diverges
the chain to a statistically-equivalent trajectory). Verified in
`test-rpem-parallel.R`. Not yet done: batching the whole (subject x sample) grid
into ONE `par_solve` (would parallelize across samples too, not just subjects per
sample); the M-step reductions remain serial.

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
- **Per-thread scalar RNG -- AVAILABLE (D21).** rxode2 now provides
  `seedEng(ncores)` and `rxNormEng(mean, sd)` in the registered `rxode2ptr.h`
  struct (indices 72-73). Pattern: `seedEng(cores)` once, then inside the OpenMP
  loop each thread calls `setRxThreadId(omp_get_thread_num())` and draws with
  `rxNormEng(0, 1)` on its own threefry engine -- thread-safe, CRAN-compliant, no
  up-front pool needed. Use this for the OpenMP E-step (task #7) and the M2
  parallel mixture chains. The M1 rxRmvn up-front-pool approach (D19) still works
  and remains fine for the single-chain M1 M-step.
- **Seeding convention (important):** seed the per-subject RNG stream as
  `seed0 + id*2` (per `id`, stride 2), so RPEM's draws are INDEPENDENT of the
  RNG streams rxode2 consumes during ODE solving. Using the plain `seed0` or a
  stride of 1 risks RPEM's draws colliding/correlating with the solver's,
  biasing the Monte Carlo. Apply this wherever RPEM seeds its own draws
  (E-step eta sampling, M-step) once they use the per-thread `rxNormEng` path.
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
