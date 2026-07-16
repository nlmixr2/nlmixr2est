# Run the NPAG adaptive-grid cycle on a set-up inner problem

Requires the FOCEi inner problem to be set up (`.npInnerSetup`). Runs
the full Yamada adaptive-grid cycle (Sobol grid -\> Psi -\> Burke IPM
-\> condensation -\> expansion -\> convergence) and returns the discrete
mixing distribution. Exposed for testing ahead of the full fit-object
wiring.

## Usage

``` r
npagCycle_(
  lower,
  upper,
  points = 2028L,
  cycles = 100L,
  cores = 1L,
  gammaOptimize = FALSE
)
```

## Arguments

- lower, upper:

  Numeric vectors, the per-eta support-point box.

- points:

  Initial Sobol grid size.

- cycles:

  Maximum cycles.

- cores:

  OpenMP threads.

- gammaOptimize:

  Optimize the residual-error magnitude (gamma) each cycle (only valid
  for uncensored normal endpoints).

## Value

A list with `support` (support points, eta space; one per row),
`weights`, `objf` (log-likelihood), `gamma`, `cycles`, and `converged`.
