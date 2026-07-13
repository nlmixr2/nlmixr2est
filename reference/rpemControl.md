# RPEM estimation control

Options for the Randomized Parametric Expectation Maximization (RPEM)
estimation method (Chen et al. 2024). This is the K=1 minimal core; see
\`design/rpem/\` for the full roadmap.

## Usage

``` r
rpemControl(
  nGauss = 1000L,
  nMH = 50000L,
  mhBurn = 5000L,
  niter = 50L,
  collect = 15L,
  seed = 42L,
  atol = 1e-08,
  rtol = 1e-08,
  cores = 1L,
  impInflate = 0,
  likLbfgs = TRUE,
  lbfgsLmm = 5L,
  lbfgsFactr = 1e+07,
  lbfgsPgtol = 0,
  lbfgsMaxIter = 20L,
  print = 1L,
  printNcol = NULL,
  useColor = NULL,
  ...
)
```

## Arguments

- nGauss:

  Number of Monte Carlo samples per subject in the E-step (\`m_Gauss\`).

- nMH:

  Number of Metropolis-Hastings trials collected in the M-step.

- mhBurn:

  Number of M-step MH burn-in trials (discarded).

- niter:

  Maximum number of E-M iterations.

- collect:

  Number of terminal iterations averaged for the final estimate.

- seed:

  RNG seed for the threefry sampler (reproducible for a fixed thread
  count).

- atol, rtol:

  ODE solver tolerances.

- cores:

  Number of cores for the threefry draw (solve threading is set by
  rxode2).

- impInflate:

  Opt-in mode-centered importance sampling for the E-step. \`0\`
  (default) keeps the paper's prior sampling (draw eta ~ N(0, Omega)). A
  value \`\>= 1\` draws instead from N(EBE, impInflate\*Omega) –
  centered at the previous iteration's posterior mean with that
  variance-inflation factor – and importance-weights, improving
  posterior-tail coverage for high-variance random effects in multi-eta
  models (whose largest Omega prior sampling under-estimates).
  Experimental: a partial mitigation, not a full fix (see
  design/rpem/04).

- likLbfgs:

  For a general log-likelihood (\`ll()\`) endpoint, refine the
  fixed-effect likelihood parameters each iteration by a box-constrained
  L-BFGS-B optimization of the importance-weighted observation
  log-likelihood (mirrors the saem/saemix ind.fix10 step), respecting
  the parameter bounds from the model, rather than the default single
  damped-Newton re-solve step. \`TRUE\` (default) for \`ll()\` models;
  ignored for standard residual-error models.

- lbfgsLmm, lbfgsFactr, lbfgsPgtol, lbfgsMaxIter:

  L-BFGS-B tuning for the \`likLbfgs\` likelihood-parameter refinement
  (number of corrections, convergence \`factr\`/\`pgtol\`, and max
  iterations).

- print:

  Iteration-print frequency: display the parameter walk (population
  estimates + omega, with the back-transformed row) every \`print\`
  iterations (saem/focei/vae style), streamed live as the loop runs.
  \`1\` (default) prints every iteration; \`0\` captures the parameter
  history silently. The walk is \*always\* saved to the fit object's
  parameter history (\`fit\$parHist\` / \`fit\$parHistStacked\`)
  regardless. May also be an \`iterPrintControl()\` object.

- printNcol, useColor:

  Iteration-print formatting (columns per row, ANSI color); passed
  through to \`iterPrintControl()\`.

- ...:

  Ignored (reserved for future options).

## Value

A list of class \`rpemControl\`.
