# Analytic outer gradient at a caller-supplied theta / eta / omega

The same C++ core the fit's own gradient uses (\`gradPooledCore\`), but
with the point passed in rather than read out of \`op_focei\`.
\`est="vae"\` with \`nonMuTheta="grad"\` needs exactly this: its M-step
evaluates the gradient at a theta and an encoder eta matrix that are not
the inner problem's, and at an omega that changes every step. Requires
\`foceiGradPooledSetupLoad\_()\` first, and a live FOCEi inner problem
(the shared pool, \`rxVaeOuter\`, and the theta/eta par_ptr maps all
come from it).

## Usage

``` r
foceiGradPooledDirect_(thVals, ebes, Oi, dOiEst, tr28, cores)
```

## Arguments

- thVals:

  natural-scale theta, in ntheta order

- ebes:

  nsub x neta matrix of etas to take the gradient at

- Oi:

  inverse of the current Omega

- dOiEst:

  list of d(Omega^-1)/d(estimation-scale omega element)

- tr28:

  the matching trace terms

- cores:

  thread count

## Value

natural-scale gradient (thetas, sigmas, omegas), or NULL if it declined
