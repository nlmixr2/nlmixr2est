# Install the pooled analytic-gradient setup for a non-focei caller

\`est="vae"\` with \`nonMuTheta="grad"\` evaluates the analytic outer
gradient once per M-step, at its own theta/eta/omega, and has no fit env
to hang the setup on. This installs the setup once so
\`foceiGradPooledDirect\_()\` can be called repeatedly.

## Usage

``` r
foceiGradPooledSetupLoad_(st)
```

## Arguments

- st:

  setup list from \`.foceiGradPooledSetup()\`

## Value

TRUE if the setup describes a shape the C++ gradient handles
