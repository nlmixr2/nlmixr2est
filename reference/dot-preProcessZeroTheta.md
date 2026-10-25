# Handle theta initial estimates that are exactly zero

A linear \`theta\` has no magnitude to be scaled by when it starts at
\`0\`, so this hook records a typical magnitude for each estimated,
non-fixed, non-residual theta that starts at \`0\` (\`1/sd()\` of the
covariate it multiplies, else 1). The scaling setup uses it as that
parameter's \`scaleC\`. For the FOCEi family (a \`foceiControl\`) the
parameter is also moved to \`+zeroTheta\` when that is within its
bounds, otherwise \`-zeroTheta\`; if neither is within the bounds it
errors. Runs before \`.preProcessBoundedTransform\` so the moved value
is what gets transformed.

## Usage

``` r
.preProcessZeroTheta(ui, est, data, control)
```

## Arguments

- ui:

  rxode2 ui model

- est:

  estimation method (all methods are shown by \`nlmixr2AllEst()\`).
  Methods can be added for other tools

- data:

  nlmixr data

- control:

  The estimation control object. These are expected to be different for
  each type of estimation method

## Value

list with the ui (possibly modified)

## Author

Matthew L. Fidler
