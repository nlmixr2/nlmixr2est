# Nudge exactly-zero theta initial estimates off zero (FOCEi family)

FOCEi scales a linear parameter by its native magnitude \`\|init\|\`,
which is \`0\` (no scale) when a population parameter is initialized at
exactly \`0\` (the parameter then freezes). This hook moves every
estimated, non-fixed \`theta\` whose initial estimate is exactly \`0\`
to \`+zeroTheta\` when that is within its bounds, otherwise
\`-zeroTheta\`; if neither is within the bounds it errors. Only runs for
the FOCEi family (a \`foceiControl\`), and runs before
\`.preProcessBoundedTransform\` so the nudged value is what gets
transformed.

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
