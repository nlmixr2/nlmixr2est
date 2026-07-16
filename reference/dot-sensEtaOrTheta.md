# Calculate d(state)/d(eta) or d(state)/d(theta) sensitivities

Calculate d(state)/d(eta) or d(state)/d(theta) sensitivities

## Usage

``` r
.sensEtaOrTheta(s, theta = FALSE)
```

## Arguments

- s:

  symengine environment (from \`.loadSymengine()\`)

- theta:

  when \`TRUE\` calculate the sensitivities with respect to
  \`THETA\[#\]\`; otherwise with respect to \`ETA\[#\]\`

## Value

the symengine environment \`s\` augmented with the sensitivity equations
(\`..sens\`, \`..ddt\`, \`..stateInfo\`, ...)

## Author

Matthew L. Fidler
