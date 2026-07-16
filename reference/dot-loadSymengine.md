# Load a model into a symengine environment

Load a model into a symengine environment

## Usage

``` r
.loadSymengine(newmod, promoteLinSens = TRUE, fullModel = FALSE)
```

## Arguments

- newmod:

  model text (normalized rxode2 model, e.g. from a prune)

- promoteLinSens:

  when \`TRUE\`, promote \`linCmt()\` to the sensitivity-based solved
  system

- fullModel:

  when \`TRUE\`, change the printed message to indicate the full model
  is being loaded

## Value

symengine environment from \`rxode2::rxS()\` with \`rx_r\_\` coerced to
a symengine object when needed

## Author

Matthew L. Fidler
