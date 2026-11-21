# Register a likelihood penalty from an rxode2 user-function expansion

Called from an \`rxUdfUi\` method while it expands a model term (for
example a smoothing spline) against the data. The method must return
\`uiUseData = TRUE\` from its data-free stage, so \`nlmixr2()\`
re-parses the model with the estimation data
(\`rxode2::rxUdfUiData()\`); the specs registered during that re-parse
are attached to the expanded model with \[nlmixr2LikPenaltyAdd()\], so
the penalty is active only while that model is estimated. Outside of an
nlmixr2 data expansion (for example in \`rxSolve()\`) the call does
nothing. See the "A model term that brings its own penalty" section of
\`vignette("likelihood-penalties")\` and the 'nlmixr2spline' package for
a complete example.

## Usage

``` r
nlmixr2LikPenaltyUdf(spec)
```

## Arguments

- spec:

  a \[nlmixr2LikPenaltySpec()\] with a unique \`name\`; a later
  registration with the same name replaces it

## Value

\`TRUE\` when the penalty was registered, invisibly

## Author

Matthew L. Fidler
