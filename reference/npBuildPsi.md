# Build the nonparametric Psi (conditional-likelihood) matrix

For an already set-up FOCEi inner problem (`vaeInnerSetup_`), evaluates
`psi[i, k] = p(y_i | support point k)` for each subject `i` (rows) and
support point `k` (columns), where each support point is an eta vector.
Exposed for testing the conditional-likelihood primitive.

## Usage

``` r
npBuildPsi(etaPoints, cores)
```

## Arguments

- etaPoints:

  Numeric matrix of support points, one per row (columns are etas).

- cores:

  Number of OpenMP threads.

## Value

Numeric matrix psi (subjects in rows, support points in columns).
