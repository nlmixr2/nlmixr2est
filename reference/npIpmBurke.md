# Burke interior-point weight solver (nonparametric maximum likelihood)

Solves the convex nonparametric-maximum-likelihood weight problem for a
fixed set of support points: given the likelihood matrix `psi` (subjects
in rows, support points in columns) it returns the maximum-likelihood
mixing weights and the objective (log-likelihood). Exposed for testing
the C++ interior-point routine against golden fixtures.

## Usage

``` r
npIpmBurke(psi)
```

## Arguments

- psi:

  Numeric matrix, `psi[i, k] = p(y_i | support point k)`, with subjects
  in rows and support points in columns.

## Value

A list with `weights` (length `ncol(psi)`, non-negative, summing to 1)
and `objective` (the maximized log-likelihood).
