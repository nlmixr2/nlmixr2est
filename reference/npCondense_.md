# Condense support points (nonparametric engines)

Condense support points (nonparametric engines)

## Usage

``` r
npCondense_(lambda, psi, ratio = 0.001, tol = 1e-08)
```

## Arguments

- lambda:

  Support-point weights.

- psi:

  Conditional-likelihood matrix (subjects x support points).

- ratio:

  Weight-threshold ratio (keep weight \> max\*ratio).

- tol:

  QR rank-revealing tolerance.

## Value

List with 1-based kept indices from the weight threshold (`weightKeep`)
and from the subsequent QR pass (`qrKeep`).
