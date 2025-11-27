# C++ implementation of Matrix's nearPD

With \`ensureSymmetry\` it makes sure it is symmetric by applying
0.5\*(t(x) + x) before using nmNearPD

## Usage

``` r
nmNearPD(
  x,
  keepDiag = FALSE,
  do2eigen = TRUE,
  doDykstra = TRUE,
  only.values = FALSE,
  ensureSymmetry = !isSymmetric(x),
  eig.tol = 1e-06,
  conv.tol = 1e-07,
  posd.tol = 1e-08,
  maxit = 100L,
  trace = FALSE
)
```

## Arguments

- x:

  numeric \\n \times n\\ approximately positive definite matrix,
  typically an approximation to a correlation or covariance matrix. If
  `x` is not symmetric (and `ensureSymmetry` is not false),
  [`symmpart`](https://rdrr.io/pkg/Matrix/man/symmpart-methods.html)`(x)`
  is used.

- keepDiag:

  logical, generalizing `corr`: if `TRUE`, the resulting matrix should
  have the same diagonal
  ([`diag`](https://rdrr.io/r/base/diag.html)`(x)`) as the input matrix.

- do2eigen:

  logical indicating if a
  [`posdefify()`](https://rdrr.io/pkg/sfsmisc/man/posdefify.html) eigen
  step should be applied to the result of the Higham algorithm.

- doDykstra:

  logical indicating if Dykstra's correction should be used; true by
  default. If false, the algorithm is basically the direct fixpoint
  iteration \\Y_k = P_U(P_S(Y\_{k-1}))\\.

- only.values:

  logical; if `TRUE`, the result is just the vector of eigenvalues of
  the approximating matrix.

- ensureSymmetry:

  logical; by default,
  [`symmpart`](https://rdrr.io/pkg/Matrix/man/symmpart-methods.html)`(x)`
  is used whenever `isSymmetric(x)` is not true. The user can explicitly
  set this to `TRUE` or `FALSE`, saving the symmetry test. *Beware*
  however that setting it `FALSE` for an **a**symmetric input `x`, is
  typically nonsense!

- eig.tol:

  defines relative positiveness of eigenvalues compared to largest one,
  \\\lambda_1\\. Eigenvalues \\\lambda_k\\ are treated as if zero when
  \\\lambda_k / \lambda_1 \le eig.tol\\.

- conv.tol:

  convergence tolerance for Higham algorithm.

- posd.tol:

  tolerance for enforcing positive definiteness (in the final
  `posdefify` step when `do2eigen` is `TRUE`).

- maxit:

  maximum number of iterations allowed.

- trace:

  logical or integer specifying if convergence monitoring should be
  traced.

## Value

unlike the matrix package, this simply returns the nearest positive
definite matrix

## Details

This implements the algorithm of Higham (2002), and then (if `do2eigen`
is true) forces positive definiteness using code from
[`posdefify`](https://rdrr.io/pkg/sfsmisc/man/posdefify.html). The
algorithm of Knol and ten Berge (1989) (not implemented here) is more
general in that it allows constraints to (1) fix some rows (and columns)
of the matrix and (2) force the smallest eigenvalue to have a certain
value.

Note that setting `corr = TRUE` just sets `diag(.) <- 1` within the
algorithm.

Higham (2002) uses Dykstra's correction, but the version by Jens
Oehlschlägel did not use it (accidentally), and still gave reasonable
results; this simplification, now only used if `doDykstra = FALSE`, was
active in `nearPD()` up to Matrix version 0.999375-40.

## References

Cheng, Sheung Hun and Higham, Nick (1998) A Modified Cholesky Algorithm
Based on a Symmetric Indefinite Factorization; *SIAM J. Matrix Anal.\\
Appl.*, **19**, 1097–1110.

Knol DL, ten Berge JMF (1989) Least-squares approximation of an improper
correlation matrix by a proper one. *Psychometrika* **54**, 53–61.

Higham, Nick (2002) Computing the nearest correlation matrix - a problem
from finance; *IMA Journal of Numerical Analysis* **22**, 329–343.

## See also

A first version of this (with non-optional `corr=TRUE`) has been
available as
[`nearcor()`](https://rdrr.io/pkg/sfsmisc/man/nearcor.html); and more
simple versions with a similar purpose
[`posdefify()`](https://rdrr.io/pkg/sfsmisc/man/posdefify.html), both
from package sfsmisc.

## Author

Jens Oehlschlägel donated a first version. Subsequent changes by the
Matrix package authors.

## Examples

``` r
set.seed(27)
m <- matrix(round(rnorm(25),2), 5, 5)
m <- m + t(m)
diag(m) <- pmax(0, diag(m)) + 1
(m <- round(cov2cor(m), 2))
#>       [,1]  [,2]  [,3]  [,4]  [,5]
#> [1,]  1.00  0.65 -0.46 -1.15 -0.76
#> [2,]  0.65  1.00  0.58  0.50 -0.90
#> [3,] -0.46  0.58  1.00 -0.45 -0.32
#> [4,] -1.15  0.50 -0.45  1.00  0.25
#> [5,] -0.76 -0.90 -0.32  0.25  1.00

near.m <- nmNearPD(m)
round(near.m, 2)
#>       [,1]  [,2]  [,3]  [,4]  [,5]
#> [1,]  1.31  0.41 -0.24 -0.85 -0.75
#> [2,]  0.41  1.19  0.41  0.27 -0.91
#> [3,] -0.24  0.41  1.15 -0.24 -0.32
#> [4,] -0.85  0.27 -0.24  1.28  0.26
#> [5,] -0.75 -0.91 -0.32  0.26  1.00
norm(m - near.m) # 1.102 / 1.08
#> [1] 1.079735

round(nmNearPD(m, only.values=TRUE), 9)
#>             [,1] [,2] [,3] [,4] [,5]
#> [1,] 2.800681404    0    0    0    0
#> [2,] 1.831722441    0    0    0    0
#> [3,] 1.229003616    0    0    0    0
#> [4,] 0.076994641    0    0    0    0
#> [5,] 0.000000028    0    0    0    0

## A longer example, extended from Jens' original,
## showing the effects of some of the options:

pr <- matrix(c(1,     0.477, 0.644, 0.478, 0.651, 0.826,
               0.477, 1,     0.516, 0.233, 0.682, 0.75,
               0.644, 0.516, 1,     0.599, 0.581, 0.742,
               0.478, 0.233, 0.599, 1,     0.741, 0.8,
               0.651, 0.682, 0.581, 0.741, 1,     0.798,
               0.826, 0.75,  0.742, 0.8,   0.798, 1),
               nrow = 6, ncol = 6)

nc  <- nmNearPD(pr)
```
