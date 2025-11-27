# Generalized Cholesky Matrix Decomposition

Performs a (modified) Cholesky factorization of the form

## Usage

``` r
cholSE(matrix, tol = (.Machine$double.eps)^(1/3))
```

## Arguments

- matrix:

  Matrix to be Factorized.

- tol:

  Tolerance; Algorithm suggests (.Machine\$double.eps) ^ (1 / 3),
  default

## Value

Generalized Cholesky decomposed matrix.

## Details

t(P) %\*% A %\*% P + E = t(R) %\*% R

As detailed in Schnabel/Eskow (1990)

## Note

This version does not pivot or return the E matrix

## References

matlab source:
http://www.dynare.org/dynare-matlab-m2html/matlab/chol_SE.html; Slightly
different return values

Robert B. Schnabel and Elizabeth Eskow. 1990. "A New Modified Cholesky
Factorization," SIAM Journal of Scientific Statistical Computing, 11, 6:
1136-58.

Elizabeth Eskow and Robert B. Schnabel 1991. "Algorithm 695 - Software
for a New Modified Cholesky Factorization," ACM Transactions on
Mathematical Software, Vol 17, No 3: 306-312

## Author

Matthew L. Fidler (translation), Johannes Pfeifer, Robert B. Schnabel
and Elizabeth Eskow
