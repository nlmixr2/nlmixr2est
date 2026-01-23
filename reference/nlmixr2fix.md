# Try to fix a nlmixr2 fit

Currently this re-evaluates the function in the current version of
rxode2.

## Usage

``` r
nlmixr2fix(fit)
```

## Arguments

- fit:

  nlmixr2 fit object from a different version of nlmixr2.

## Value

A nlmixr2 fit that has been (possibly) adjusted to work with the current
version of nlmixr2.

## Author

Matthew L. Fidler

## Examples

``` r
# \donttest{

if (requireNamespace("qs", quietly = TRUE)) {

  # This is a nlmixr2 v3 fit
  fit <- system.file("testfit_nlmixr3.rds", package = "nlmixr2est")
  fit <- readRDS(fit)

  # While it prints well, it can't be used in all functions because
  # Language features (like +var()) are not supported in the v3 version

  print(fit)

  try(rxSolve(fit)) # should error, but with try it will just display the error

  # This function attempts to fix it by regenerating the rxode2 model with the
  # new features

  # This function also prints out the information on how this fit was created

  fit <- nlmixr2fix(fit)

  # Now solving and other functions work

  rxSolve(fit)

}
#> Warning: decompression of an rxUi object from rxode2 < 4.0 requires qs which is not on CRAN
#> Warning: decompression of an rxUi object from rxode2 < 4.0 requires qs which is not on CRAN
#> ── nlmixr² SAEM OBJF by FOCEi approximation ──
#> 
#>  Gaussian/Laplacian Likelihoods: AIC() or $objf etc. 
#>  FOCEi CWRES & Likelihoods: addCwres() 
#> 
#> ── Time (sec $time): ──
#> 
#>            setup covariance  saem table compress  other
#> elapsed 0.001395   0.010005 4.071 0.063     0.02 1.7026
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.     SE %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tka           Ka 0.465  0.195 41.8       1.59 (1.09, 2.33)     71.0    -0.403% 
#> tcl           Cl  1.01 0.0846  8.4       2.74 (2.32, 3.23)     27.3      4.51% 
#> tv             V  3.46 0.0455 1.32         31.7 (29, 34.6)     13.4      8.80% 
#> add.sd           0.695                               0.695                     
#>  
#>   Covariance Type ($covMethod): linFim
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 132 × 19
#>   ID     TIME    DV  PRED    RES IPRED   IRES  IWRES eta.ka eta.cl   eta.v    cp
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl>
#> 1 1      0     0.74  0     0.74   0     0.74   1.06  0.0869 -0.481 -0.0832  0   
#> 2 1      0.25  2.84  3.28 -0.437  3.83 -0.995 -1.43  0.0869 -0.481 -0.0832  3.83
#> 3 1      0.57  6.57  5.85  0.715  6.76 -0.194 -0.280 0.0869 -0.481 -0.0832  6.76
#> # ℹ 129 more rows
#> # ℹ 7 more variables: depot <dbl>, center <dbl>, ka <dbl>, cl <dbl>, v <dbl>,
#> #   tad <dbl>, dosenum <dbl>
#> Warning: decompression of an rxUi object from rxode2 < 4.0 requires qs which is not on CRAN
#> Warning: decompression of an rxUi object from rxode2 < 4.0 requires qs which is not on CRAN
#> Warning: decompression of an rxUi object from rxode2 < 4.0 requires qs which is not on CRAN
#> Error in if (pred1$variance) { : argument is of length zero
#> # This function is meant to load nlmixr2 fits from other versions
#> # To reproduce the fit, you need to use the same version of nlmixr2
#> ## ==============================
#> ## nlmixr2est Session Information
#> ## ==============================
#> # OS: Debian GNU/Linux 12 (bookworm)
#> # LAPACK: /cosbi/home/boaretti/miniconda3/envs/r_env3/lib/libopenblasp-r0.3.29.so
#> # LAPACK Version: 3.12.0
#> # R Version: R version 4.3.3 (2024-02-29)
#> 
#> 
#> # Install Package version 1.3.1-13 of 'dparser'
#> remotes::install_version("dparser", version = "1.3.1-13") 
#> 
#> # Install Package version 1.0.0 of 'lotri'
#> remotes::install_version("lotri", version = "1.0.0") 
#> 
#> # Install Package version 0.7 of 'PreciseSums'
#> remotes::install_version("PreciseSums", version = "0.7") 
#> 
#> # Install Package version 2.0.11 of 'rxode2ll'
#> remotes::install_version("rxode2ll", version = "2.0.11") 
#> 
#> # Install Package version 3.0.4 of 'rxode2'
#> remotes::install_version("rxode2", version = "3.0.4") 
#> 
#> # Install Package version 2024-3.5 of 'lbfgsb3c'
#> remotes::install_version("lbfgsb3c", version = "2024-3.5") 
#> 
#> # Install Package version 6.0.1-12 of 'n1qn1'
#> remotes::install_version("n1qn1", version = "6.0.1-12") 
#> 
#> # Install Package version 3.0.4 of 'nlmixr2est'
#> remotes::install_version("nlmixr2est", version = "3.0.4") 
#> 
#> # Install Package version 3.0.0 of 'nlmixr2extra'
#> remotes::install_version("nlmixr2extra", version = "3.0.0") 
#> 
#> # Package 'nlmixr2lib' is not installed, but known to enhance nlmixr2/babelmixr2
#> 
#> # Install Package version 3.0.1 of 'nlmixr2'
#> remotes::install_version("nlmixr2", version = "3.0.1") 
#> 
#> # Package 'nonemem2rx' is not installed, but known to enhance nlmixr2/babelmixr2
#> 
#> # Package 'monolix2rx' is not installed, but known to enhance nlmixr2/babelmixr2
#> 
#> # Package 'babelmixr2' is not installed, but known to enhance nlmixr2/babelmixr2
#> 
#> # Package 'PopED' is not installed, but known to enhance nlmixr2/babelmixr2
#> 
#> # Install Package version 0.11.0 of 'PKNCA'
#> remotes::install_version("PKNCA", version = "0.11.0") 
#> 
#> # If all else fails you can try to install the version of nlmixr2 used to create the fit
#> Warning: decompression of an rxUi object from rxode2 < 4.0 requires qs which is not on CRAN
#>  
#>  
#> ℹ using original fit data for simulation
#> ── Solved rxode2 object ──
#> ── Parameters (value$params): ──
#> # A tibble: 12 × 8
#>    id      tka   tcl    tv add.sd  eta.ka  eta.cl   eta.v
#>    <fct> <dbl> <dbl> <dbl>  <dbl>   <dbl>   <dbl>   <dbl>
#>  1 1     0.465  1.01  3.46  0.695 -0.534   0.0622 -0.110 
#>  2 2     0.465  1.01  3.46  0.695  0.0223 -0.0423  0.172 
#>  3 3     0.465  1.01  3.46  0.695  0.128   0.183  -0.0706
#>  4 4     0.465  1.01  3.46  0.695  0.364  -0.475   0.0514
#>  5 5     0.465  1.01  3.46  0.695 -0.246  -0.0216 -0.0411
#>  6 6     0.465  1.01  3.46  0.695  0.732  -0.417   0.0144
#>  7 7     0.465  1.01  3.46  0.695  0.881  -0.148   0.248 
#>  8 8     0.465  1.01  3.46  0.695 -0.182   0.134  -0.0749
#>  9 9     0.465  1.01  3.46  0.695 -0.0343  0.756   0.0456
#> 10 10    0.465  1.01  3.46  0.695  0.254  -0.178  -0.151 
#> 11 11    0.465  1.01  3.46  0.695  0.740   0.239   0.204 
#> 12 12    0.465  1.01  3.46  0.695  0.0644  0.0680 -0.0701
#> ── Initial Conditions (value$inits): ──
#>  depot center 
#>      0      0 
#> ── First part of data (object): ──
#> # A tibble: 132 × 10
#>      id  time    ka    cl     v    cp ipredSim   sim  depot center
#>   <int> <dbl> <dbl> <dbl> <dbl> <dbl>    <dbl> <dbl>  <dbl>  <dbl>
#> 1     1  0    0.933  2.91  28.4  0        0    0.865 320.      0  
#> 2     1  0.25 0.933  2.91  28.4  2.31     2.31 0.970 253.     65.7
#> 3     1  0.57 0.933  2.91  28.4  4.50     4.50 5.11  188.    128. 
#> 4     1  1.12 0.933  2.91  28.4  6.84     6.84 6.19  113.    194. 
#> 5     1  2.02 0.933  2.91  28.4  8.37     8.37 9.18   48.6   238. 
#> 6     1  3.82 0.933  2.91  28.4  8.20     8.20 8.52    9.06  233. 
#> # ℹ 126 more rows

# }
```
