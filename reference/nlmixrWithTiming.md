# Time a part of a nlmixr operation and add to nlmixr object

Time a part of a nlmixr operation and add to nlmixr object

## Usage

``` r
nlmixrWithTiming(name, code, envir = NULL)
```

## Arguments

- name:

  Name of the timing to be integrated

- code:

  Code to be evaluated and timed

- envir:

  nlmixr2 fit data, fit environment, or NULL (timing is added when the
  fit is finalized); supply this if called after a fit already exists

## Value

Result of code

## Author

Matthew L. Fidler

## Examples

``` r

# \donttest{

one.cmt <- function() {
 ini({
   ## You may label each parameter with a comment
   tka <- 0.45 # Ka
   tcl <- log(c(0, 2.7, 100)) # Log Cl
   ## This works with interactive models
   ## You may also label the preceding line with label("label text")
   tv <- 3.45; label("log V")
   ## the label("Label name") works with all models
   eta.ka ~ 0.6
   eta.cl ~ 0.3
   eta.v ~ 0.1
   add.sd <- 0.7
 })
 model({
   ka <- exp(tka + eta.ka)
   cl <- exp(tcl + eta.cl)
   v <- exp(tv + eta.v)
   linCmt() ~ add(add.sd)
 })
}
fit <- nlmixr(one.cmt, theo_sd, est="saem")
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#>  
#>  
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of saem model...
#> ✔ done
#> → finding duplicate expressions in saem model...
#> ✔ done
#> ℹ calculate uninformed etas
#> ℹ done
#> Calculating covariance matrix
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of saem model...
#> ✔ done
#> → finding duplicate expressions in saem predOnly model 0...
#> → finding duplicate expressions in saem predOnly model 1...
#> → finding duplicate expressions in saem predOnly model 2...
#> → optimizing duplicate expressions in saem predOnly model 2...
#> ✔ done
#>  
#>  
#> → Calculating residuals/tables
#> ✔ done
#> → compress origData in nlmixr2 object, save 6584
#> → compress parHistData in nlmixr2 object, save 8848
#> → compress phiM in nlmixr2 object, save 443520

nlmixrWithTiming("time1", {
   Sys.sleep(1)
   # note this can be nested, time1 will exclude the timing from time2
   nlmixrWithTiming("time2", {
      Sys.sleep(1)
   }, envir=fit)
}, envir=fit)

print(fit)
#> ── nlmixr² SAEM OBJF by FOCEi approximation ──
#> 
#>  Gaussian/Laplacian Likelihoods: AIC() or $objf etc. 
#>  FOCEi CWRES & Likelihoods: addCwres() 
#> 
#> ── Time (sec $time): ──
#> 
#>            setup   optimize covariance preprocess configure  saem postprocess
#> elapsed 0.105639 6.6005e-05 0.01000576      0.078     0.274 1.633       0.495
#>         table compress     other time2 time1
#> elapsed 0.062    0.135 0.1392893 1.002 1.001
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.        SE      %RSE Back-transformed(95%CI) BSV(CV%)
#> tka              0.452     0.192      42.4       1.57 (1.08, 2.29)     69.7
#> tcl               1.04    0.0864      8.30       2.83 (2.39, 3.36)     28.0
#> tv         log V  3.45    0.0445      1.29       31.5 (28.9, 34.4)     13.0
#> add.sd           0.699 2.12e-311 3.03e-309    0.699 (0.699, 0.699)         
#>        Shrink(SD)%
#> tka        -0.830 
#> tcl          3.23 
#> tv           12.2 
#> add.sd            
#>  
#>   Covariance Type ($covMethod): linFim
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 132 × 18
#>   ID     TIME    DV  PRED    RES IPRED   IRES  IWRES eta.ka eta.cl   eta.v depot
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl>
#> 1 1      0     0.74  0     0.74   0     0.74   1.06   0.118 -0.513 -0.0739  320.
#> 2 1      0.25  2.84  3.22 -0.378  3.84 -1.00  -1.44   0.118 -0.513 -0.0739  206.
#> 3 1      0.57  6.57  5.65  0.916  6.67 -0.101 -0.144  0.118 -0.513 -0.0739  117.
#> # ℹ 129 more rows
#> # ℹ 6 more variables: central <dbl>, ka <dbl>, cl <dbl>, v <dbl>, tad <dbl>,
#> #   dosenum <dbl>

# }
```
