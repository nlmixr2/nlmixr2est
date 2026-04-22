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

  can be either the nlmixr2 fit data, the nlmixr2 fit environment or
  NULL, which implies it is going to be added to the nlmixr fit when it
  is finalized. If the function is being called after a fit is created,
  please supply this environmental variable

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
#> params:  tka rxBoundedTr.tcl tv  V(eta.ka)   V(eta.v)    V(eta.cl)   add.sd
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
#> → compress origData in nlmixr2 object, save 6592
#> → compress parHistData in nlmixr2 object, save 8280
#> → compress phiM in nlmixr2 object, save 429432
#> Warning:  mu-reference transform (exp) for `tcl` lost since bounded (and performance degraded)
#> Warning: to keep mu-referencing remove bounds or use control=list(boundedTransform=FALSE)

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
#>           setup covariance  saem table compress    other time2 time1
#> elapsed 0.00204   0.009016 2.733 0.074    0.058 0.705944 1.002 1.002
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.     SE %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tka              0.459  0.193 42.1       1.58 (1.08, 2.31)     70.3    -0.958% 
#> tcl               0.99 0.0851  8.6       2.69 (2.28, 3.18)     27.5      4.30% 
#> tv         log V  3.45 0.0452 1.31         31.6 (29, 34.6)     13.3      10.2% 
#> add.sd           0.699                               0.699                     
#>  
#>   Covariance Type ($covMethod): linFim
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Information about run found ($runInfo):
#>    • mu-reference transform (exp) for `tcl` lost since bounded (and performance degraded) 
#>    • to keep mu-referencing remove bounds or use control=list(boundedTransform=FALSE) 
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 132 × 18
#>   ID     TIME    DV  PRED    RES IPRED    IRES   IWRES eta.ka eta.cl   eta.v
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>
#> 1 1      0     0.74  0     0.74   0     0.74    1.06    0.102 -0.463 -0.0793
#> 2 1      0.25  2.84  3.21 -0.371  3.81 -0.967  -1.38    0.102 -0.463 -0.0793
#> 3 1      0.57  6.57  5.61  0.964  6.59 -0.0191 -0.0273  0.102 -0.463 -0.0793
#> # ℹ 129 more rows
#> # ℹ 7 more variables: depot <dbl>, central <dbl>, ka <dbl>, cl <dbl>, v <dbl>,
#> #   tad <dbl>, dosenum <dbl>

# }
```
