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
#> params:  tka tcl tv  V(eta.ka)   V(eta.cl)   V(eta.v)    add.sd
#> Calculating covariance matrix
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of saem model...
#> ✔ done
#> → finding duplicate expressions in saem predOnly model 0...
#> → finding duplicate expressions in saem predOnly model 1...
#> → finding duplicate expressions in saem predOnly model 2...
#> ✔ done
#>  
#>  
#> → Calculating residuals/tables
#> ✔ done
#> → compress origData in nlmixr2 object, save 6584
#> → compress parHistData in nlmixr2 object, save 8280
#> → compress phiM in nlmixr2 object, save 312896

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
#>            setup covariance  saem table compress    other time2 time1
#> elapsed 0.001742   0.008012 1.922  0.08    0.047 0.559246 1.002 1.002
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.     SE %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tka               0.46  0.195 42.3       1.58 (1.08, 2.32)     71.0    -0.606% 
#> tcl               1.01 0.0838 8.31       2.74 (2.33, 3.23)     27.0      4.04% 
#> tv         log V  3.45 0.0453 1.31       31.6 (28.9, 34.6)     13.4      9.53% 
#> add.sd           0.695                               0.695                     
#>  
#>   Covariance Type ($covMethod): linFim
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Information about run found ($runInfo):
#>    • 'one.cmt' has the following user-defined boundaries: tcl which are ignored in 'saem' 
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 132 × 18
#>   ID     TIME    DV  PRED    RES IPRED   IRES  IWRES eta.ka eta.cl   eta.v depot
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl>
#> 1 1      0     0.74  0     0.74   0     0.74   1.06  0.0883 -0.477 -0.0846  320.
#> 2 1      0.25  2.84  3.27 -0.431  3.84 -0.996 -1.43  0.0883 -0.477 -0.0846  208.
#> 3 1      0.57  6.57  5.85  0.722  6.77 -0.201 -0.290 0.0883 -0.477 -0.0846  119.
#> # ℹ 129 more rows
#> # ℹ 6 more variables: central <dbl>, ka <dbl>, cl <dbl>, v <dbl>, tad <dbl>,
#> #   dosenum <dbl>

# }
```
