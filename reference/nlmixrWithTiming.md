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
#> → compress parHistData in nlmixr2 object, save 8984
#> → compress phiM in nlmixr2 object, save 448192

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
#>              setup   optimize covariance preprocess configure  saem postprocess
#> elapsed 0.09928175 2.4175e-05 0.01100475      0.068     0.215 1.652       0.399
#>         table compress     other time2 time1
#> elapsed  0.06    0.129 0.1256893 1.002 1.001
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.     SE     %RSE Back-transformed(95%CI) BSV(CV%)
#> tka              0.444  0.194     43.6       1.56 (1.07, 2.28)    70.48
#> tcl               1.12 0.0259     2.31       3.07 (2.92, 3.23)    29.54
#> tv         log V  3.45 0.0452     1.31       31.4 (28.7, 34.3)    13.18
#> add.sd           0.698   69.5 9.96e+03       0.698 (-136, 137)         
#>        Shrink(SD)%
#> tka       -0.467% 
#> tcl         8.37% 
#> tv          12.2% 
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
#> 1 1      0     0.74  0     0.74   0     0.74   1.06   0.109 -0.580 -0.0762  320.
#> 2 1      0.25  2.84  3.23 -0.388  3.83 -0.994 -1.42   0.109 -0.580 -0.0762  207.
#> 3 1      0.57  6.57  5.72  0.846  6.72 -0.149 -0.214  0.109 -0.580 -0.0762  119.
#> # ℹ 129 more rows
#> # ℹ 6 more variables: central <dbl>, ka <dbl>, cl <dbl>, v <dbl>, tad <dbl>,
#> #   dosenum <dbl>

# }
```
