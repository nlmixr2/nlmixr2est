# Add CWRES

This returns a new fit object with CWRES attached

## Usage

``` r
addCwres(fit, focei = TRUE, updateObject = TRUE, envir = parent.frame(1))
```

## Arguments

- fit:

  nlmixr2 fit without WRES/CWRES

- focei:

  Boolean indicating if the focei objective function is added. If not
  the foce objective function is added.

- updateObject:

  Boolean indicating if the original fit object should be updated. By
  default this is true.

- envir:

  Environment that should be checked for object to update. By default
  this is the global environment.

## Value

fit with CWRES

## Author

Matthew L. Fidler

## Examples

``` r

# \donttest{

one.cmt <- function() {
  ini({
    ## You may label each parameter with a comment
    tka <- 0.45 # Log Ka
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

f <- try(nlmixr2(one.cmt, theo_sd, "saem"))
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of saem model...
#> ✔ done
#> → finding duplicate expressions in saem model...
#> ✔ done
#>  
#>  
#> ℹ calculate uninformed etas
#> ℹ done
#> Calculating covariance matrix
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
#> → compress parHistData in nlmixr2 object, save 9560
#> → compress phiM in nlmixr2 object, save 64384

print(f)
#> ── nlmixr² SAEM OBJF by FOCEi approximation ──
#> 
#>  Gaussian/Laplacian Likelihoods: AIC() or $objf etc. 
#>  FOCEi CWRES & Likelihoods: addCwres() 
#> 
#> ── Time (sec $time): ──
#> 
#>             setup   optimize covariance preprocess configure  saem postprocess
#> elapsed 0.1165397 3.7189e-05 0.04401476      0.082     0.871 2.814       0.948
#>         table compress    other
#> elapsed 0.093    0.076 2.724408
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.     SE %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tka              0.462  0.188 40.8       1.59 (1.10, 2.30)     70.3       1.89 
#> tcl              0.983 0.0847 8.62       2.67 (2.26, 3.16)     27.3       4.88 
#> tv         log V  3.45 0.0430 1.24       31.7 (29.1, 34.4)     13.4       9.35 
#> add.sd           0.697 0.0468 6.72    0.697 (0.605, 0.789)                     
#>  
#>   Covariance Type ($covMethod): sa
#>   Fixed parameter correlations in $cor
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Information about run found ($runInfo):
#>    • IWRES shrinkage is strongly negative; check the residual model 
#>    • closed form: tka, tv; numerical (optimize): tcl 
#>    • IWRES positively autocorrelated within subjects: Durbin-Watson 0.25 
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 132 × 18
#>   ID     TIME    DV  PRED    RES IPRED    IRES   IWRES eta.ka eta.cl   eta.v
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>
#> 1 1      0     0.74  0     0.74   0     0.74    1.06   0.0911 -0.452 -0.0847
#> 2 1      0.25  2.84  3.22 -0.377  3.80 -0.961  -1.38   0.0911 -0.452 -0.0847
#> 3 1      0.57  6.57  5.61  0.963  6.58 -0.0122 -0.0175 0.0911 -0.452 -0.0847
#> # ℹ 129 more rows
#> # ℹ 7 more variables: depot <dbl>, central <dbl>, ka <dbl>, cl <dbl>, v <dbl>,
#> #   tad <dbl>, dosenum <dbl>

# even though you may have forgotten to add the cwres, you can add it to the data.frame:

if (!inherits(f, "try-error")) {
  f <- try(addCwres(f))
  print(f)
}
#> → Calculating residuals/tables
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → calculate ∂(f)/∂(η)
#> → calculate ∂(R²)/∂(η)
#> → finding duplicate expressions in inner model...
#> → finding duplicate expressions in EBE model...
#> → compiling inner model...
#>  
#>  
#> ✔ done
#> → finding duplicate expressions in FD model...
#> → compiling EBE model...
#>  
#>  
#> ✔ done
#> → compiling events FD model...
#>  
#>  
#> ✔ done
#> ✔ done
#>  
#>  
#> ── nlmixr² SAEM OBJF by FOCEi approximation ──
#> 
#>           OBJF     AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> FOCEi 116.9652 373.565 393.7446      -179.7825        432.9067        1.568121
#> 
#> ── Time (sec $time): ──
#> 
#>             setup   optimize covariance preprocess configure  saem postprocess
#> elapsed 0.1165397 3.7189e-05 0.04401476      0.082     0.871 2.814       0.948
#>         table compress    other foceiLik
#> elapsed 0.484    0.076 2.724408    0.279
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.     SE %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tka              0.462  0.188 40.8       1.59 (1.10, 2.30)     70.3       1.89 
#> tcl              0.983 0.0847 8.62       2.67 (2.26, 3.16)     27.3       4.88 
#> tv         log V  3.45 0.0430 1.24       31.7 (29.1, 34.4)     13.4       9.35 
#> add.sd           0.697 0.0468 6.72    0.697 (0.605, 0.789)                     
#>  
#>   Covariance Type ($covMethod): sa
#>   Fixed parameter correlations in $cor
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Information about run found ($runInfo):
#>    • IWRES shrinkage is strongly negative; check the residual model 
#>    • closed form: tka, tv; numerical (optimize): tcl 
#>    • IWRES positively autocorrelated within subjects: Durbin-Watson 0.25 
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 132 × 22
#>   ID     TIME    DV  PRED    RES IPRED    IRES   IWRES eta.ka eta.cl   eta.v
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>
#> 1 1      0     0.74  0     0.74   0     0.74    1.06   0.0911 -0.452 -0.0847
#> 2 1      0.25  2.84  3.22 -0.377  3.80 -0.961  -1.38   0.0911 -0.452 -0.0847
#> 3 1      0.57  6.57  5.61  0.963  6.58 -0.0122 -0.0175 0.0911 -0.452 -0.0847
#> # ℹ 129 more rows
#> # ℹ 11 more variables: depot <dbl>, central <dbl>, ka <dbl>, cl <dbl>, v <dbl>,
#> #   tad <dbl>, dosenum <dbl>, WRES <dbl>, CPRED <dbl>, CRES <dbl>, CWRES <dbl>

# Note this also adds the FOCEi objective function
# }
```
