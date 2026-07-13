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
#>  
#>  
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of saem model...
#> ✔ done
#> → finding duplicate expressions in saem model...
#> ✔ done
#> ℹ calculate uninformed etas
#> ℹ done
#> rxode2 5.1.3 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
#> 
#> Attaching package: ‘rxode2’
#> The following objects are masked from ‘package:nlmixr2est’:
#> 
#>     boxCox, yeoJohnson
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

print(f)
#> ── nlmixr² SAEM OBJF by FOCEi approximation ──
#> 
#>  Gaussian/Laplacian Likelihoods: AIC() or $objf etc. 
#>  FOCEi CWRES & Likelihoods: addCwres() 
#> 
#> ── Time (sec $time): ──
#> 
#>            setup   optimize covariance preprocess configure  saem postprocess
#> elapsed 0.693426 3.1439e-05 0.01301097      0.079     1.064 1.624       1.044
#>         table compress     other
#> elapsed 0.337    0.134 0.4045316
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.        SE      %RSE Back-transformed(95%CI) BSV(CV%)
#> tka              0.452     0.192      42.4       1.57 (1.08, 2.29)     69.7
#> tcl               1.04    0.0242      2.33        2.83 (2.7, 2.97)     28.0
#> tv         log V  3.45    0.0445      1.29       31.5 (28.9, 34.4)     13.0
#> add.sd           0.699 6.91e-310 9.88e-308    0.699 (0.699, 0.699)         
#>        Shrink(SD)%
#> tka       -0.830% 
#> tcl         3.23% 
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
#> 1 1      0     0.74  0     0.74   0     0.74   1.06   0.118 -0.513 -0.0739  320.
#> 2 1      0.25  2.84  3.22 -0.378  3.84 -1.00  -1.44   0.118 -0.513 -0.0739  206.
#> 3 1      0.57  6.57  5.65  0.916  6.67 -0.101 -0.144  0.118 -0.513 -0.0739  117.
#> # ℹ 129 more rows
#> # ℹ 6 more variables: central <dbl>, ka <dbl>, cl <dbl>, v <dbl>, tad <dbl>,
#> #   dosenum <dbl>

# even though you may have forgotten to add the cwres, you can add it to the data.frame:

if (!inherits(f, "try-error")) {
  f <- try(addCwres(f))
  print(f)
}
#>  
#>  
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → calculate ∂(f)/∂(η)
#> → calculate ∂(R²)/∂(η)
#> → finding duplicate expressions in inner model...
#> → optimizing duplicate expressions in inner model...
#> → finding duplicate expressions in EBE model...
#> → optimizing duplicate expressions in EBE model...
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
#> → Calculating residuals/tables
#> ✔ done
#> ── nlmixr² SAEM OBJF by FOCEi approximation ──
#> 
#>           OBJF      AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> FOCEi 117.0711 373.6709 393.8505      -179.8354        64.21327        1.427592
#> 
#> ── Time (sec $time): ──
#> 
#>            setup   optimize covariance preprocess configure  saem postprocess
#> elapsed 0.693426 3.1439e-05 0.01301097      0.079     1.064 1.624       1.044
#>         table compress     other
#> elapsed 0.337    0.134 0.4045316
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.        SE      %RSE Back-transformed(95%CI) BSV(CV%)
#> tka              0.452     0.192      42.4       1.57 (1.08, 2.29)     69.7
#> tcl               1.04    0.0242      2.33        2.83 (2.7, 2.97)     28.0
#> tv         log V  3.45    0.0445      1.29       31.5 (28.9, 34.4)     13.0
#> add.sd           0.699 6.91e-310 9.88e-308    0.699 (0.699, 0.699)         
#>        Shrink(SD)%
#> tka       -0.830% 
#> tcl         3.23% 
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
#> # A tibble: 132 × 22
#>   ID     TIME    DV  PRED    RES IPRED   IRES  IWRES eta.ka eta.cl   eta.v depot
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>  <dbl>  <dbl>   <dbl> <dbl>
#> 1 1      0     0.74  0     0.74   0     0.74   1.06   0.118 -0.513 -0.0739  320.
#> 2 1      0.25  2.84  3.22 -0.378  3.84 -1.00  -1.44   0.118 -0.513 -0.0739  206.
#> 3 1      0.57  6.57  5.65  0.916  6.67 -0.101 -0.144  0.118 -0.513 -0.0739  117.
#> # ℹ 129 more rows
#> # ℹ 10 more variables: central <dbl>, ka <dbl>, cl <dbl>, v <dbl>, tad <dbl>,
#> #   dosenum <dbl>, WRES <dbl>, CPRED <dbl>, CRES <dbl>, CWRES <dbl>

# Note this also adds the FOCEi objective function
# }
```
