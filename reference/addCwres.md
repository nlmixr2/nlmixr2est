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
#> → compress parHistData in nlmixr2 object, save 8584
#> → compress phiM in nlmixr2 object, save 446912

print(f)
#> ── nlmixr² SAEM OBJF by FOCEi approximation ──
#> 
#>  Gaussian/Laplacian Likelihoods: AIC() or $objf etc. 
#>  FOCEi CWRES & Likelihoods: addCwres() 
#> 
#> ── Time (sec $time): ──
#> 
#>             setup   optimize covariance preprocess configure  saem postprocess
#> elapsed 0.7806567 2.5959e-05 0.01200737      0.063      0.88 3.211       0.825
#>         table compress   other
#> elapsed 0.064    0.128 0.43431
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.     SE %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tka              0.454  0.192 42.2       1.57 (1.08, 2.29)     69.6     -0.514 
#> tcl               1.01 0.0850 8.40       2.75 (2.33, 3.25)     27.5       3.02 
#> tv         log V  3.45 0.0451 1.31       31.6 (28.9, 34.5)     13.2       11.2 
#> add.sd           0.700      0    0    0.700 (0.700, 0.700)                     
#>  
#>   Covariance Type ($covMethod): linFim
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 132 × 18
#>   ID     TIME    DV  PRED    RES IPRED    IRES   IWRES eta.ka eta.cl   eta.v
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>
#> 1 1      0     0.74  0     0.74   0     0.74    1.06    0.103 -0.486 -0.0781
#> 2 1      0.25  2.84  3.21 -0.371  3.81 -0.968  -1.38    0.103 -0.486 -0.0781
#> 3 1      0.57  6.57  5.62  0.946  6.61 -0.0394 -0.0563  0.103 -0.486 -0.0781
#> # ℹ 129 more rows
#> # ℹ 7 more variables: depot <dbl>, central <dbl>, ka <dbl>, cl <dbl>, v <dbl>,
#> #   tad <dbl>, dosenum <dbl>

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
#> FOCEi 116.9138 373.5136 393.6932      -179.7568        18.80749         1.42463
#> 
#> ── Time (sec $time): ──
#> 
#>             setup   optimize covariance preprocess configure  saem postprocess
#> elapsed 0.7806567 2.5959e-05 0.01200737      0.063      0.88 3.211       0.825
#>         table compress   other
#> elapsed 0.064    0.128 0.43431
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.     SE %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tka              0.454  0.192 42.2       1.57 (1.08, 2.29)     69.6     -0.514 
#> tcl               1.01 0.0850 8.40       2.75 (2.33, 3.25)     27.5       3.02 
#> tv         log V  3.45 0.0451 1.31       31.6 (28.9, 34.5)     13.2       11.2 
#> add.sd           0.700      0    0    0.700 (0.700, 0.700)                     
#>  
#>   Covariance Type ($covMethod): linFim
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 132 × 22
#>   ID     TIME    DV  PRED    RES IPRED    IRES   IWRES eta.ka eta.cl   eta.v
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>
#> 1 1      0     0.74  0     0.74   0     0.74    1.06    0.103 -0.486 -0.0781
#> 2 1      0.25  2.84  3.21 -0.371  3.81 -0.968  -1.38    0.103 -0.486 -0.0781
#> 3 1      0.57  6.57  5.62  0.946  6.61 -0.0394 -0.0563  0.103 -0.486 -0.0781
#> # ℹ 129 more rows
#> # ℹ 11 more variables: depot <dbl>, central <dbl>, ka <dbl>, cl <dbl>, v <dbl>,
#> #   tad <dbl>, dosenum <dbl>, WRES <dbl>, CPRED <dbl>, CRES <dbl>, CWRES <dbl>

# Note this also adds the FOCEi objective function
# }
```
