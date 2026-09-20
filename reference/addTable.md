# Add table information to nlmixr2 fit object without tables

Add table information to nlmixr2 fit object without tables

## Usage

``` r
addTable(
  object,
  updateObject = FALSE,
  data = object$dataSav,
  thetaEtaParameters = object$foceiThetaEtaParameters,
  table = tableControl(),
  keep = NULL,
  drop = NULL,
  envir = parent.frame(1)
)
```

## Arguments

- object:

  nlmixr2 family of objects

- updateObject:

  Update the object (default FALSE)

- data:

  Saved data from

- thetaEtaParameters:

  Internal theta/eta parameters

- table:

  a \`tableControl()\` list of options

- keep:

  Character Vector of items to keep

- drop:

  Character Vector of items to drop or NULL

- envir:

  Environment to search for updating

## Value

Fit with table information attached

## Author

Matthew Fidler

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

# run without tables step
f <- nlmixr2(one.cmt, theo_sd, "saem", control=list(calcTables=FALSE))
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
#> covMethod="sa" could not be computed; using the linearized FIM
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
#> → compress origData in nlmixr2 object, save 6584
#> → compress parHistData in nlmixr2 object, save 8616
#> → compress phiM in nlmixr2 object, save 439912

print(f)
#> ── nlmixr² SAEM OBJF by FOCEi approximation ──
#> 
#>  Gaussian/Laplacian Likelihoods: AIC() or $objf etc. 
#>  FOCEi CWRES & Likelihoods: addCwres() 
#> 
#> ── Time (sec $time): ──
#> 
#>             setup   optimize covariance preprocess configure  saem postprocess
#> elapsed 0.1096606 2.8182e-05 0.01400444      0.085     0.147 4.266       0.338
#>         compress     other
#> elapsed    0.134 0.1473068
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.     SE %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tka              0.459  0.193 42.0       1.58 (1.08, 2.31)     70.0     -0.106 
#> tcl               1.01 0.0843 8.35       2.74 (2.33, 3.24)     27.2       4.20 
#> tv         log V  3.45 0.0453 1.31       31.6 (28.9, 34.6)     13.3       9.01 
#> add.sd           0.696 0.0500 7.18    0.696 (0.598, 0.794)                     
#>  
#>   Covariance Type ($covMethod): linFim
#>   Fixed parameter correlations in $cor
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Censoring ($censInformation): No censoring

# Now add the tables

f <- addTable(f)
#> → Calculating residuals/tables
#> ✔ done

print(f)
#> ── nlmixr² SAEM OBJF by FOCEi approximation ──
#> 
#>  Gaussian/Laplacian Likelihoods: AIC() or $objf etc. 
#>  FOCEi CWRES & Likelihoods: addCwres() 
#> 
#> ── Time (sec $time): ──
#> 
#>             setup   optimize covariance preprocess configure  saem postprocess
#> elapsed 0.1096606 2.8182e-05 0.01400444      0.085     0.147 4.266       0.338
#>         compress     other
#> elapsed    0.134 0.1473068
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.     SE %RSE Back-transformed(95%CI) BSV(CV%) Shrink(SD)%
#> tka              0.459  0.193 42.0       1.58 (1.08, 2.31)     70.0     -0.106 
#> tcl               1.01 0.0843 8.35       2.74 (2.33, 3.24)     27.2       4.20 
#> tv         log V  3.45 0.0453 1.31       31.6 (28.9, 34.6)     13.3       9.01 
#> add.sd           0.696 0.0500 7.18    0.696 (0.598, 0.794)                     
#>  
#>   Covariance Type ($covMethod): linFim
#>   Fixed parameter correlations in $cor
#>   No correlations in between subject variability (BSV) matrix
#>   Full BSV covariance ($omega) or correlation ($omegaR; diagonals=SDs) 
#>   Distribution stats (mean/skewness/kurtosis/p-value) available in $shrink 
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 132 × 18
#>   ID     TIME    DV  PRED    RES IPRED    IRES   IWRES eta.ka eta.cl   eta.v
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>
#> 1 1      0     0.74  0     0.74   0     0.74    1.06   0.0971 -0.479 -0.0834
#> 2 1      0.25  2.84  3.22 -0.377  3.81 -0.975  -1.40   0.0971 -0.479 -0.0834
#> 3 1      0.57  6.57  5.63  0.942  6.62 -0.0507 -0.0728 0.0971 -0.479 -0.0834
#> # ℹ 129 more rows
#> # ℹ 7 more variables: depot <dbl>, central <dbl>, ka <dbl>, cl <dbl>, v <dbl>,
#> #   tad <dbl>, dosenum <dbl>

# }
```
