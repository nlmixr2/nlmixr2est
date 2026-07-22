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
#> → optimizing duplicate expressions in saem predOnly model 2...
#> ✔ done
#>  
#>  
#> → compress origData in nlmixr2 object, save 6584
#> → compress parHistData in nlmixr2 object, save 8592
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
#> elapsed 0.1006042 2.8452e-05 0.01100467      0.072     0.227 7.082       0.404
#>         compress     other
#> elapsed    0.134 0.1343627
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.        SE      %RSE Back-transformed(95%CI) BSV(CV%)
#> tka              0.454     0.192      42.3       1.57 (1.08, 2.29)     69.6
#> tcl               1.01    0.0851      8.40       2.75 (2.33, 3.25)     27.5
#> tv         log V  3.45    0.0451      1.31       31.6 (28.9, 34.5)     13.2
#> add.sd           0.700 4.66e-310 6.66e-308    0.700 (0.700, 0.700)         
#>        Shrink(SD)%
#> tka         -1.53 
#> tcl          3.98 
#> tv           9.72 
#> add.sd            
#>  
#>   Covariance Type ($covMethod): linFim
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
#> elapsed 0.1006042 2.8452e-05 0.01100467      0.072     0.227 7.082       0.404
#>         compress     other
#> elapsed    0.134 0.1343627
#> 
#> ── Population Parameters ($parFixed or $parFixedDf): ──
#> 
#>        Parameter  Est.        SE      %RSE Back-transformed(95%CI) BSV(CV%)
#> tka              0.454     0.192      42.3       1.57 (1.08, 2.29)     69.6
#> tcl               1.01    0.0851      8.40       2.75 (2.33, 3.25)     27.5
#> tv         log V  3.45    0.0451      1.31       31.6 (28.9, 34.5)     13.2
#> add.sd           0.700 4.66e-310 6.66e-308    0.700 (0.700, 0.700)         
#>        Shrink(SD)%
#> tka         -1.53 
#> tcl          3.98 
#> tv           9.72 
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
#>   ID     TIME    DV  PRED    RES IPRED    IRES   IWRES eta.ka eta.cl   eta.v
#>   <fct> <dbl> <dbl> <dbl>  <dbl> <dbl>   <dbl>   <dbl>  <dbl>  <dbl>   <dbl>
#> 1 1      0     0.74  0     0.74   0     0.74    1.06    0.107 -0.484 -0.0782
#> 2 1      0.25  2.84  3.21 -0.371  3.82 -0.977  -1.40    0.107 -0.484 -0.0782
#> 3 1      0.57  6.57  5.62  0.946  6.62 -0.0514 -0.0734  0.107 -0.484 -0.0782
#> # ℹ 129 more rows
#> # ℹ 7 more variables: depot <dbl>, central <dbl>, ka <dbl>, cl <dbl>, v <dbl>,
#> #   tad <dbl>, dosenum <dbl>

# }
```
