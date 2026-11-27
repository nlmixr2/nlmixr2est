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
#> → compress origData in nlmixr2 object, save 6584
#> → compress parHistData in nlmixr2 object, save 9560
#> → compress phiM in nlmixr2 object, save 64160

print(f)
#> ── nlmixr² SAEM OBJF by FOCEi approximation ──
#> 
#>  Gaussian/Laplacian Likelihoods: AIC() or $objf etc. 
#>  FOCEi CWRES & Likelihoods: addCwres() 
#> 
#> ── Time (sec $time): ──
#> 
#>              setup   optimize covariance preprocess configure  saem postprocess
#> elapsed 0.08161661 2.4807e-05  0.0320132       0.07     0.149 2.179       0.211
#>         compress    other
#> elapsed    0.057 2.124345
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
#>    • closed form: tka, tv; numerical (optimize): tcl 
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
#>              setup   optimize covariance preprocess configure  saem postprocess
#> elapsed 0.08161661 2.4807e-05  0.0320132       0.07     0.149 2.179       0.211
#>         compress    other
#> elapsed    0.057 2.124345
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
#>    • closed form: tka, tv; numerical (optimize): tcl 
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

# }
```
