# NPDE calculation for nlmixr2

NPDE calculation for nlmixr2

## Usage

``` r
addNpde(
  object,
  updateObject = TRUE,
  table = tableControl(),
  ...,
  envir = parent.frame(1)
)
```

## Arguments

- object:

  nlmixr2 fit object

- updateObject:

  Boolean indicating if original object should be updated. By default
  this is TRUE.

- table:

  \`tableControl()\` list of options

- ...:

  Other ignored parameters.

- envir:

  Environment that should be checked for object to update. By default
  this is the global environment.

## Value

New nlmixr2 fit object

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

f <- nlmixr2(one.cmt, theo_sd, "saem")
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
#> → compress parHistData in nlmixr2 object, save 8296
#> → compress phiM in nlmixr2 object, save 312896

# even though you may have forgotten to add the NPDE, you can add it to the data.frame:

f <- addNpde(f)
#> → Add NPDE
#>  
#>  
#> ✔ done

# }
```
