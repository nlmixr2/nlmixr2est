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

## Details

Besides \`EPRED\`, \`ERES\`, \`NPDE\`, \`NPD\`, \`PDE\` and \`PD\`, the
fit gets \`EWRES\`, the Monte Carlo weighted residual from \`EPRED\` (as
NONMEM's \`EWRES\`; Monolix calls it \`PWRES\`, and \`fit\$PWRES\` is an
alias): each subject's \`Sigma^(-1/2) (y - EPRED)\`, with \`EPRED\` and
\`Sigma\` the mean and covariance of the same simulations, decorrelated
as \`tableControl(ewresMethod=)\` asks. It is computed on the scale of
the residual error transformation, as NPDE. It keeps the residual's
magnitude, which the rank-based NPDE caps near \`qnorm(1/nsim)\`, but
uses only the first two moments: a skewed marginal (log-normal IIV,
proportional error) can show apparent bias even under the true model, so
NPDE remains the primary simulation diagnostic. A subject with a
censored observation, or with \`nsim\` no larger than its number of
observations, gets \`NA\`.

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
#> → compress phiM in nlmixr2 object, save 64160

# even though you may have forgotten to add the NPDE, you can add it to the data.frame:

f <- addNpde(f)
#> → Add NPDE
#>  
#>  
#> ✔ done

# }
```
