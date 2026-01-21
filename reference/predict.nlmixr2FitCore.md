# Predict method for nlmixr2 fit core objects

This function generates predictions from an \`nlmixr2FitCore\` object.
It allows for both population-level and individual-level predictions
based on the specified \`level\` parameter.

## Usage

``` r
# S3 method for class 'nlmixr2FitCore'
predict(object, ..., level = c("population", "individual"))
```

## Arguments

- object:

  nlmixr2 fit core object to predict

- ...:

  additional arguments passed to rxode2::rxSolve or nlmixr2; matching
  other \`predict\` methods, these can include \`newdata\` and
  \`rxControl\` settings

- level:

  the prediction level; one of \`"population"\` (default) or
  \`"individual"\`; numeric values \`0\` and \`1\` are also accepted

## Value

A data frame with predictions

## Examples

``` r
# \donttest{

one.compartment <- function() {
 ini({
  tka <- log(1)
  tcl <- log(10)
  tv <- log(35)
  eta.ka ~ 0.1
  eta.cl ~ 0.1
  eta.v ~ 0.1
  add.sd <- 0.1
 })
 model({
  ka <- exp(tka + eta.ka)
  cl <- exp(tcl + eta.cl)
  v <- exp(tv + eta.v)
  d/dt(depot) = -ka * depot
  d/dt(center) = ka * depot - cl / v * center
  cp = center / v
  cp ~ add(add.sd)
 })
}

# The fit is performed by the function nlmixr/nlmix2 specifying
# the model, data and estimate
fit <- nlmixr2(one.compartment, theo_sd, est = "focei",
               foceiControl(maxOuterIterations = 0L))
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → calculate jacobian
#> → calculate sensitivities
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
#> → optimizing duplicate expressions in FD model...
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

# Population predictions
ppred <- predict(fit, theo_sd, level="population")
#> ℹ population predictions requested (`level="population"`)
#> ℹ using new data for predictions
#>  
#>  
#>  
#>  
#>  
#>  

# Individual predictions
ipred <- predict(fit, theo_sd, level="individual")
#> ℹ individual predictions requested (`level="individual"`)
#> ℹ using new data for predictions
#>  
#>  

# }
```
