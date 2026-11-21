# Attach likelihood penalties to a model

Used by model extensions (for example from an rxode2 user-function UI
handler) so the penalty is active only while this model is fit.

## Usage

``` r
nlmixr2LikPenaltyAdd(ui, ..., replace = FALSE)
```

## Arguments

- ui:

  rxode2 model (or function)

- ...:

  \`nlmixr2LikPenaltySpec()\` objects

- replace:

  replace the penalties already attached instead of adding to them

## Value

the rxode2 ui with the penalties stored in \`ui\$nlmixr2LikPenalty\`

## Author

Matthew L. Fidler

## Examples

``` r
one.cmt <- function() {
  ini({
    tka <- 0.45; tcl <- 1; tv <- 3.45
    eta.ka ~ 0.6
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl)
    v <- exp(tv)
    linCmt() ~ add(add.sd)
  })
}
ui <- nlmixr2LikPenaltyAdd(one.cmt,
  nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = 10))
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
ui$nlmixr2LikPenalty
#> [[1]]
#> likelihood penalty thetaL2 on tka par: 10 inObjf: NA 
#> 
```
