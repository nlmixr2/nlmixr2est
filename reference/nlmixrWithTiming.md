# Time a part of a nlmixr operation and add to nlmixr object

Time a part of a nlmixr operation and add to nlmixr object

## Usage

``` r
nlmixrWithTiming(name, code, envir = NULL)
```

## Arguments

- name:

  Name of the timing to be integrated

- code:

  Code to be evaluated and timed

- envir:

  can be either the nlmixr2 fit data, the nlmixr2 fit environment or
  NULL, which implies it is going to be added to the nlmixr fit when it
  is finalized. If the function is being called after a fit is created,
  please supply this environmental variable

## Value

Result of code

## Author

Matthew L. Fidler

## Examples

``` r
# \donttest{

one.cmt <- function() {
 ini({
   ## You may label each parameter with a comment
   tka <- 0.45 # Ka
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
fit <- nlmixr(one.cmt, theo_sd, est="saem")
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
#> Error in .inParsLookup[[n]] : subscript out of bounds
#> Error: subscript out of bounds

nlmixrWithTiming("time1", {
   Sys.sleep(1)
   # note this can be nested, time1 will exclude the timing from time2
   nlmixrWithTiming("time2", {
      Sys.sleep(1)
   }, envir=fit)
}, envir=fit)
#> Error: object 'fit' not found

print(fit)
#> Error: object 'fit' not found

# }
```
