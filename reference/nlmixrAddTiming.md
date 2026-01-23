# Manually add time to a nlmixr2 object

Manually add time to a nlmixr2 object

## Usage

``` r
nlmixrAddTiming(object, name, time)
```

## Arguments

- object:

  nlmixr2 object

- name:

  string of the timing name

- time:

  time (in seconds)

## Value

Nothing, called for side effects

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

# will add to the current setup
nlmixrAddTiming(fit, "setup", 3)
#> Error: object 'fit' not found

# Add a new item to the timing dataframe
nlmixrAddTiming(fit, "new", 3)
#> Error: object 'fit' not found

# }
```
