# VPC simulation

VPC simulation

## Usage

``` r
vpcSim(
  object,
  ...,
  keep = NULL,
  n = 300,
  pred = FALSE,
  seed = 1009,
  nretry = 50,
  minN = 10,
  normRelated = TRUE
)
```

## Arguments

- object:

  This is the nlmixr2 fit object

- ...:

  Other arguments sent to \`rxSolve()\`

- keep:

  Column names to keep in the output simulated dataset

- n:

  Number of simulations

- pred:

  Should predictions be added to the simulation

- seed:

  Seed to set for the VPC simulation

- nretry:

  Number of times to retry the simulation if there is NA values in the
  simulation

- minN:

  With retries, the minimum number of studies to restimulate (by default
  10)

- normRelated:

  should the VPC style simulation be for normal related variables only

## Value

data frame of the VPC simulation

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

fit <- nlmixr(one.cmt, theo_sd, est="focei")
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → Calculating residuals/tables
#> ✔ done

head(vpcSim(fit, pred=TRUE))
#>  
#>  
#>   sim.id id time     ipred        sim  tad        depot  central nlmixrRowNums
#> 1      1  1 0.00  0.000000  0.9485653 0.00 3.199920e+02   0.0000             2
#> 2      1  1 0.25 10.041901 10.2000012 0.25 4.583738e+01 269.4264             3
#> 3      1  1 0.57 11.234563 11.4319253 0.57 3.810684e+00 301.4257             4
#> 4      1  1 1.12 10.727509 10.3367536 1.12 5.301276e-02 287.8214             5
#> 5      1  1 2.02  9.746305 10.6974795 2.02 4.855911e-05 261.4954             6
#> 6      1  1 3.82  8.041929  8.5018306 3.82 4.074284e-11 215.7667             7
#>   rxLambda rxYj rxLow rxHi     pred
#> 1        1    2     0    1 0.000000
#> 2        1    2     0    1 3.257260
#> 3        1    2     0    1 5.823935
#> 4        1    2     0    1 7.863516
#> 5        1    2     0    1 8.510299
#> 6        1    2     0    1 7.628016

# }
```
