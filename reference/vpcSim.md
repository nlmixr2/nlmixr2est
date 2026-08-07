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
#> 1      1  1 0.00  0.000000  0.9517617 0.00 3.199920e+02   0.0000             2
#> 2      1  1 0.25 10.021812 10.1804443 0.25 4.529562e+01 269.9694             3
#> 3      1  1 0.57 11.194311 11.3923384 0.57 3.708772e+00 301.5544             4
#> 4      1  1 1.12 10.687702 10.2956291 1.12 5.026292e-02 287.9072             5
#> 5      1  1 2.02  9.713337 10.6677170 2.02 4.411125e-05 261.6596             6
#> 6      1  1 3.82  8.020164  8.4816159 3.82 3.397446e-11 216.0486             7
#>   rxLambda rxYj rxLow rxHi     pred
#> 1        1    2     0    1 0.000000
#> 2        1    2     0    1 3.263024
#> 3        1    2     0    1 5.830051
#> 4        1    2     0    1 7.864517
#> 5        1    2     0    1 8.504253
#> 6        1    2     0    1 7.619392

# }
```
