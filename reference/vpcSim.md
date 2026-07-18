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
#> covType="analytic": a linCmt() model is out of analytic-covariance scope; using the finite-difference covariance instead
#> → Calculating residuals/tables
#> ✔ done

head(vpcSim(fit, pred=TRUE))
#>  
#>  
#>   sim.id id time     ipred        sim  tad        depot  central nlmixrRowNums
#> 1      1  1 0.00  0.000000  0.9527709 0.00 3.199920e+02   0.0000             2
#> 2      1  1 0.25 10.024480 10.1832802 0.25 4.511370e+01 270.1492             3
#> 3      1  1 0.57 11.191146 11.3893827 0.57 3.674897e+00 301.5897             4
#> 4      1  1 1.12 10.683905 10.2914170 1.12 4.936484e-02 287.9201             5
#> 5      1  1 2.02  9.710472 10.6658638 2.02 4.269994e-05 261.6871             6
#> 6      1  1 3.82  8.018815  8.4807565 3.82 3.194820e-11 216.0987             7
#>   rxLambda rxYj rxLow rxHi     pred
#> 1        1    2     0    1 0.000000
#> 2        1    2     0    1 3.266456
#> 3        1    2     0    1 5.834730
#> 4        1    2     0    1 7.868329
#> 5        1    2     0    1 8.505958
#> 6        1    2     0    1 7.619941

# }
```
