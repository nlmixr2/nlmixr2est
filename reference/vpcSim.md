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
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
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

head(vpcSim(fit, pred=TRUE))
#>  
#>  
#>   sim.id id time     ipred        sim  tad        depot  central nlmixrRowNums
#> 1      1  1 0.00  0.000000  0.9524295 0.00 3.199920e+02   0.0000             2
#> 2      1  1 0.25 10.010985 10.1697288 0.25 4.549704e+01 269.7767             3
#> 3      1  1 0.57 11.189486 11.3876518 0.57 3.746480e+00 301.5350             4
#> 4      1  1 1.12 10.684817 10.2924693 1.12 5.127199e-02 287.9352             5
#> 5      1  1 2.02  9.711361 10.6664104 2.02 4.572130e-05 261.7024             6
#> 6      1  1 3.82  8.019542  8.4813176 3.82 3.635761e-11 216.1112             7
#>   rxLambda rxYj rxLow rxHi     pred
#> 1        1    2     0    1 0.000000
#> 2        1    2     0    1 3.262285
#> 3        1    2     0    1 5.829159
#> 4        1    2     0    1 7.864093
#> 5        1    2     0    1 8.504682
#> 6        1    2     0    1 7.620541

# }
```
