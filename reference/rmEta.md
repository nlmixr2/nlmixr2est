# Remove an eta from the model

Remove an eta from the model

## Usage

``` r
rmEta(ui, eta)
```

## Arguments

- ui:

  rxode2 user interface

- eta:

  eta to remove

## Value

ui model with eta removed

## Author

Matthew L. Fidler

## Examples

``` r
mod <- function ()  {
 description <- "One compartment PK model with linear clearance"
 ini({
   lka <- 0.45
   lcl <- 1
   lvc <- 3.45
    propSd <- c(0, 0.5)
    etaKa ~ 0.1
  })
 model({
   ka <- exp(lka + etaKa)
   cl <- exp(lcl)
   vc <- exp(lvc)
   Cc <- linCmt()
   Cc ~ prop(propSd)
 })
}

mod |> rmEta("etaKa")
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> Error in zeroEtasLookup[[.n]]: subscript out of bounds

# This can also remove more than one eta

mod <- function ()  {
 description <- "One compartment PK model with linear clearance"
 ini({
   lka <- 0.45
   lcl <- 1
   lvc <- 3.45
   propSd <- c(0, 0.5)
   etaKa ~ 0.1
   etaCl ~ 0.2
   etaVc ~ 0.3
  })
 model({
   ka <- exp(lka + etaKa)
   cl <- exp(lcl + etaCl)
   vc <- exp(lvc + etaVc)
   Cc <- linCmt()
   Cc ~ prop(propSd)
 })
}

mod |> rmEta(c("etaKa", "etaCl"))
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> Error in zeroEtasLookup[[.n]]: subscript out of bounds
```
