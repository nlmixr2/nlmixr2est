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

mod %>% rmEta("etaKa")
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#>  ── rxode2-based solved PK 1-compartment model ────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka    lcl    lvc propSd 
#>   0.45   1.00   3.45   0.50 
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name Rate   Off Internal #
#> 1                  1            depot TRUE FALSE          1
#> 2                  2          central TRUE FALSE          2
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     description <- "One compartment PK model with linear clearance"
#>     ini({
#>         lka <- 0.45
#>         lcl <- 1
#>         lvc <- 3.45
#>         propSd <- c(0, 0.5)
#>     })
#>     model({
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc)
#>         Cc <- linCmt()
#>         Cc ~ prop(propSd)
#>     })
#> }

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

mod %>% rmEta(c("etaKa", "etaCl"))
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#>  ── rxode2-based solved PK 1-compartment model ────────────────────────────────── 
#>  ── Initalization: ──  
#> Fixed Effects ($theta): 
#>    lka    lcl    lvc propSd 
#>   0.45   1.00   3.45   0.50 
#> 
#> Omega ($omega): 
#>       etaVc
#> etaVc   0.3
#> 
#> States ($state or $stateDf): 
#>   Compartment Number Compartment Name Rate   Off Internal #
#> 1                  1            depot TRUE FALSE          1
#> 2                  2          central TRUE FALSE          2
#>  ── μ-referencing ($muRefTable): ──  
#>   theta   eta level
#> 1   lvc etaVc    id
#> 
#>  ── Model (Normalized Syntax): ── 
#> function() {
#>     description <- "One compartment PK model with linear clearance"
#>     ini({
#>         lka <- 0.45
#>         lcl <- 1
#>         lvc <- 3.45
#>         propSd <- c(0, 0.5)
#>         etaVc ~ 0.3
#>     })
#>     model({
#>         ka <- exp(lka)
#>         cl <- exp(lcl)
#>         vc <- exp(lvc + etaVc)
#>         Cc <- linCmt()
#>         Cc ~ prop(propSd)
#>     })
#> }
```
