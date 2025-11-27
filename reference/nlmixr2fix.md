# Try to fix a nlmixr2 fit

Currently this re-evaluates the function in the current version of
rxode2.

## Usage

``` r
nlmixr2fix(fit)
```

## Arguments

- fit:

  nlmixr2 fit object from a different version of nlmixr2.

## Value

A nlmixr2 fit that has been (possibly) adjusted to work with the current
version of nlmixr2.

## Author

Matthew L. Fidler

## Examples

``` r
# \donttest{

# This is a nlmixr2 v3 fit
fit <- system.file("testfit_nlmixr3.rds", package = "nlmixr2est")
fit <- readRDS(fit)

# While it prints well, it can't be used in all functions because
# Language features (like +var()) are not supported in the v3 version

print(fit)
#> Error in list2env(ui, parent = emptyenv()): first argument must be a named list

try(rxSolve(fit)) # should error, but with try it will just display the error
#> Error in list2env(ui, parent = emptyenv()) : 
#>   first argument must be a named list

# This function attempts to fix it by regenerating the rxode2 model with the
# new features

# This function also prints out the information on how this fit was created

fit <- nlmixr2fix(fit)
#> # This function is meant to load nlmixr2 fits from other versions
#> # To reproduce the fit, you need to use the same version of nlmixr2
#> ## ==============================
#> ## nlmixr2est Session Information
#> ## ==============================
#> # OS: Debian GNU/Linux 12 (bookworm)
#> # LAPACK: /cosbi/home/boaretti/miniconda3/envs/r_env3/lib/libopenblasp-r0.3.29.so
#> # LAPACK Version: 3.12.0
#> # R Version: R version 4.3.3 (2024-02-29)
#> 
#> 
#> # Install Package version 1.3.1-13 of 'dparser'
#> remotes::install_version("dparser", version = "1.3.1-13") 
#> 
#> # Install Package version 1.0.0 of 'lotri'
#> remotes::install_version("lotri", version = "1.0.0") 
#> 
#> # Install Package version 0.7 of 'PreciseSums'
#> remotes::install_version("PreciseSums", version = "0.7") 
#> 
#> # Install Package version 2.0.11 of 'rxode2ll'
#> remotes::install_version("rxode2ll", version = "2.0.11") 
#> 
#> # Install Package version 3.0.4 of 'rxode2'
#> remotes::install_version("rxode2", version = "3.0.4") 
#> 
#> # Install Package version 2024-3.5 of 'lbfgsb3c'
#> remotes::install_version("lbfgsb3c", version = "2024-3.5") 
#> 
#> # Install Package version 6.0.1-12 of 'n1qn1'
#> remotes::install_version("n1qn1", version = "6.0.1-12") 
#> 
#> # Install Package version 3.0.4 of 'nlmixr2est'
#> remotes::install_version("nlmixr2est", version = "3.0.4") 
#> 
#> # Install Package version 3.0.0 of 'nlmixr2extra'
#> remotes::install_version("nlmixr2extra", version = "3.0.0") 
#> 
#> # Package 'nlmixr2lib' is not installed, but known to enhance nlmixr2/babelmixr2
#> 
#> # Install Package version 3.0.1 of 'nlmixr2'
#> remotes::install_version("nlmixr2", version = "3.0.1") 
#> 
#> # Package 'nonemem2rx' is not installed, but known to enhance nlmixr2/babelmixr2
#> 
#> # Package 'monolix2rx' is not installed, but known to enhance nlmixr2/babelmixr2
#> 
#> # Package 'babelmixr2' is not installed, but known to enhance nlmixr2/babelmixr2
#> 
#> # Package 'PopED' is not installed, but known to enhance nlmixr2/babelmixr2
#> 
#> # Install Package version 0.11.0 of 'PKNCA'
#> remotes::install_version("PKNCA", version = "0.11.0") 
#> 
#> # If all else fails you can try to install the version of nlmixr2 used to create the fit
#> Error in list2env(ui, parent = emptyenv()): first argument must be a named list

# Now solving and other functions work

rxSolve(fit)
#> Error in list2env(ui, parent = emptyenv()): first argument must be a named list

# }
```
