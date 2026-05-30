# This function is to set the number of threads to 2

In general it is a CRAN requirement that packages not use more than 2
threads. This function is to set the number of threads to 2 for CRAN
testing. It is not intended for general use.

## Usage

``` r
aaaCranNlmixrThreads()
```

## Value

nothing, called for side effect of setting the number of threads to 2
for CRAN testing

## Details

When testing with devtools::test() or testthat::test_package(), the
NOT_CRAN environment variable is set to "true", so the number of threads
will not be limited to 2.

## Author

Matthew L. Fidler

## Examples

``` r

# Set the number of threads to 2 for CRAN testing
aaaCranNlmixrThreads()
```
