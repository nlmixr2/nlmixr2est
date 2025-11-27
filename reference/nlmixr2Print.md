# Print x using the message facility

This allows the suppressMessages to work on print functions. This
captures the output function sends it through the message routine.

## Usage

``` r
nlmixr2Print(x, ...)
```

## Arguments

- x:

  object to print

- ...:

  Other things output

## Value

Nothing, called for its side effects

## Details

catpureOutput was used since it is much faster than the internal
capture.output see
https://www.r-bloggers.com/performance-captureoutput-is-much-faster-than-capture-output/

## Author

Matthew L. Fidler
