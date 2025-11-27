# Respect suppress messages for nlmixr2 C functions

This turns on the silent REprintf in C when \`suppressMessages()\` is
turned on. This makes the \`REprintf\` act like \`messages\` in R, they
can be suppressed with \`suppressMessages()\`

## Usage

``` r
nmSuppressMsg()
```

## Value

Nothing

## Author

Matthew Fidler

## Examples

``` r
# nmSupressMsg() is called with nlmixr2()

# In nlmixr2, we use REprintf so that interrupted threads do not crash R
# if there is a user interrupt. This isn't captured by R's messages, but
# This interface allows the `suppressMessages()` to suppress the C printing
# as well

# If you  want to suppress messages from nlmixr2 in other packages, you can use
# this function
```
