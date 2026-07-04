# Collect warnings and just warn once.

Collect warnings and just warn once.

## Usage

``` r
.collectWarn(expr, lst = FALSE, collectErr = FALSE)
```

## Arguments

- expr:

  R expression

- lst:

  When `TRUE` return `list(object, warning = ws, error = es)` instead of
  issuing the warnings.

- collectErr:

  When `TRUE`, also record errors raised during evaluation instead of
  letting them propagate; used by `nlmixr2Est0()` so all errors from a
  failed run are reported together rather than only the last one.

## Value

The value of the expression, or when `lst = TRUE` a list
`list(object, warning = ws, error = es)` of unique warning/error
messages (`es` is `NULL` unless `collectErr = TRUE` and an error
escaped).

## Author

Matthew L. Fidler
