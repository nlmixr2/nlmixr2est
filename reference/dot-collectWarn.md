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

  When `TRUE` return a list with
  `list(object, warning = ws, error = es)` instead of issuing the
  warnings. Otherwise, when `FALSE` issue the warnings and return the
  object.

- collectErr:

  When `TRUE`, errors raised during evaluation of `expr` are recorded in
  addition to warnings. A calling handler captures every error message
  as it is signalled, then lets the condition continue to propagate so
  that inner
  [`try()`](https://rdrr.io/r/base/try.html)/[`tryCatch()`](https://rdrr.io/r/base/conditions.html)
  blocks in `expr` keep working as usual. An outer
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html) catches errors
  that escape all inner handling, so the call always returns instead of
  stopping. If `expr` evaluated to a result (i.e. no error escaped), the
  recorded errors were caught by inner handlers and are discarded. If an
  error did escape, every message observed along the error chain
  (including follow-up errors raised by `on.exit` handlers, see
  issue 607) is returned in the `error` element of the result list (when
  `lst = TRUE`) or re-raised joined by newlines (when `lst = FALSE`).
  When `FALSE` (the default) errors propagate normally and `es` is
  always `NULL`. This is used by `nlmixr2Est0()` so that all errors from
  a failed estimation run are reported together rather than only the
  last one.

## Value

The value of the expression, or when `lst = TRUE` a list of the form
`list(object, warning = ws, error = es)` where `ws` and `es` are
character vectors of unique warning and error messages (`es` is always
`NULL` when `collectErr = FALSE`, and is also `NULL` when the expression
evaluated successfully under `collectErr = TRUE`).

## Author

Matthew L. Fidler
