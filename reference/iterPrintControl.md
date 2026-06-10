# Control iteration-time print formatting

Bundles every option that controls the formatting of the iteration
progress output emitted by \`nlmixr2\` estimators (focei, saem, bobyqa,
nlm, optim, nls, nlminb, lbfgsb3c, n1qn1, newuoa, uobyqa). The returned
list is consumed by the shared C++ helper
(\`scaleApplyIterPrintControl\` in \`src/scale.h\`) so every estimator
formats its iteration trace through one code path.

## Usage

``` r
iterPrintControl(
  every = 1L,
  ncol = NULL,
  headerEvery = NULL,
  useColor = NULL,
  simple = FALSE
)
```

## Arguments

- every:

  Integer. Print one iteration row every \`every\` parameter
  evaluations. \`0\` suppresses iteration output entirely. Defaults to
  \`1L\`. Inside \[iterPrintControl()\] this is the canonical name for
  what the outer \`\*Control()\` functions call \`print\`.

- ncol:

  Integer or \`NULL\`. Number of parameter columns to emit per row
  before wrapping to a continuation row. \`NULL\` (default) uses
  \`floor((getOption("width") - 23) / 12)\`, which fits an 80-column
  terminal. Inside \[iterPrintControl()\] this is the canonical name for
  what the outer \`\*Control()\` functions call \`printNcol\`.

- headerEvery:

  Integer or \`NULL\`. Re-emit the column header every \`headerEvery\`
  parameter-print events (not raw iterations). With \`every = 5\` and
  \`headerEvery = 10\`, the header re-prints every 50 iterations. \`0\`
  prints the header once at fit start only. \`NULL\` (default) uses
  \`10L\`.

- useColor:

  Logical (or \`NULL\`) — whether to emit ANSI bold/color escapes in the
  iteration print. \`NULL\` (the default) defers to
  \[iterPrintControl()\]'s default (\[crayon::has_color()\]). Equivalent
  to \`print = iterPrintControl(useColor = ...)\`.

- simple:

  Logical. When \`TRUE\`, the printer emits a single row per iteration
  (just the optimizer-scale parameters) and suppresses the unscaled
  (\`U\`) / back-transformed (\`X\`) follow-up rows. Used by estimators
  (like saem) that have no internal optimizer scaling, where U and X
  would be degenerate copies of the first row. Defaults to \`FALSE\`
  (full three-row output).

## Value

A list with the validated, defaulted iteration-print options. Has class
\`"iterPrintControl"\` so the outer \`\*Control()\` functions can
distinguish a pre-built object from a scalar \`print = N\`.

## Details

Pass this as the \`print\` argument to any of the \`\*Control()\`
functions, e.g. \`foceiControl(print = iterPrintControl(every = 5,
headerEvery = 20))\`. Equivalently, the existing scalar form
\`foceiControl(print = 5, printNcol = 8)\` continues to work — the outer
\`\*Control()\` wraps your individual arguments into an
\`iterPrintControl()\` call internally.

All argument descriptions live in \[iterPrintParams\] so the same text
is shared with every \`\*Control()\` function.

## Author

Bill Denney, Matthew L. Fidler

## Examples

``` r
iterPrintControl()
#> $every
#> [1] 1
#> 
#> $ncol
#> [1] 4
#> 
#> $headerEvery
#> [1] 10
#> 
#> $useColor
#> [1] TRUE
#> 
#> $simple
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "iterPrintControl" "list"            
iterPrintControl(every = 5, headerEvery = 0)
#> $every
#> [1] 5
#> 
#> $ncol
#> [1] 4
#> 
#> $headerEvery
#> [1] 0
#> 
#> $useColor
#> [1] TRUE
#> 
#> $simple
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "iterPrintControl" "list"            
```
