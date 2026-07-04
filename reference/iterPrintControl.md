# Control iteration-time print formatting

Bundles the options controlling the iteration progress output emitted by
\`nlmixr2\` estimators. Pass as the \`print\` argument to any
\`\*Control()\` function; the scalar form (\`print = N\`) still works
and is wrapped into an \`iterPrintControl()\` internally.

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
  evaluations; \`0\` suppresses output. Defaults to \`1L\`.

- ncol:

  Integer or \`NULL\`. Parameter columns per row before wrapping.
  \`NULL\` (default) uses \`floor((getOption("width") - 23) / 12)\`.

- headerEvery:

  Integer or \`NULL\`. Re-emit the column header every \`headerEvery\`
  parameter-print events; \`0\` prints it once at fit start. \`NULL\`
  (default) uses \`10L\`.

- useColor:

  Logical (or \`NULL\`) emit ANSI bold/color escapes in the iteration
  print. \`NULL\` (default) defers to \[crayon::has_color()\].

- simple:

  Logical. When \`TRUE\`, print a single row per iteration, suppressing
  the unscaled (\`U\`) / back-transformed (\`X\`) rows. Defaults to
  \`FALSE\`.

## Value

A list with the validated, defaulted iteration-print options. Has class
\`"iterPrintControl"\` so the outer \`\*Control()\` functions can
distinguish a pre-built object from a scalar \`print = N\`.

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
