# Wrap scalar or list arguments into an iterPrintControl object

Absorbs the scalar \`print\`/\`printNcol\`/\`useColor\` arguments into a
single \[iterPrintControl()\] object, or passes through an already
pre-built one.

## Usage

``` r
.absorbIterPrintControl(
  print = 1L,
  printNcol = NULL,
  useColor = NULL,
  iterPrintControl = NULL
)
```

## Arguments

- print:

  Either an integer print-frequency or an \`iterPrintControl\` object.

- printNcol, useColor:

  Scalar \`\*Control()\` arguments forwarded to \[iterPrintControl()\]
  only when \`print\` is a scalar.

- iterPrintControl:

  Optional pre-built \[iterPrintControl()\] object. Wins over \`print\`
  and the other scalars when supplied.

## Value

An \`iterPrintControl\` list.
