# Wrap scalar or list arguments into an iterPrintControl object

Package-set internal helper used by every \`\*Control()\` function in
\`nlmixr2est\` and downstream packages (e.g. \`babelmixr2\`) to absorb
the scalar \`print\` / \`printNcol\` / \`useColor\` arguments into a
single \[iterPrintControl()\] sub-list. If the user already passed a
pre-built \`iterPrintControl()\` object via the \`print\` argument (or,
on round-trip, via an \`iterPrintControl =\` slot in \`...\`), return it
directly.

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
  and the other scalars when supplied — used by the round-trip case
  where a returned control list is passed back through
  \`do.call(\*Control, .ctl)\`.

## Value

An \`iterPrintControl\` list.

## Details

Exported under a leading-dot name to mark it as a package-set internal —
callable from sibling packages in the nlmixr2 family but not meant for
end users.
