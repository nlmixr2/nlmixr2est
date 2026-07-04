# Iteration-print configuration parameters (documentation stub)

Shared \`@param\` docs for iteration-print formatting, used by both the
scalar \`\*Control()\` arguments and \[iterPrintControl()\] via
\`@inheritParams\`.

## Arguments

- print:

  Either a scalar print-frequency (\`0\` = suppress, \`1\` (default) =
  every evaluation, \`N\` = every Nth), OR a pre-built
  \[iterPrintControl()\] object. Equivalent to \`iterPrintControl(every
  = print, ncol = printNcol, useColor = useColor)\`.

- printNcol:

  Integer (or \`NULL\`) parameter columns per row before wrapping.
  \`NULL\` (default) uses \`floor((getOption("width") - 23) / 12)\`.

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

Nothing; this is a documentation-only helper.
