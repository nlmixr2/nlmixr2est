# Iteration-print configuration parameters (documentation stub)

Documentation-only stub. This is the single canonical location for every
argument that controls iteration-time print formatting in
\`nlmixr2est\`. Two parallel sets of names live here:

## Arguments

- print:

  Either a scalar print-frequency (\`0\` = suppress iteration output;
  \`1\` (default) = print every parameter evaluation; \`N\` = print
  every Nth evaluation), OR a pre-built \[iterPrintControl()\] object
  bundling all iteration-print options (column wrap, header cadence,
  color, simple/three-row mode). The scalar form is equivalent to
  \`iterPrintControl(every = print, ncol = printNcol, useColor =
  useColor)\`.

- printNcol:

  Integer (or \`NULL\`) — number of parameter columns emitted per row
  before wrapping to a continuation row. \`NULL\` (the default) defers
  to \[iterPrintControl()\]'s default (\`floor((getOption("width") - 23)
  / 12)\`, which fits an 80-column terminal). Equivalent to \`print =
  iterPrintControl(ncol = ...)\`.

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

Nothing; this is a documentation-only helper.

## Details

\- The scalar arguments accepted by every \`\*Control()\` function
(\`print\`, \`printNcol\`, \`useColor\`), and - The arguments of
\[iterPrintControl()\] itself (\`every\`, \`ncol\`, \`headerEvery\`,
\`useColor\`, \`simple\`).

Both sets are documented here so a reader looking up any one argument
finds all of them — and downstream functions pull just the entries they
need via \`@inheritParams iterPrintParams\`.
