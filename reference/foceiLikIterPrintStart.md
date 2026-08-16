# Iteration printing / parameter history over a loaded likelihood

The standard nlmixr2est iteration table (the shared \`scale.h\`
machinery behind every estimation method's printout) exposed for an
EXTERNAL sampler driving \[foceiLikLoad()\]'s conditional likelihood –
e.g. nlmixr2stan's \`est="stan"\`. \`foceiLikIterPrintStart()\` prints
the header and arms the row entry (entry 8 of the FOCEi C API table,
\`nlmixr2FoceiIterPrintRow\`), which the sampler's gradient evaluations
call from C; every \`every\`-th call prints a row AND records it, so the
returned history holds exactly the printed rows.
\`foceiLikIterPrintEnd()\` prints the closing line and returns the
history as a \`parHistData\`-style data frame.

## Usage

``` r
foceiLikIterPrintStart(
  every,
  initPar,
  names,
  iterPrintControl = NULL,
  xform = NULL
)

foceiLikIterPrintEnd()
```

## Arguments

- every:

  print/record cadence in row-entry calls (0 arms nothing; the row entry
  then records nothing and returns immediately)

- initPar:

  display vector at the starting point (defines the width)

- names:

  column names, same length as \`initPar\`

- iterPrintControl:

  optional list applied via the shared \`scaleApplyIterPrintControl\`
  (e.g. \`list(useColor=FALSE, printNcol=6L)\`)

- xform:

  optional back-transform list (the \`.iterPrintXParFromUi\` shape)
  driving the \`X\` row

## Value

\`foceiLikIterPrintStart()\`: invisibly \`NULL\`;
\`foceiLikIterPrintEnd()\`: the recorded history data frame (or \`NULL\`
when printing was never started)

## Details

The display vector is CALLER-defined, not the internal focei parameter
vector: pass natural-scale thetas plus (optionally) the sampler's
current actual omega entries, named \`om.\<eta\>\` for variances and
\`cov.\<eta1\>.\<eta2\>\` for covariances. (The internal vector's omega
tail is \`chol(Omega^-1)\` in the \`diagXform\` parameterization, which
a Bayesian sampler neither uses nor understands, so it is NOT printed.)

## Author

Matthew L. Fidler
