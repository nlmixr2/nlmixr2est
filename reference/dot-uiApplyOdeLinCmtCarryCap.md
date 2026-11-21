# Solve a \`linCmt()\` model past the carry cap as ODEs

The eta-gradient carry has \`.rxFoceiLinCmtCarryMaxPairs()\` columns.
The etas past them would keep the standard gradient, which is wrong when
a covariate changes a \`linCmt()\` parameter within a subject and makes
the objective depend on the order the etas are declared in (#1318). The
ODE translation has exact sensitivities for every eta. Named to run
after \`.uiApplyIov()\` and \`.uiApplyMu2()\`, so it sees the model the
build does.

## Usage

``` r
.uiApplyOdeLinCmtCarryCap(ui, est, data, control)
```

## Arguments

- ui:

  rxode2 ui

- est:

  estimation method (all methods are shown by \`nlmixr2AllEst()\`).
  Methods can be added for other tools

- data:

  nlmixr data

- control:

  The estimation control object. These are expected to be different for
  each type of estimation method

## Value

list with the translated ui and renumbered data, or \`NULL\` when
nothing to do

## Author

Matthew L. Fidler
