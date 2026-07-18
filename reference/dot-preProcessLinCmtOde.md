# Translate a mixed \`linCmt()\`/ODE model to all-ODEs

rxode2 requires the \`linCmt()\` compartments to be the last states of
the solve (\`op\$linOffset = neq - numLin - numLinSens\`), so their
compartment number is one past the ODE states. The FOCEi inner model
adds an ODE state per eta, and the nlm family one per theta, which
pushes depot/central past the compartment numbers the data was
translated against (\`.foceiPreProcessData()\` uses the plain model). A
dose then silently lands in a sensitivity state – every prediction comes
back 0. Solving the linear part as ODEs removes the \`linCmt()\` block,
so the numbering agrees again. SAEM builds no sensitivity states and is
left alone.

## Usage

``` r
.preProcessLinCmtOde(ui, est, data, control)
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

list with the translated ui, or \`NULL\` when nothing to do

## Author

Matthew L. Fidler
