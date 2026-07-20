# Preprocessing hook: mu-expand non-mu structural fixed-effect thetas for the nonparametric engines, before the rest of the pipeline builds on the ui. Doing it here (rather than mutating the ui mid-setup) keeps the injected pseudo-etas consistent through covariate/mu processing and model compilation. Off when control\$muExpand is FALSE.

Preprocessing hook: mu-expand non-mu structural fixed-effect thetas for
the nonparametric engines, before the rest of the pipeline builds on the
ui. Doing it here (rather than mutating the ui mid-setup) keeps the
injected pseudo-etas consistent through covariate/mu processing and
model compilation. Off when control\$muExpand is FALSE.

## Usage

``` r
.nlmixr0preProcessNpMuExpand(ui, est, data, control)
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

list(ui=) when the model was expanded, else NULL

## Author

Matthew L. Fidler
