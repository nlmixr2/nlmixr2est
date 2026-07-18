# VAE preprocessing hook: inject etas for non-mu-referenced thetas per vaeControl(nonMuTheta=).

VAE preprocessing hook: inject etas for non-mu-referenced thetas per
vaeControl(nonMuTheta=).

## Usage

``` r
.preProcessVaeNonMuTheta(ui, est, data, control)
```

## Arguments

- est:

  estimation method (all methods are shown by \`nlmixr2AllEst()\`).
  Methods can be added for other tools

- data:

  nlmixr data

- control:

  The estimation control object. These are expected to be different for
  each type of estimation method

## Value

list(ui=) possibly with injected etas

## Author

Matthew L. Fidler
