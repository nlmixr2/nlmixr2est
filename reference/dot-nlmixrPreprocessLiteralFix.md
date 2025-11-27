# This literally fixes parameters in the model

Whenever there is a fixed parameter in the model, the parameter is
replaced with the literal value inside of the model and dropped from the
\`ini\` block. This only occurs when the \`control\$literalFix=TRUE\`.

## Usage

``` r
.nlmixrPreprocessLiteralFix(ui, est, data, control)
```

## Arguments

- ui:

  model function/object

- est:

  estimation method (all methods are shown by \`nlmixr2AllEst()\`).
  Methods can be added for other tools

- data:

  nlmixr data

- control:

  The estimation control object. These are expected to be different for
  each type of estimation method

## Value

list with possibly updated ui

## Author

Matthew L. Fidler
