# Second-order eta sensitivities of the prediction for the exact log-likelihood (\`ll()\`/generalized) inner Hessian under \`fast=TRUE\`.

For a generalized endpoint \`rx_pred\_\` is the per-observation
log-density, so its second eta-derivatives \`d2(rx_pred\_)/deta_i
deta_j\` let \`calcEtaHessian\` assemble the exact inner Hessian \`H =
Omega^-1 - sum_obs d2(logLik)/deta2\` analytically – mirroring the
Gaussian Gauss-Newton \`sum(cHff\*a\*a)+Omega^-1\` assembly – instead of
the Shi21 finite difference of the inner gradient. Reuses the
augmented-model second-order chain (\`.g2\`, see
\[.foceiAnalyticAugModelDirs\]). Stores on the symengine env
\`..HdEta2\` (the \`rx\_\_d2pred_i_j\_\_\` lhs lines, upper triangle
i\<=j) and \`..sens2\` (the second-order state-sensitivity ODEs). Only
built for a fast generalized fit; the ordinary inner model is unchanged.

## Usage

``` r
# S3 method for class 'foceiHdEta2'
rxUiGet(x, ...)
```

## Arguments

- x:

  list of rxode2 UI

- ...:

  ignored

## Value

symengine env with \`..HdEta2\` and \`..sens2\` added
