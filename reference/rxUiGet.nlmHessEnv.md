# Build the symengine env with second-order theta sensitivities + hess lines

Starts from \`rxUiGet.nlmHdTheta\` (the same forward first-order env the
thetaGrad model uses; never the adjoint substitution in
\`rxUiGet.nlmEnv\`). Since \`rx_pred\_\` is the per-observation -llik,
the emitted \`rx_hess\_\<i\>\_\<j\>\_\` columns (i\<=j triangle, j
outer/i inner) are the objective-Hessian contributions directly.

## Usage

``` r
# S3 method for class 'nlmHessEnv'
rxUiGet(x, ...)
```

## Value

symengine env with \`..sens2\` and \`..hessLines\`, or a character
reason when out of scope
