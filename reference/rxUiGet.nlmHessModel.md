# Compile the analytic Hessian model for the nlm family

Compile the analytic Hessian model for the nlm family

## Usage

``` r
# S3 method for class 'nlmHessModel'
rxUiGet(x, ...)
```

## Value

\`list(thetaHess=\<rxode2 model\>, ntri=\<lhs triangle length\>)\` when
the model can be built; \`list(thetaHess=NULL, reason=\<why\>)\`
otherwise (feeds \`.nlmHessAnalyticFallback\`)
