# Finalizes output list

Finalizes output list

## Usage

``` r
.nlmFinalizeList(env, lst, par = "par", printLine = TRUE, hessianCov = TRUE)
```

## Arguments

- env:

  nlm environment

- lst:

  output list

- par:

  parameter name of final estimate in output

- printLine:

  Print the final line when print is nonzero

- hessianCov:

  boolean indicating a hessian should be used/calculated for covariance

## Value

modified list with \`\$cov\`

## Author

Matthew L. Fidler
