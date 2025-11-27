# Get an option for the estimation method

By default it gets the focei models if available.

## Usage

``` r
nmObjGetRxSolve(x, what)

# Default S3 method
nmObjGetRxSolve(x, what)
```

## Arguments

- x:

  nlmixr fit object in a list. The class is the estimation method used.

- what:

  What part of the rx solve are you attempting to get?

## Value

The estimation option based on \`what\`, for example
\`nlmixrObjGetRxSolve(x, "atol")\` will get a double vector of absolute
tolerances
