# Get the ipred model for a fit object depending on the object type

By default it gets the focei models if available.

## Usage

``` r
nmObjGetIpredModel(x)

# S3 method for class 'saem'
nmObjGetIpredModel(x)

# Default S3 method
nmObjGetIpredModel(x)

# S3 method for class 'saem'
nmObjGetEstimationModel(x)

# Default S3 method
nmObjGetEstimationModel(x)
```

## Arguments

- x:

  nlmixr fit object

## Value

ipred \`rxode2\` model
