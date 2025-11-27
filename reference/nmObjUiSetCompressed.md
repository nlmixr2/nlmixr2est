# Set if the nlmixr2 object will return a compressed ui

Set if the nlmixr2 object will return a compressed ui

## Usage

``` r
nmObjUiSetCompressed(type)
```

## Arguments

- type:

  is a boolean indicating if the compressed ui will be returned
  (\`TRUE\`) or not be returned (\`FALSE\`)

## Value

invisible logical type

## Author

Matthew L. Fidler

## Examples

``` r
nmObjUiSetCompressed(FALSE) # now the $ui will return an environment
nmObjUiSetCompressed(TRUE) # now the $ui will return a compressed value
```
