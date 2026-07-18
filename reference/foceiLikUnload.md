# Unload the general FOCE-family likelihood from memory

Frees the FOCEi inner problem set up by \[foceiLikLoad()\]. A no-op
(returns \`FALSE\`) if nothing is loaded.

## Usage

``` r
foceiLikUnload()
```

## Value

Invisibly \`TRUE\` if a system was freed, \`FALSE\` if none was loaded.

## See also

\[foceiLikLoad()\], \[foceiLikRun()\]

## Author

Matthew L. Fidler
