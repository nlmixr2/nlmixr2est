# FOCEi conditional-likelihood C API pointers

Returns the small external-pointer table (apiVersion, dims, setTheta,
condBatch, setOmegaInv, thetaSensIdx, condThetaGrad) over the
\[foceiLikLoad()\]-ed problem that downstream packages (e.g.
nlmixr2stan) install via \`inst/include/nlmixr2estFoceiPtr.h\`, which
also documents the calling contract (#937).

## Usage

``` r
.nlmixr2estFoceiPtrs()
```

## Value

a named list of external pointers.
