# The nlm population-likelihood C entry-point table (#953)

Returns the external-pointer table a downstream package installs with
\`iniNlmixr2estNlm()\` from \`inst/include/nlmixr2estNlmPtr.h\`:
apiVersion, dims, and the value + analytic theta-gradient evaluation
over the \[nlmObjectiveSetup()\]-loaded problem.

## Usage

``` r
.nlmixr2estNlmPtrs()
```

## Value

a named list of external pointers.
