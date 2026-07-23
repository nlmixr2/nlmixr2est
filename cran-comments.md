# nlmixr2est 7.0.1

Why are we resubmitting -- There was a feature that we sped up vae covariate selection
1648.9x.  Since CRAN is still not in a working state, and this improvement is so vast
we are resubmitting with a different number

## Windows lazy-loading failure on the incoming checks

A submission failed on the Windows builder while byte-compiling and preparing
the package for lazy loading:

```
Error in inDL(x, as.logical(local), as.logical(now), ...) :
  unable to load shared object '.../stringfish/libs/x64/stringfish.dll':
  LoadLibrary failure:  Das angegebene Modul wurde nicht gefunden.
Calls: <Anonymous> ... asNamespace -> loadNamespace -> library.dynam -> dyn.load
ERROR: lazy loading failed for package 'nlmixr2est'
```

This is a failure to load the `stringfish` DLL on the check machine rather than
an error in this package's code.  `nlmixr2est` does not import `qs2` or
`stringfish` (`qs2` is only in Suggests here); they are reached through
`rxode2`, which imports `qs2`, which imports `stringfish`.  Loading the
`nlmixr2est` namespace therefore loads `stringfish`, and the failure happens
before any of this package's code runs.

Both `stringfish` and `qs2` link against `RcppParallel`; "Das angegebene Modul
wurde nicht gefunden" indicates a dependent DLL of `stringfish.dll` could not be
resolved on that machine.  We have not been able to reproduce it on our own
Windows checks.
