# nlmixr2est 7.0.2

## Submitted together with rxode2 5.1.6

`rxode2` 5.1.6 is being submitted at the same time as this package.  **This
submission does not depend on it.**  `nlmixr2est` builds and runs against the
currently published `rxode2` 5.1.5 as well as 5.1.6, and `DESCRIPTION` requires
only `rxode2 (>= 5.1.5)`, so the two can be reviewed and published in either
order.

`rxode2` 5.1.6 adds entry points to its linked C function-pointer table.  Where
they are present this package uses them; where they are not it uses the code it
shipped in 7.0.1.  Which path is taken is decided twice, deliberately:

* at compile time, by `configure` probing the installed `rxode2ptr.h`, so the
  newer calls are only compiled when the header declares them; and
* at run time, by checking the function pointer is non-`NULL` before it is
  called, because a binary compiled against 5.1.6 may later be loaded against
  5.1.5, where those table slots are absent.

We verified both directions on our own machines: the package builds and its test
suite passes against 5.1.5 and against 5.1.6, with identical numerical results.

## Previous Windows lazy-loading failure

The 7.0.1 submission failed on the Windows builder while preparing the package
for lazy loading, in `stringfish`:

```
unable to load shared object '.../stringfish/libs/x64/stringfish.dll':
  LoadLibrary failure:  Das angegebene Modul wurde nicht gefunden.
ERROR: lazy loading failed for package 'nlmixr2est'
```

That was not this package's code -- `stringfish` was reached through
`rxode2` -> `qs2` -> `stringfish`, and the failure happened before any of our
code ran.  It is resolved rather than worked around: the `qs2` dependency (and
with it `stringfish`) has been **dropped**.  The compiled-model disk cache now
uses RDS files, and compressed fit components use base R serialization.  Neither
`qs2` nor `stringfish` appears in `DESCRIPTION` any more.

## Test environments

* local: Ubuntu 24.04, R 4.6.1 (x86_64-pc-linux-gnu)
* GitHub Actions: Windows, macOS and Linux, R release and devel

## R CMD check results

<!-- filled in from the --as-cran run before submission -->

## Reverse dependencies

We checked the reverse dependencies of `nlmixr2est`; see
`revdep/` for the results.
