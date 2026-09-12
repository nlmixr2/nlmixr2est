# nlmixr2est 7.0.3

## Requires rxode2 5.1.7

`DESCRIPTION` requires `rxode2 (>= 5.1.7)`.  The version currently published on
CRAN is 5.1.6, so **this submission should be reviewed after `rxode2` 5.1.7**.
This differs from the 7.0.2 submission, which was deliberately written to build
and run against either of two `rxode2` versions; 7.0.3 uses entry points that
5.1.7 adds and does not carry a fallback for them.

## Test environments

* local: Ubuntu 24.04, R 4.6.1 (x86_64-pc-linux-gnu), gcc/g++ 14.2.0
* GitHub Actions: Windows, macOS and Linux, R release and devel

## R CMD check results

`R CMD check --as-cran` on the submitted tarball: **Status: 1 NOTE**.

```
* checking compilation flags used ... NOTE
Compilation used the following non-portable flag(s):
  '-mno-omit-leaf-frame-pointer'
```

That flag is not set by this package.  It comes from the R installation used for
the check (it is in that platform's `Makeconf` `CFLAGS`), so it reflects the local
build environment rather than anything in `src/Makevars`.  We expect it not to
appear on CRAN's builders.

There are no other NOTEs, WARNINGs or ERRORs.  `checking CRAN incoming
feasibility`, `checking compiled code`, `checking examples`, `checking examples
with --run-donttest`, `checking tests`, `checking package vignettes` and
`checking re-building of vignette outputs` are all OK.

## Reverse dependencies

All 16 CRAN reverse dependencies were checked against this version: `admixr2`,
`babelmixr2`, `ggPMX`, `nlmixr2`, `nlmixr2auto`, `nlmixr2autoinit`,
`nlmixr2extra`, `nlmixr2lib`, `nlmixr2plot`, `nlmixr2rpt`, `nlmixr2save`,
`nlmixr2targets`, `nlmixr2utils`, `shinyMixR`, `xpose.nlmixr2` and
`xpose.xtras`.

None of them broke.  The few NOTEs and WARNINGs seen are properties of those
packages' own sources or of how our harness invoked the check (it builds and
checks without vignettes), not of this update.
