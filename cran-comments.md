# nlmixr2est 7.2.0

This release fixes the r-devel-linux-x86_64-debian-clang ERROR in the CRAN checks
of 7.1.0 (a segfault in `rxOptExpr()` during the tests).  The cause was a data
race in 'dparser' when 'rxode2' parsed in parallel; 'rxode2' 5.1.8 (now on CRAN)
only parses in parallel when 'dparser' is race-free, and this version requires
`rxode2 (>= 5.1.8)`.  The failing test file was run repeatedly with a clang 21 /
libomp 21 build against the CRAN 'rxode2' and 'dparser' without a crash.

## Test environments

* local: Ubuntu 24.04, R 4.6.1 (x86_64-pc-linux-gnu), gcc/g++ 14.2.0, checked
  against the CRAN versions of all dependencies (including `rxode2` 5.1.8)
* local: clang 21.1.8 + libomp 21 build of the package and of the CRAN
  `rxode2`/`dparser`
* GitHub Actions: Windows, macOS and Linux, R release and devel

## R CMD check results

`R CMD check --as-cran` on the submitted tarball: 0 ERRORs, 0 WARNINGs, 1 NOTE:

* `checking compilation flags used ... NOTE` (`-mno-omit-leaf-frame-pointer`):
  set by the local R installation's `Makeconf`, not by `src/Makevars`.

## Reverse dependencies

All 18 CRAN reverse dependencies were checked (CRAN versions, CRAN-mode
`R CMD check`) against this version: `admixr2`, `babelmixr2`, `ggPMX`,
`nlmixr2`, `nlmixr2auto`, `nlmixr2autoinit`, `nlmixr2extra`, `nlmixr2lib`,
`nlmixr2plot`, `nlmixr2rpt`, `nlmixr2save`, `nlmixr2scm`, `nlmixr2targets`,
`nlmixr2utils`, `nonmem2rx`, `shinyMixR`, `xpose.nlmixr2` and `xpose.xtras`.
All returned `Status: OK`.
