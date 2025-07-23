# nlmixr2est 4.0.2

- The loading and unloading of DLLs has been minimized in this version
  of nlmixr2est. This avoids loading/reloading the same DLLs and causing the
  CRAN mac m1 ASAN/USBAN false positive issue observed in CRAN.

- Additionally a new function `nlmixr2fix(fit)` has been added to
 `nlmixr2est`.  It attempts to make the fit loaded from a different
 version of nlmixr2 compatible with nlmixr2 4.0.  It also prints out
 the versions of `nlmixr2` that were used when creating this fit.
 With this information you are more likely to find a way to use the
 fit in your current session (or in an old session). (Issue #562)
