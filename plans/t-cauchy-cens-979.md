# M2/M3/M4 censoring ignored for t()/cauchy() endpoints (#979)

## Root cause

`t()`/`cauchy()` compile through the general-log-likelihood path
(`rx_pred_ ~ llikT(dv, nu, mean, sd)` / `llikCauchy(dv, mean, sd)`).  Every
censoring guard in the C++ layer is `dist == rxDistributionNorm` (or
`Dnorm`); non-normal dist branches take `f` (the already-computed log
density) verbatim, so `CENS`/`LIMIT` never reach the objective for t/cauchy.
`doCensT1` did not exist yet. Cauchy is Student-t with `nu=1`, so one
function serves both.

## Scope (revised during implementation)

- **nlm-family** (`nlm`/`bobyqa`/`newuoa`/`uobyqa`/`n1qn1`/`lbfgsb3c`/
  `optim`/`nlminb`, `src/nlm.cpp`): full M2/M3/M4 support for t()/cauchy(),
  value AND gradient (nlm already forces FD for any censored row's theta
  gradient regardless of distribution). Verified end-to-end against a
  censoring-ignored baseline and against the t(nu=1)==cauchy() identity.
- **FOCEi/FOCE/AGQ/Laplace (`src/inner.cpp`) -- NOT fixed in this PR.**
  Giving the llik-forced t()/cauchy() endpoint a real (nonzero) `rx_r_`
  (needed for the CDF correction) crashes rxode2's eta-sensitivity
  generation once real etas are present ("user function 'get' requires 5
  arguments") -- reproduces with plain `cauchy()`, no censoring at all, and
  is caused by making `rx_r_` nonzero, not by anything censoring-specific.
  `censEst.h`'s `doCensT1()`/`dCensT1()` and `src/inner.cpp`'s
  `focei_tCensLl()` are written and validated (nlm-family exercises the
  identical `doCensT1()` path) but are deliberately kept INERT for FOCEi:
  `.fixCensRNuLine` (R/focei.R) only emits the real `rx_r_`/`rx_nu_` when
  `nlmixr2global$rxCensNuFix` is set, and only `rxUiGet.nlmModel0` sets it.
  A follow-up issue should track the rxode2-side interaction.
- SAEM forces any `ll()`-family endpoint (t/cauchy included) onto
  `distribution==4`, whose M-step (`augmentCensY`) would need a truncated
  Student-t sampler -- real new machinery, out of scope here.
- Every other generalized-likelihood distribution (pois/binom/beta/chisq/
  dexp/f/geom/unif/weibull/dgamma/ordinal/general `ll()`), and t()/cauchy()
  under every method except nlm-family, get a concise runInfo warning that
  censoring is ignored, instead of silently wrong results.

## Phases

1. [x] `src/censEst.h`: `doCensT1()` (M2/M3/M4 value correction via
   `R::pt()`, nu=1 for cauchy) + `dCensT1()` (first derivative). No 2nd/3rd
   order partials -- censored T/Cauchy falls through to the existing
   Gauss-Newton curvature fallback (documented as a KNOWN GAP). Validated
   standalone: nu=1 matches pcauchy exactly, large-nu matches
   doCensNormal1, dCensT1 matches finite differences to <1e-7.
2. [x] R-side plumbing: `.fixCensRNuLine` (R/focei.R) squares `rx_rll_`
   into `rx_r_` and captures `llikT()`/`llikXT()`'s `nu` argument (1 for
   cauchy) into a new `rx_nu_` line, gated to nlm-family only
   (`nlmixr2global$rxCensNuFix`) -- see Scope above for why FOCEi is
   excluded. `.nlmGetFRLines`/`rxUiGet.nlmRxModel` (R/nlm.R) re-splice
   `rx_nu_` from the symengine environment the same way they already do
   for `rx_pred_f_`/`rx_r_` (a lone `rx_nu_ ~ nu` does not survive
   rxOptExpr's dead-code elimination on its own).
3. [x] `src/nlm.cpp`: widened the three `dist ==` guards (objective,
   likelihood-contrib cotangent, gradient) to include T/Cauchy via a new
   name-looked-up `predNuOffset`/`gradNuOffset`; calls `doCensT1`.
4. [x] `src/inner.cpp`: added `predFOffset`/`predROffset`/`predNuOffset`
   and `focei_tCensLl()` (mirrors nlm's value correction) at the three
   `llikObs[kk] = f` sites -- currently inert for FOCEi (see Scope); ready
   for the follow-up. Eta-gradient/curvature stay a documented KNOWN GAP
   even once FOCEi's R-side gap closes (no `d(rx_pred_f_)/d(eta)` column).
5. [x] Catch-all warning: `R/preProcessCensDistWarn.R` hook, runs for every
   estimation method, warns (<75 chars, no method prefix) when CENS/LIMIT
   data meets an endpoint outside {norm, dnorm} (plus {t, cauchy} for
   nlm-family only). Skips `est="nls"` (already hard stops).
6. [x] Tests: `test-nlm-cens-t.R` (nlm t()/cauchy() censoring vs a
   censoring-ignored baseline, censInformation asserts the mechanism ran,
   t(nu=1)==cauchy() identity, no spurious warning) and
   `test-cens-dist-warn.R` (pois()+nlm warns, t()+focei warns but
   t()+nlm does not, message-length check). Both essential-tier (small,
   fast fits) -- not added to `.slowBatches`.
7. [x] NEWS.md bullet under `## Bug fixes`, referencing #979.
8. [ ] Antigravity review until clean.
9. [ ] Push, open PR closing #979.
