# M2/M3/M4 censoring ignored for t()/cauchy() endpoints (#979)

## Root cause

`t()`/`cauchy()` compile through the general-log-likelihood path
(`rx_pred_ ~ llikT(dv, nu, mean, sd)` / `llikCauchy(dv, mean, sd)`).  Every
censoring guard in the C++ layer is `dist == rxDistributionNorm` (or
`Dnorm`); non-normal dist branches take `f` (the already-computed log
density) verbatim, so `CENS`/`LIMIT` never reach the objective for t/cauchy.
`doCensT1` does not exist yet. Cauchy is Student-t with `nu=1`, so one
function serves both.

## Scope

- FOCEI/FOCE/Laplace/AGQ (`src/inner.cpp`) and nlm-family (`src/nlm.cpp`)
  get real M2/M3/M4 support for t()/cauchy().
- SAEM forces any `ll()`-family endpoint (t/cauchy included) onto
  `distribution==4`, whose M-step (`augmentCensY`) would need a truncated
  Student-t sampler -- real new machinery, out of scope here.  SAEM falls
  under the catch-all warning below; a follow-up issue tracks real support.
- Every OTHER generalized-likelihood distribution (pois/binom/beta/chisq/
  dexp/f/geom/hyper/unif/weibull/gamma/ordinal/general `ll()`) gets a
  concise runInfo warning that censoring is ignored, instead of silently
  wrong results.

## Phases

1. [ ] `src/censEst.h`: `doCensT1()` (M2/M3/M4 value correction via
   `R::pt()`, nu=1 for cauchy) + `dCensT1()` (first derivative). No 2nd/3rd
   order partials -- censored T/Cauchy falls through to the existing
   Gauss-Newton curvature fallback (document as a KNOWN GAP).
2. [ ] R-side plumbing: emit `rx_pred_f_`/`rx_r_`/`rx_nu_` as consecutive
   lhs columns for t/cauchy in both the nlm predOnly model (extend
   `.nlmFixCensRLine`) and the FOCEi inner model (new sibling injection).
3. [ ] `src/nlm.cpp`: widen the three `dist ==` guards to include
   T/Cauchy; read nu from the new slot; call `doCensT1`.
4. [ ] `src/inner.cpp`: widen the ~6 `dist == rxDistributionNorm` sites to
   add a parallel T/Cauchy branch calling `doCensT1`/`dCensT1` for censored
   observations only (uncensored stays bit-identical). Decline
   `fast=TRUE` analytic outer gradient to FD for a subject with a censored
   T/Cauchy observation.
5. [ ] Catch-all warning: new `R/preProcessCensDistWarn.R` hook, runs for
   every estimation method, warns (<75 chars, no method prefix) when
   CENS/LIMIT data meets a non-{norm,dnorm,t,cauchy} endpoint. Skip
   `est="nls"` (already hard stops).
6. [ ] Tests: doCensT1/dCensT1 correctness (numDeriv, nu=1==cauchy
   identity, large-nu==normal limit); FOCEI + nlm regression using the
   issue's repro model (t() and cauchy()), asserting the mechanism
   actually ran; catch-all warning test.
7. [ ] NEWS.md bullet under `## Bug fixes`, referencing #979.
8. [ ] Antigravity review until clean.
9. [ ] Push, open PR closing #979.
