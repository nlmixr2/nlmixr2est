# calc.2LL transform-both-sides Jacobian sign (#903)

## Status: measured -- the defect is REAL

## Measurement

Model: 1-compartment IV bolus, ONE eta (on CL), `cp ~ lnorm(prop.sd)`, 12
subjects x 7 observations, simulated.  A single eta makes the marginal
likelihood a 1-D integral, and an IV bolus makes the prediction closed form,
so the reference is `stats::integrate()` over plain R arithmetic -- it shares
no code with saem.

Reference (density of the ORIGINAL DV, i.e. including the `1/DV` Jacobian):

```
-2LL(truth)                 = 113.860799
-2LL(transformed scale)     = -50.255372
sum(log DV)                 =  82.058085   (sum(powerL) = -sum(log DV))
```

`calc.2LL` at the fit's estimates, swept over quadrature settings:

| nnodes | nsd | calc.2LL | calc.2LL + 4*sum(log DV) |
|---|---|---|---|
| 15 | 4 | -214.328753 | 113.903589 |
| 21 | 5 | -214.348208 | 113.884134 |
| 25 | 5 | -214.371728 | 113.860613 |

At 25 nodes / 5 sd the *corrected* value agrees with the reference to 2e-4,
while the shipped value is off by `4*sum(powerL)` = -328.23.  The residual at
coarser settings is Gauss-Legendre error, not a second defect.

Cross-check that the reference itself is right: FOCEi/FOCE/Laplace on the same
model at the same parameters return -2LL = 113.8733 (Laplace error 0.012), and
a no-random-effect FOCEi fit -- where the likelihood is exact -- returns
LL = -173.4529 against a hand-computed -173.5005 with the Jacobian ADDED
(-9.3843 if subtracted).  So the rest of the package adds `powerL`; `calc.2LL`
alone subtracts it.

## Blast radius

Every saem fit with a transformed endpoint.  `saemControl(logLik=)` defaults
to FALSE, but that only DEFERS the calculation: reading `$objf`/`$logLik`/
`$AIC`/`$BIC` off the fit routes through `.nmObjEnsureObjective`, which calls
`setOfv(fit, "gauss3_1.6")` -- i.e. `calc.2LL`.  With `logLik=TRUE` it is
computed eagerly instead (confirmed: the reported row equalled
`calc.2LL(9, 1.6)` exactly).  `add()`/`prop()`/`ll()` endpoints give
`powerL == 0` and are unaffected.

After the fix the same default fit reports -2LL = 113.734 against the
reference 113.861 (gauss3_1.6 quadrature error); before the fix it reported
-214.499.

## Found while measuring, NOT fixed here (separate defects)

- `saem` never propagates a `boxCox()`/`yeoJohnson()` lambda into `transMat`.
  `src/saem.cpp` sets the `lambda` member once from the input list (line 1340)
  and never updates it from `lres`, where the estimated -- or even `fixed()` --
  lambda actually lives (`resMat[, 4]`).  So `transMat[, 1]` stays 1 and
  `calc.2LL` evaluates the likelihood at lambda = 1 regardless of what the fit
  reports.  A `boxCox(lambda = fixed(0.5))` fit on data simulated at
  lambda = 0.5 recovers `tcl = 2.006` against a true 1.386 and
  `omega = 0.0036` against a true 0.09, which is what a fit with the transform
  effectively switched off looks like.
- A `logitNorm(sd, 0, 10)` saem fit dies with "missing value where TRUE/FALSE
  needed" during fit finalization.  Verified to fail identically on `main`.

Both are why the new test covers `lnorm()` and `add()` rather than a wider
sweep of transforms; neither is caused by this change.

## Phases

1. [x] Measure (above).
2. [x] Fix the sign in `R/saem_fit_aux.R`.
3. [x] Accumulate `Q` in the log domain (the issue's blocker: `exp(ltot)`
   overflows for a small-SD lognormal, so the natural twin test returns `Inf`).
4. [x] Fix `.saemCalcLikelihood` reading control `nsd.gq`; `saemControl()`
   stores `nsdGq`, so a user-supplied `nsdGq` was silently ignored.
5. [x] Regression test with the closed-form reference.
6. [x] NEWS.
7. [x] Antigravity review until clean.
8. [x] PR closing #903.
