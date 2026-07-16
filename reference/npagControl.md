# Control for the npag (nonparametric adaptive grid) method

A wrapper around \[impmapControl()\] that reuses the shared FOCEI family
plumbing for the nonparametric adaptive grid engine. The nonparametric
support-point knobs are added in a later milestone.

## Usage

``` r
npagControl(
  points = NULL,
  cycles = 100L,
  gammaOptimize = TRUE,
  residOptimize = c("alternate", "final", "none"),
  muExpand = FALSE,
  gridWidth = 4,
  gridBounds = c("auto", "ini", "both"),
  ...
)
```

## Arguments

- points:

  Initial Sobol grid size (support points). \`NULL\` (default) picks it
  automatically from the number of support-point dimensions (etas):
  \`max(2028, 512 \* n_eta)\` – a fixed grid (Pmetrics uses 2028) covers
  a low-dimensional model but grows sparse and can collapse in high
  dimensions, so the auto size floors at 2028 and scales up per added
  eta. Supply an integer to override.

- cycles:

  Maximum adaptive-grid cycles.

- gammaOptimize:

  Use a global assay-error multiplier (gamma) as a per-cycle warm start
  for the overall residual magnitude, folded into the variance-scale
  coefficients (\`add\`/\`prop\`/\`lnorm\`). The per-endpoint values,
  the add/prop ratio, and the transform/autocorrelation parameters come
  from `residOptimize`. Only valid for normal endpoints; censoring and
  transform-both-sides are supported.

- residOptimize:

  How to estimate the residual-error thetas (every endpoint's
  \`add\`/\`prop\`/\`lnorm\`, each transform \`lambda\`, each \`ar\`)
  with the support points and weights held fixed, using bounded
  [`minqa::bobyqa`](https://rdrr.io/pkg/minqa/man/bobyqa.html) on the
  EXTENDED LEAST SQUARES objective \`sum_obs((f-dv)^2/r + log(r))\` at
  the posterior-mean etas. The \`log(r)\` term keeps the residual from
  drifting to zero on a flexible support (which the marginal likelihood
  would reward), giving the saem/focei residual; each variance scale is
  warm-started from the per-endpoint moment (additive SD from
  \`sqrt(mean(err^2))\`, proportional from \`sqrt(mean((err/f)^2))\`, on
  the transform- both-sides scale). `"alternate"` (default) optimizes
  every cycle (block- coordinate ascent); `"final"` optimizes once at
  the converged support; `"none"` holds them at their initial values.
  Fixed residual parameters are always held. After the residual thetas
  converge, a final adaptive-grid pass re-optimizes the support with
  them held constant so the support remains the nonparametric MLE (D(F)
  ~ 0) for the fitted residual.

- muExpand:

  how to estimate non-mu structural fixed-effect parameters (a theta
  with no eta, e.g. \`ke \<- exp(tke)\`; npag's grid otherwise covers
  only mu-referenced and residual/likelihood parameters). \`FALSE\`
  (default) optimizes them directly as "regressors" in the residual step
  – \`bobyqa\` moves them alongside the residual parameters, re-deriving
  the posterior-mean etas each candidate (so the eta grid cannot
  stale-absorb the structural shift) – which identifies them well (e.g.
  recovering a clearance from a poor start). Not available for mix()
  models (the ELS step is not mixture-aware; component parameters are
  held). \`TRUE\` instead uses the saem-style mu-expansion: inject a
  pseudo-eta (\`ke \<- exp(tke + eta.tke)\`), grid-estimate, and recover
  it as a fixed effect at finalization (support-mean folded into the
  theta, injected random effect collapsed). The regressor default
  usually identifies these parameters more sharply than the grid.

- gridWidth:

  support-point box half-width, in initial-eta SDs, for the
  \`gridBounds="auto"\` grid (default 4). A narrower box focuses the
  initial Sobol grid on the plausible region – useful for
  high-dimensional models where a wide box wastes points on
  near-zero-density support (which can collapse the fit).

- gridBounds:

  how to set the initial support-point box: \`"auto"\` (default) uses
  \`+/- gridWidth \* initial eta SD\`; \`"ini"\` uses each mu-referenced
  parameter's ini-block lower/upper bounds where finite (else auto);
  \`"both"\` uses the ini bounds when present and auto otherwise. For a
  high-dimensional model, bounded ini estimates + \`"ini"\` keep the
  grid in range.

- ...:

  Parameters passed to \[impmapControl()\].

## Value

An \`impmapControl\` object tagged for the npag engine.

## Details

Note: the npag objective is the nonparametric marginal log-likelihood
and uses a different constant convention than NONMEM/FOCEI, so its
\`-2LL\` is NOT comparable to nlmixr2's FOCEI/SAEM/FOCE \`-2LL\`.
Compare npag runs to each other or to Pmetrics NPAG.

Note on residual error with a flexible support distribution: the
residual parameters are estimated against the nonparametric objective
(see `residOptimize`) with the support-point distribution held fixed.
Because that distribution is flexible, it can absorb variability a
parametric model (FOCEI/SAEM) would attribute to residual error –
especially the additive term of a combined additive+proportional model
at low concentrations. As a result the additive coefficient of a
combined error model may be estimated smaller (sometimes toward zero)
than the corresponding parametric fit, while the proportional term and
per-endpoint magnitudes are recovered well. This is an expected property
of nonparametric estimation, not a convergence failure; use
`residOptimize = "none"` to hold the residual parameters at their
initial values if a fixed error model is desired.

## Author

Matthew L. Fidler

## Examples

``` r

npagControl()
#> $maxOuterIterations
#> [1] 5000
#> 
#> $maxInnerIterations
#> [1] 1000
#> 
#> $n1qn1nsim
#> [1] 10001
#> 
#> $iterPrintControl
#> $every
#> [1] 1
#> 
#> $ncol
#> [1] 4
#> 
#> $headerEvery
#> [1] 10
#> 
#> $useColor
#> [1] TRUE
#> 
#> $simple
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "iterPrintControl" "list"            
#> 
#> $lbfgsLmm
#> [1] 7
#> 
#> $lbfgsPgtol
#> [1] 0
#> 
#> $lbfgsFactr
#> [1] 4.5036e+11
#> 
#> $scaleTo
#> [1] 1
#> 
#> $epsilon
#> [1] 1e-04
#> 
#> $derivEps
#> [1] 2.980232e-07 2.980232e-07
#> 
#> $derivMethod
#> [1] 3
#> 
#> $covDerivMethod
#> [1] 1
#> 
#> $covMethod
#> [1] 2
#> 
#> $covType
#> [1] "analytic"
#> 
#> $covSolveTol
#> NULL
#> 
#> $covFull
#> [1] TRUE
#> 
#> $fast
#> [1] FALSE
#> 
#> $centralDerivEps
#> [1] 2.980232e-07 2.980232e-07
#> 
#> $eigen
#> [1] 1
#> 
#> $diagXform
#> [1] "sqrt"
#> 
#> $iovXform
#> [1] "sd"
#> 
#> $sumProd
#> [1] FALSE
#> 
#> $optExpression
#> [1] TRUE
#> 
#> $literalFix
#> [1] TRUE
#> 
#> $literalFixRes
#> [1] TRUE
#> 
#> $outerOpt
#> [1] -1
#> 
#> $ci
#> [1] 0.95
#> 
#> $sigdig
#> [1] 3
#> 
#> $sigdigTable
#> [1] 3
#> 
#> $scaleObjective
#> [1] 0
#> 
#> $boundTol
#> [1] 0.05
#> 
#> $calcTables
#> [1] TRUE
#> 
#> $noAbort
#> [1] 1
#> 
#> $interaction
#> [1] 1
#> 
#> $foce
#> [1] "nonmem"
#> 
#> $foceType
#> [1] 0
#> 
#> $cholSEtol
#> [1] 6.055454e-06
#> 
#> $hessEps
#> [1] 6.055454e-06
#> 
#> $hessEpsLlik
#> [1] 6.055454e-06
#> 
#> $optimHessType
#> [1] 1
#> 
#> $optimHessCovType
#> [1] 1
#> 
#> $censOption
#> [1] 0
#> 
#> $cholAccept
#> [1] 0.001
#> 
#> $resetEtaSize
#> [1] 1.439531
#> 
#> $resetThetaSize
#> [1] 1.959964
#> 
#> $resetThetaFinalSize
#> [1] 1.439531
#> 
#> $diagOmegaBoundUpper
#> [1] 5
#> 
#> $diagOmegaBoundLower
#> [1] 100
#> 
#> $cholSEOpt
#> [1] 0
#> 
#> $cholSECov
#> [1] 0
#> 
#> $fo
#> [1] 0
#> 
#> $covTryHarder
#> [1] 0
#> 
#> $outerOptFun
#> function (par, fn, gr, lower = -Inf, upper = Inf, control = list(), 
#>     ...) 
#> {
#>     .ctl <- .controlIterMax(control)
#>     .ctl <- .ctl[names(.ctl) %in% c("eval.max", "iter.max", "trace", 
#>         "abs.tol", "rel.tol", "x.tol", "xf.tol", "step.min", 
#>         "step.max", "sing.tol", "scale.inti", "diff.g")]
#>     .ctl$trace <- 0
#>     .ret <- stats::nlminb(start = par, objective = fn, gradient = gr, 
#>         hessian = NULL, control = .ctl, lower = lower, upper = upper)
#>     .ret$x <- .ret$par
#>     .ret
#> }
#> <bytecode: 0x5572b3270188>
#> <environment: namespace:nlmixr2est>
#> 
#> $rhobeg
#> [1] 0.2
#> 
#> $rhoend
#> [1] 1e-04
#> 
#> $npt
#> NULL
#> 
#> $rel.tol
#> [1] 1e-04
#> 
#> $x.tol
#> [1] 1e-04
#> 
#> $eval.max
#> [1] 4000
#> 
#> $iter.max
#> [1] 2000
#> 
#> $innerOpt
#> [1] 1
#> 
#> $abstol
#> [1] 1e-04
#> 
#> $reltol
#> [1] 1e-04
#> 
#> $derivSwitchTol
#> [1] 2e-04
#> 
#> $resetHessianAndEta
#> [1] 0
#> 
#> $muModel
#> [1] "lin"
#> 
#> $muRefCovAlg
#> [1] TRUE
#> 
#> $muModelTol
#> [1] 1e-05
#> 
#> $muModelMaxCycles
#> [1] 20
#> 
#> $muModelClampRetries
#> [1] 10
#> 
#> $stateTrim
#> [1] Inf
#> 
#> $gillK
#> [1] 10
#> 
#> $gillKcov
#> [1] 10
#> 
#> $gillKcovLlik
#> [1] 10
#> 
#> $gillRtol
#> [1] 1.490116e-08
#> 
#> $gillStep
#> [1] 4
#> 
#> $gillStepCov
#> [1] 2
#> 
#> $gillStepCovLlik
#> [1] 4.5
#> 
#> $scaleType
#> [1] 2
#> 
#> $normType
#> [1] 1
#> 
#> $scaleC
#> NULL
#> 
#> $scaleCmin
#> [1] 1e-05
#> 
#> $scaleCmax
#> [1] 1e+05
#> 
#> $scaleC0
#> [1] 1e+05
#> 
#> $outerOptTxt
#> [1] "nlminb"
#> 
#> $outerOptDefault
#> [1] TRUE
#> 
#> $rmatNorm
#> [1] 1
#> 
#> $rmatNormLlik
#> [1] 1
#> 
#> $smatNorm
#> [1] 1
#> 
#> $smatNormLlik
#> [1] 1
#> 
#> $covGillF
#> [1] 1
#> 
#> $optGillF
#> [1] 1
#> 
#> $gillFtol
#> [1] 0
#> 
#> $gillFtolCov
#> [1] 0
#> 
#> $gillFtolCovLlik
#> [1] 0
#> 
#> $covSmall
#> [1] 1e-05
#> 
#> $adjLik
#> [1] TRUE
#> 
#> $gradTrim
#> [1] Inf
#> 
#> $gradCalcCentralSmall
#> [1] 1e-04
#> 
#> $gradCalcCentralLarge
#> [1] 10000
#> 
#> $etaNudge
#> [1] 1.131586
#> 
#> $etaNudge2
#> [1] 1.518182
#> 
#> $maxOdeRecalc
#> [1] 5
#> 
#> $odeRecalcFactor
#> [1] 3.162278
#> 
#> $nRetries
#> [1] 3
#> 
#> $seed
#> [1] 42
#> 
#> $resetThetaCheckPer
#> [1] 0.1
#> 
#> $etaMat
#> NULL
#> 
#> $repeatGillMax
#> [1] 1
#> 
#> $stickyRecalcN
#> [1] 4
#> 
#> $indTolRelax
#> [1] TRUE
#> 
#> $eventType
#> [1] 2
#> 
#> $eventSens
#> [1] "jump"
#> 
#> $gradProgressOfvTime
#> [1] 10
#> 
#> $addProp
#> [1] "combined2"
#> 
#> $badSolveObjfAdj
#> [1] 100
#> 
#> $compress
#> [1] FALSE
#> 
#> $rxControl
#> $scale
#> NULL
#> 
#> $method
#> liblsoda 
#>        2 
#> 
#> $atol
#> [1] 5e-06
#> 
#> $rtol
#> [1] 5e-06
#> 
#> $maxsteps
#> [1] 500000
#> 
#> $hmin
#> [1] 0
#> 
#> $hmax
#> [1] NA
#> 
#> $hini
#> [1] 0
#> 
#> $maxordn
#> [1] 12
#> 
#> $maxords
#> [1] 5
#> 
#> $covsInterpolation
#> locf 
#>    1 
#> 
#> $addCov
#> [1] TRUE
#> 
#> $returnType
#> rxSolve 
#>       0 
#> 
#> $sigma
#> NULL
#> 
#> $sigmaDf
#> NULL
#> 
#> $nCoresRV
#> [1] 1
#> 
#> $sigmaIsChol
#> [1] FALSE
#> 
#> $sigmaSeparation
#> [1] "auto"
#> 
#> $sigmaXform
#> identity 
#>        4 
#> 
#> $nDisplayProgress
#> [1] 10000
#> 
#> $amountUnits
#> [1] NA
#> 
#> $timeUnits
#> [1] "hours"
#> 
#> $addDosing
#> [1] FALSE
#> 
#> $stateTrim
#> [1] Inf
#> 
#> $updateObject
#> [1] FALSE
#> 
#> $omega
#> NULL
#> 
#> $omegaDf
#> NULL
#> 
#> $omegaIsChol
#> [1] FALSE
#> 
#> $omegaSeparation
#> [1] "auto"
#> 
#> $omegaXform
#> variance 
#>        6 
#> 
#> $nSub
#> [1] 1
#> 
#> $thetaMat
#> NULL
#> 
#> $thetaDf
#> NULL
#> 
#> $thetaIsChol
#> [1] FALSE
#> 
#> $nStud
#> [1] 1
#> 
#> $dfSub
#> [1] 0
#> 
#> $dfObs
#> [1] 0
#> 
#> $seed
#> NULL
#> 
#> $nsim
#> NULL
#> 
#> $minSS
#> [1] 10
#> 
#> $maxSS
#> [1] 10000
#> 
#> $strictSS
#> [1] 1
#> 
#> $infSSstep
#> [1] 12
#> 
#> $istateReset
#> [1] TRUE
#> 
#> $subsetNonmem
#> [1] TRUE
#> 
#> $hmaxSd
#> [1] 0
#> 
#> $maxAtolRtolFactor
#> [1] 0.1
#> 
#> $from
#> NULL
#> 
#> $to
#> NULL
#> 
#> $by
#> NULL
#> 
#> $length.out
#> NULL
#> 
#> $iCov
#> NULL
#> 
#> $keep
#> NULL
#> 
#> $keepF
#> character(0)
#> 
#> $drop
#> NULL
#> 
#> $warnDrop
#> [1] TRUE
#> 
#> $omegaLower
#> [1] -Inf
#> 
#> $omegaUpper
#> [1] Inf
#> 
#> $sigmaLower
#> [1] -Inf
#> 
#> $sigmaUpper
#> [1] Inf
#> 
#> $thetaLower
#> [1] -Inf
#> 
#> $thetaUpper
#> [1] Inf
#> 
#> $indLinPhiM
#> [1] 0
#> 
#> $indLinPhiTol
#> [1] 1e-07
#> 
#> $indLinMatExpType
#> expokit 
#>       2 
#> 
#> $indLinMatExpOrder
#> [1] 6
#> 
#> $idFactor
#> [1] TRUE
#> 
#> $mxhnil
#> [1] 0
#> 
#> $hmxi
#> [1] 0
#> 
#> $warnIdSort
#> [1] TRUE
#> 
#> $ssAtol
#> [1] 5e-04
#> 
#> $ssRtol
#> [1] 5e-04
#> 
#> $safeZero
#> [1] 1
#> 
#> $sumType
#> pairwise 
#>        1 
#> 
#> $prodType
#> long double 
#>           1 
#> 
#> $resample
#> NULL
#> 
#> $resampleID
#> [1] TRUE
#> 
#> $maxwhile
#> [1] 100000
#> 
#> $cores
#> [1] 0
#> 
#> $atolSens
#> [1] 1.581139e-05
#> 
#> $rtolSens
#> [1] 1.581139e-05
#> 
#> $ssAtolSens
#> [1] 0.002108483
#> 
#> $ssRtolSens
#> [1] 0.002108483
#> 
#> $simVariability
#> [1] NA
#> 
#> $nLlikAlloc
#> NULL
#> 
#> $useStdPow
#> [1] 0
#> 
#> $naTimeHandle
#> ignore 
#>      1 
#> 
#> $addlKeepsCov
#> [1] FALSE
#> 
#> $addlDropSs
#> [1] TRUE
#> 
#> $ssAtDoseTime
#> [1] TRUE
#> 
#> $ss2cancelAllPending
#> [1] FALSE
#> 
#> $naInterpolation
#> locf 
#>    1 
#> 
#> $keepInterpolation
#> na 
#>  2 
#> 
#> $safeLog
#> [1] 1
#> 
#> $safePow
#> [1] 1
#> 
#> $ssSolved
#> [1] TRUE
#> 
#> $linCmtSensType
#> auto 
#>  100 
#> 
#> $linCmtSensH
#> [1] 1e-04
#> 
#> $linCmtGillFtol
#> [1] 0
#> 
#> $linCmtGillK
#> [1] 20
#> 
#> $linCmtGillStep
#> [1] 4
#> 
#> $linCmtGillRtol
#> [1] 1.490116e-08
#> 
#> $linCmtShiErr
#> [1] 1.490116e-08
#> 
#> $linCmtShiMax
#> [1] 20
#> 
#> $linCmtScale
#> [1] 0 0 0 0 0 0 0
#> 
#> $linCmtHcmt
#> [1] 1
#> 
#> $linCmtHmeanI
#> geometric 
#>         2 
#> 
#> $linCmtHmeanO
#> geometric 
#>         2 
#> 
#> $linCmtSuspect
#> [1] 1e-06
#> 
#> $linCmtForwardMax
#> [1] 2
#> 
#> $indOwnAlloc
#> [1] -1
#> 
#> $maxExtra
#> [1] 1000
#> 
#> $tolFactor
#> NULL
#> 
#> $serializeFile
#> NULL
#> 
#> $dense
#> [1] FALSE
#> 
#> $cvodeLinSolver
#> dense 
#>     1 
#> 
#> $stiff2
#> [1] 0
#> 
#> $autoSwitchMaxStiff
#> [1] 10
#> 
#> $autoSwitchMaxNonstiff
#> [1] 3
#> 
#> $autoSwitchStiffFirst
#> [1] 0
#> 
#> $autoSwitchNonstifftol
#> [1] 0.9
#> 
#> $autoSwitchStifftol
#> [1] 0.9
#> 
#> $autoSwitchDtfac
#> [1] 2
#> 
#> $autoSwitchSwitchMax
#> [1] 5
#> 
#> $useLinCmt
#> [1] TRUE
#> 
#> $file
#> NULL
#> 
#> $chunkSize
#> NULL
#> 
#> $parallel
#> [1] 0
#> 
#> $.zeros
#> NULL
#> 
#> attr(,"class")
#> [1] "rxControl"
#> 
#> $genRxControl
#> [1] TRUE
#> 
#> $skipCov
#> NULL
#> 
#> $fallbackFD
#> [1] FALSE
#> 
#> $shi21maxOuter
#> [1] 0
#> 
#> $shi21maxInner
#> [1] 20
#> 
#> $shi21maxInnerCov
#> [1] 20
#> 
#> $shi21maxFD
#> [1] 20
#> 
#> $smatPer
#> [1] 0.6
#> 
#> $sdLowerFact
#> [1] 0.001
#> 
#> $zeroGradFirstReset
#> [1] TRUE
#> 
#> $zeroGradRunReset
#> [1] TRUE
#> 
#> $zeroGradBobyqa
#> [1] TRUE
#> 
#> $mceta
#> [1] -2
#> 
#> $warm
#> [1] 1
#> 
#> $nAGQ
#> [1] 0
#> 
#> $agqHi
#> [1] Inf
#> 
#> $agqLow
#> [1] -Inf
#> 
#> $sensMethod
#> [1] "default"
#> 
#> $boundedTransform
#> [1] TRUE
#> 
#> $isample
#> [1] 300
#> 
#> $nIter
#> [1] 100
#> 
#> $mapIter
#> [1] 1
#> 
#> $gamma
#> [1] 1
#> 
#> $iscaleMin
#> [1] 0.1
#> 
#> $iscaleMax
#> [1] 10
#> 
#> $iaccept
#> [1] 0.4
#> 
#> $nConvWindow
#> [1] 10
#> 
#> $impSeed
#> [1] 42
#> 
#> $impCov
#> [1] FALSE
#> 
#> $qr
#> [1] FALSE
#> 
#> $qrShift
#> [1] TRUE
#> 
#> $qrRefresh
#> [1] TRUE
#> 
#> $sir
#> [1] FALSE
#> 
#> $sirSample
#> [1] 30
#> 
#> $est
#> [1] "npag"
#> 
#> $points
#> [1] NA
#> 
#> $cycles
#> [1] 100
#> 
#> $gammaOptimize
#> [1] TRUE
#> 
#> $residOptimize
#> [1] "alternate"
#> 
#> $muExpand
#> [1] FALSE
#> 
#> $gridWidth
#> [1] 4
#> 
#> $gridBounds
#> [1] "auto"
#> 
#> attr(,"class")
#> [1] "impmapControl"
```
