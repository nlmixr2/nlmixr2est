# Control options for the impmap (importance-sampling EM) estimation method

A NONMEM-style Monte Carlo importance-sampling EM built on the
mu-referenced FOCEI MAP. The proposal density for each subject is
centered at the MAP mode (\`muModel="lin"\`); mu-referenced population
parameters are updated by the EM gradient, while non-mu parameters
(structural and residual error) are updated by a symbolic-sensitivity
Newton step – the importance-sampling-weighted score and Gauss-Newton
information built from the analytic \`d(f)/d(theta)\` and
\`d(V)/d(theta)\` (exact censored partials for BLQ/M2/M3/M4 points).

## Usage

``` r
impmapControl(
  sigdig = 4,
  ...,
  isample = 300L,
  nIter = 100L,
  mapIter = 1L,
  gamma = 1,
  iscaleMin = 0.1,
  iscaleMax = 10,
  iaccept = 0.4,
  ctol = NULL,
  nConvWindow = 10L,
  impSeed = 42L,
  covMethod = c("imp", "analytic", "r,s", "r", "s", "sa", ""),
  qr = FALSE,
  qrShift = TRUE,
  qrRefresh = TRUE,
  sir = FALSE,
  sirSample = NULL,
  muModel = c("lin", "none")
)
```

## Arguments

- sigdig:

  Optimization significant digits. One value drives, with a single
  consistent formula, the inner/outer optimizer convergence tolerance
  (`10^-sigdig`), the boundary check tolerance (`5*10^(-sigdig+1)`), and
  the ODE solver tolerances: the `rtol` exponent IS `sigdig` and `atol`
  sits three orders below, so `rtol = 10^-sigdig`,
  `atol = 10^(-sigdig-3)` for every solver (stiff, non-stiff or
  auto-switching). The sensitivity (`atolSens`/`rtolSens`) tolerances
  match the main solve (the outer gradient and covariance are built from
  them); the steady-state (`ssAtol`/`ssRtol`) tolerances run one order
  looser. Keying the optimizer to the same `10^-sigdig` means it
  converges to exactly the precision the solve supports. At the default
  `sigdig = 4` this is `atol = 1e-7`, `rtol = 1e-4`.

- ...:

  Parameters used in the default \`foceiControl()\`

- isample:

  Number of importance samples drawn per subject per iteration (NONMEM
  ISAMPLE).

- nIter:

  Maximum number of importance-sampling EM iterations.

- mapIter:

  Number of MAP re-centering iterations per EM step; \`\> 0\` re-centers
  the proposal at the MAP mode each iteration.

- gamma:

  Initial proposal-variance inflation factor (NONMEM ISCALE); the
  proposal covariance is \`gamma\` times the inverse of the inner
  information matrix at the mode.

- iscaleMin, iscaleMax:

  Lower/upper bounds for the adapted \`gamma\` (NONMEM ISCALE_MIN /
  ISCALE_MAX).

- iaccept:

  Minimum importance-sampling effective-sample fraction (NONMEM
  IACCEPT). The proposal scale \`gamma\` is kept at its efficient
  starting value while the achieved fraction stays at or above
  \`iaccept\`, and is inflated (toward \`iscaleMax\`) only when it drops
  below this floor.

- ctol:

  Convergence tolerance on the windowed objective-function change;
  \`NULL\` derives it from \`sigdig\`.

- nConvWindow:

  Length of the trailing iteration window used to average the
  objective-function change for convergence (NONMEM-style CTYPE).

- impSeed:

  Base seed for the per-subject thread-safe (threefry) RNG streams;
  results are reproducible and independent of the thread count.

- covMethod:

  Covariance method. \`"imp"\` (default) computes the Monte-Carlo
  importance-sampling observed-information covariance for the estimated
  thetas and Omega parameters (a finite-difference Hessian of the
  importance-sampling objective over fixed common-random-number
  samples), stashed as \`\$impCov\` / \`\$impSe\` and installed as the
  fit covariance; the theta standard errors match the Hessian-based
  FOCEI covariance, though the variance of a tightly-determined random
  effect (an Omega diagonal) can be over-estimated because the fixed
  samples barely span its prior variation. \`"analytic"\`, \`"r,s"\`,
  \`"r"\`, \`"s"\` instead compute the FOCEI covariance post-fit at the
  converged estimates (see \[foceiControl()\]); \`""\` skips the
  covariance step.

- qr:

  When \`TRUE\`, draw quasi-random (Sobol low-discrepancy) importance
  samples instead of pseudo-random Gaussian samples (QRPEM, Leary &
  Dunlavey PAGE 2012); the E-step integrals converge at O(1/N) instead
  of O(1/sqrt(N)).

- qrShift:

  Only used with \`qr=TRUE\`. When \`TRUE\` each (iteration, subject)
  applies a random Cranley-Patterson shift to the Sobol points (seeded,
  thread-count independent); \`FALSE\` reuses one fixed Sobol point set
  everywhere (fully deterministic E-step, no RNG in the draw).

- qrRefresh:

  Only used with \`qr=TRUE\` and \`qrShift=TRUE\`. When \`TRUE\` the
  shift is redrawn each iteration so residual quasi-random error
  averages out over the EM; \`FALSE\` draws one shift per subject at the
  fit start, making each EM iteration a deterministic map (smoothest
  objective trace).

- sir:

  When \`TRUE\`, accelerate the non-mu / residual-error M-step by SIR
  (sampling-importance-resampling): the theta-sensitivity Newton step
  uses \`sirSample\` equal-weight resampled points per subject instead
  of all \`isample\` weighted samples.

- sirSample:

  Number of SIR resampled points per subject; \`NULL\` uses \`max(25,
  ceiling(isample/10))\`. Must be at most \`isample\`.

- muModel:

  Mu-referencing variant for the MAP inner problem; for
  \`impmapControl()\` this is always \`"lin"\` and cannot be changed.

## Value

impmapControl object

## Author

Matthew L. Fidler

## Examples

``` r

impmapControl()
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
#> $covMethodDeferred
#> [1] NA
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
#> [1] 4
#> 
#> $sigdigTable
#> [1] 4
#> 
#> $scaleObjective
#> [1] 0
#> 
#> $boundTol
#> [1] 0.005
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
#>     .ctl <- .controlMaxfun(control)
#>     if (is.null(.ctl$npt)) 
#>         .ctl$npt <- length(par) * 2 + 1
#>     .ctl$iprint <- 0L
#>     .ctl <- .ctl[names(.ctl) %in% c("npt", "rhobeg", "rhoend", 
#>         "iprint", "maxfun")]
#>     .ret <- minqa::bobyqa(par, fn, control = .ctl, lower = lower, 
#>         upper = upper)
#>     .ret$x <- .ret$par
#>     .ret$message <- .ret$msg
#>     .ret$convergence <- .ret$ierr
#>     .ret$value <- .ret$fval
#>     .ret
#> }
#> <bytecode: 0x55bfc6a7f330>
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
#> $scaleCband
#> [1]  0.1 10.0
#> 
#> $scaleCmax
#> [1] 1e+05
#> 
#> $scaleC0
#> [1] 1e+05
#> 
#> $outerOptTxt
#> [1] "bobyqa"
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
#> [1] 1e-07
#> 
#> $rtol
#> [1] 1e-04
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
#> [1] 1e-06
#> 
#> $ssRtol
#> [1] 0.001
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
#> [1] 1e-07
#> 
#> $rtolSens
#> [1] 1e-04
#> 
#> $ssAtolSens
#> [1] 1e-06
#> 
#> $ssRtolSens
#> [1] 0.001
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
#> $shi21hMax
#> [1] 2
#> 
#> $shi21hMin
#> [1] 1e-04
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
#> $zeroTheta
#> [1] 0.001
#> 
#> $impCov
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
#> attr(,"class")
#> [1] "impmapControl"
```
