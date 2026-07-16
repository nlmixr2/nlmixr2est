# Control options for the magq estimation method

Mu-referenced-FOCEI-family closed-form-regression (\`"lin"\`) variant of
adaptive Gauss-Hermite quadrature; see \`foceiControl(muModel=)\`.

## Usage

``` r
magqControl(
  sigdig = 3,
  nAGQ = 2,
  ...,
  interaction = TRUE,
  agqLow = -Inf,
  agqHi = Inf,
  muModel = c("lin", "irls", "none")
)
```

## Arguments

- sigdig:

  Optimization significant digits; controls the inner/outer optimization
  tolerance (`10^-sigdig`), ODE solver tolerance (`0.5*10^(-sigdig-2)`,
  or `0.5*10^(-sigdig-1.5)` for sensitivity/steady-state with liblsoda),
  and boundary check tolerance (`5*10^(-sigdig+1)`).

- nAGQ:

  Number of Gauss-Hermite adaptive quadrature points. \`0\` disables
  AGQ; \`1\` is equivalent to Laplace. Cost grows quickly with ETAs:
  once the EBE is found, expect \`nAGQ^neta\` (even \`nAGQ\`) or
  \`(nAGQ^neta)-1\` (odd \`nAGQ\`) additional evaluations per subject.

- ...:

  Parameters used in the default \`foceiControl()\`

- interaction:

  boolean, Interaction term for the model, in this case the default is
  \`TRUE\`; For adaptive quadrature, with normal distribution the
  Hessian is calculated with the foce(i) approximation

- agqLow:

  The lower bound for adaptive quadrature log-likelihood. By default
  this is -Inf; in the original nlmixr's gnlmm it was -700.

- agqHi:

  The upper bound for adaptive quadrature log-likelihood. By default
  this is Inf; in the original nlmixr's gnlmm was 400.

- muModel:

  Selects the regression variant; for \`magqControl()\` this is always
  \`"lin"\` and cannot be changed – use \`iagqControl()\` for the IRLS
  variant.

## Value

magqControl object

## Author

Matthew L. Fidler

## Examples

``` r

magqControl()
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
#> <bytecode: 0x5615a07d3fd8>
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
#> [1] 2
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
#> attr(,"class")
#> [1] "magqControl"
```
