# Control options for the foce estimation method

This is the first order option without the interaction between residuals
and etas.

## Usage

``` r
foceControl(sigdig = 3, ..., interaction = FALSE)
```

## Arguments

- sigdig:

  Optimization significant digits. This controls:

  - The tolerance of the inner and outer optimization is `10^-sigdig`

  - The tolerance of the ODE solvers is `0.5*10^(-sigdig-2)`; For the
    sensitivity equations and steady-state solutions the default is
    `0.5*10^(-sigdig-1.5)` (sensitivity changes only applicable for
    liblsoda)

  - The tolerance of the boundary check is `5 * 10 ^ (-sigdig + 1)`

- ...:

  Parameters used in the default \`foceiControl()\`

- interaction:

  Interaction term for the model, in this case the default is \`FALSE\`;
  it cannot be changed, use \`focei\` instead

## Value

foceControl object

## Author

Matthew L. Fidler

## Examples

``` r
foceControl()
#> $maxOuterIterations
#> [1] 5000
#> 
#> $maxInnerIterations
#> [1] 1000
#> 
#> $n1qn1nsim
#> [1] 10001
#> 
#> $print
#> [1] 1
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
#> [1] 1
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
#> $useColor
#> [1] TRUE
#> 
#> $boundTol
#> [1] 0.05
#> 
#> $calcTables
#> [1] TRUE
#> 
#> $printNcol
#> [1] 4
#> 
#> $noAbort
#> [1] 1
#> 
#> $interaction
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
#>     .ctl <- control
#>     .ctl <- .ctl[names(.ctl) %in% c("eval.max", "iter.max", "trace", 
#>         "abs.tol", "rel.tol", "x.tol", "xf.tol", "step.min", 
#>         "step.max", "sing.tol", "scale.inti", "diff.g")]
#>     .ctl$trace <- 0
#>     .ret <- stats::nlminb(start = par, objective = fn, gradient = gr, 
#>         hessian = NULL, control = .ctl, lower = lower, upper = upper)
#>     .ret$x <- .ret$par
#>     .ret
#> }
#> <bytecode: 0x56296a8e0dd8>
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
#> $eventType
#> [1] 2
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
#> [1] -1
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
#> attr(,"class")
#> [1] "foceControl"
```
