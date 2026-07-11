# Control options for the fsaem (fast SAEM) estimation method

\`fsaem\` is the fast-SAEM (f-SAEM) of Karimi, Lavielle and Moulines
(2020). It is sugar over \`saem\` with \`saemControl(fast=TRUE)\`
forced: the MCMC simulation step samples the individual random effects
from an independent Metropolis-Hastings proposal centered at each
subject's conditional MAP estimate, accelerating the early SAEM
iterations. All other options are the \`saemControl()\` options; see
there (in particular \`fastKernel\`, \`fastCov\`, \`fastIter\` and
\`fastLik\`) for the fast-specific tuning knobs.

## Usage

``` r
fsaemControl(..., fast = TRUE)
```

## Arguments

- ...:

  Parameters used in the default \`saemControl()\`

- fast:

  Always \`TRUE\` for \`fsaemControl()\` and cannot be changed – use
  \`saemControl(fast=FALSE)\` (or \`est="saem"\`) for standard SAEM.

## Value

fsaemControl object (a \`saemControl\` with \`fast=TRUE\`)

## Author

Matthew L. Fidler

## Examples

``` r

fsaemControl()
#> $mcmc
#> $mcmc$niter
#> [1] 200 300
#> 
#> $mcmc$nmc
#> [1] 3
#> 
#> $mcmc$nu
#> [1] 2 2 2
#> 
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
#> [1] 1e-08
#> 
#> $rtol
#> [1] 1e-06
#> 
#> $maxsteps
#> [1] 70000
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
#> [1] 1e-08
#> 
#> $ssRtol
#> [1] 1e-06
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
#> [1] 1e-08
#> 
#> $rtolSens
#> [1] 1e-06
#> 
#> $ssAtolSens
#> [1] 1e-08
#> 
#> $ssRtolSens
#> [1] 1e-06
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
#> $seed
#> [1] 99
#> 
#> $censOption
#> [1] 0
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
#> $DEBUG
#> [1] 0
#> 
#> $optExpression
#> [1] TRUE
#> 
#> $literalFix
#> [1] FALSE
#> 
#> $sumProd
#> [1] FALSE
#> 
#> $nnodesGq
#> [1] 3
#> 
#> $nsdGq
#> [1] 1.6
#> 
#> $adjObf
#> [1] TRUE
#> 
#> $addProp
#> [1] "combined2"
#> 
#> $itmax
#> [1] 30
#> 
#> $tol
#> [1] 1e-06
#> 
#> $type
#> [1] "newuoa"
#> 
#> $powRange
#> [1] 10
#> 
#> $lambdaRange
#> [1] 3
#> 
#> $odeRecalcFactor
#> [1] 3.162278
#> 
#> $maxOdeRecalc
#> [1] 5
#> 
#> $indTolRelax
#> [1] TRUE
#> 
#> $perSa
#> [1] 0.75
#> 
#> $perNoCor
#> [1] 0.75
#> 
#> $perFixOmega
#> [1] 0.1
#> 
#> $perFixResid
#> [1] 0.1
#> 
#> $compress
#> [1] TRUE
#> 
#> $genRxControl
#> [1] TRUE
#> 
#> $sigdigTable
#> [1] 3
#> 
#> $ci
#> [1] 0.95
#> 
#> $covMethod
#> [1] "linFim"
#> 
#> $covFull
#> [1] TRUE
#> 
#> $nSaCov
#> [1] 500
#> 
#> $logLik
#> [1] FALSE
#> 
#> $calcTables
#> [1] TRUE
#> 
#> $muRefCov
#> [1] TRUE
#> 
#> $muRefCovAlg
#> [1] TRUE
#> 
#> $handleUninformativeEtas
#> [1] TRUE
#> 
#> $iovXform
#> [1] "sd"
#> 
#> $boundedTransform
#> [1] TRUE
#> 
#> $eventSens
#> [1] "jump"
#> 
#> $mixProbMethod
#> [1] "regularized"
#> 
#> $mixProbStepExp
#> [1] 1
#> 
#> $mixProbPriorN
#> [1] 20
#> 
#> $mixSampleMethod
#> [1] "parallel"
#> 
#> $fast
#> [1] TRUE
#> 
#> $fastKernel
#> [1] "firstN"
#> 
#> $fastCov
#> [1] "auto"
#> 
#> $fastIter
#> [1] 20
#> 
#> $fastLik
#> [1] "focei"
#> 
#> $lbfgsLmm
#> [1] 5
#> 
#> $lbfgsFactr
#> [1] 1e+07
#> 
#> $lbfgsPgtol
#> [1] 0
#> 
#> $lbfgsMaxIter
#> [1] 20
#> 
#> $nRetry
#> [1] 10
#> 
#> attr(,"class")
#> [1] "fsaemControl"
```
