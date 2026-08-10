# Control for the fbvi (full-Bayes variational inference) method

A convenience wrapper around \[emviControl()\] with \`pointEstimate =
FALSE\`, i.e. the variational posterior covers the unconstrained
population vector (under flat priors) as well as the per-subject etas,
rather than point-estimating the population parameters by an M-step. See
\[emviControl()\] for the full parameter list.

## Usage

``` r
fbviControl(..., pointEstimate = FALSE)
```

## Arguments

- ...:

  Parameters passed to \[emviControl()\].

- pointEstimate:

  Which of the two methods to run, normally left at its \`NULL\` default
  so it follows \`est\`: \`est="emvi"\` implies \`TRUE\` and
  \`est="fbvi"\` implies \`FALSE\`. \`est\` wins over a contradicting
  value, and says so – it has to, because re-estimating a fit with the
  other method pipes the completed fit's control forward. \`TRUE\` runs
  the variational-EM hybrid: the variational posterior covers the
  per-subject etas only, and the population parameters (thetas / omega /
  residual error) are point estimates maximized by the ELBO gradient;
  output semantics match FOCEi/SAEM. \`FALSE\` runs full Bayes: the
  variational posterior also covers the unconstrained population vector,
  with flat priors.

  Two things about "flat" are worth being explicit about, because they
  define the prior rather than merely describe the implementation. (1) A
  BOUNDED theta is fitted on its unconstrained scale, and the
  log-determinant of that constraining transform IS added to the
  full-Bayes objective, so the flat prior is flat on the NATURAL
  parameter, as in Stan. It is deliberately NOT added when
  \`pointEstimate=TRUE\`: a maximum-likelihood estimate has to stay
  invariant to reparameterization, which is why Stan's own \`optimize\`
  defaults to \`jacobian=0\`. (2) The between-subject variances are
  carried as per-eta LOG-variances and no Jacobian is applied to them,
  so the prior is flat on \`log(omega)\` – the conventional
  weakly-informative choice for a scale parameter, but a choice, not an
  accident: it is not flat on \`omega\`.

  These point estimates maximize the ELBO, NOT the likelihood, and the
  difference is not merely cosmetic. Since \`ELBO = log p(y\|theta) -
  KL(q \|\| p(eta\|y,theta))\`, any dependence of that KL on \`theta\`
  displaces the maximizer from the MLE – "variational maximum
  likelihood". For VARIANCE components the displacement has a known
  direction: a variational family that understates posterior spread
  makes the omega M-step, \`Omega = mean_i(mu_i mu_i' + Sigma_i)\`,
  inherit that understatement, so between-subject variability is biased
  DOWNWARD. The bias is worst for \`viFamily="meanField"\`, which cannot
  represent within-subject posterior correlation at all; \`"fullRank"\`
  can, which is why it is the default. Structural (typical-value)
  parameters are far less affected. If the between-subject variances are
  themselves the quantity of interest, prefer \`"fullRank"\` and
  cross-check against \`est="focei"\` or \`est="saem"\`.

## Value

An \`emviControl\` object with \`pointEstimate = FALSE\`.

## Author

Matthew L. Fidler

## Examples

``` r

fbviControl()
#> $seed
#> [1] 42
#> 
#> $iters
#> [1] 300
#> 
#> $nMc
#> [1] 1
#> 
#> $viFamily
#> [1] "fullRank"
#> 
#> $pointEstimate
#> [1] FALSE
#> 
#> $optim
#> [1] "advi"
#> 
#> $adaptEta
#> [1] TRUE
#> 
#> $perNoCor
#> [1] 0.75
#> 
#> $etaCandidates
#> [1] 0.010 0.025 0.050 0.100 0.250
#> 
#> $tau
#> [1] 1
#> 
#> $alpha
#> [1] 0.1
#> 
#> $tol
#> [1] 0.001
#> 
#> $evalElbo
#> [1] 100
#> 
#> $klWarmup
#> [1] 0
#> 
#> $temperInit
#> [1] 10
#> 
#> $likelihood
#> [1] "focei"
#> 
#> $returnVi
#> [1] FALSE
#> 
#> $resume
#> NULL
#> 
#> $covMethod
#> [1] "vi"
#> 
#> $optExpression
#> [1] TRUE
#> 
#> $sumProd
#> [1] FALSE
#> 
#> $literalFix
#> [1] TRUE
#> 
#> $literalFixRes
#> [1] TRUE
#> 
#> $addProp
#> [1] "combined2"
#> 
#> $calcTables
#> [1] TRUE
#> 
#> $compress
#> [1] FALSE
#> 
#> $adjObf
#> [1] TRUE
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
#> $stickyRecalcN
#> [1] 4
#> 
#> $maxOdeRecalc
#> [1] 5
#> 
#> $odeRecalcFactor
#> [1] 3.162278
#> 
#> $indTolRelax
#> [1] TRUE
#> 
#> $eventSens
#> [1] "jump"
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
#> $rxControl
#> $scale
#> NULL
#> 
#> $method
#> liblsoda 
#>        2 
#> 
#> $atol
#> [1] 1e-06
#> 
#> $rtol
#> [1] 0.001
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
#> Al-Mohy 
#>       3 
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
#> [1] 1e-05
#> 
#> $ssRtol
#> [1] 0.01
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
#> [1] 1e-06
#> 
#> $rtolSens
#> [1] 0.001
#> 
#> $ssAtolSens
#> [1] 1e-05
#> 
#> $ssRtolSens
#> [1] 0.01
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
#> $zeroVarParamHandle
#> [1] "warn"
#> 
#> $indLinStepSearch
#> [1] 1
#> 
#> $indLinMaxIter
#> [1] 20
#> 
#> $indLinRichardson
#> [1] 2
#> 
#> $indLinIteration
#> [1] 3
#> 
#> $indLinJac
#> [1] 0
#> 
#> attr(,"class")
#> [1] "rxControl"
#> 
#> $genRxControl
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "emviControl"
```
