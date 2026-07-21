# Control Values for nlme Fit with extra options for nlmixr

The values supplied in the function call replace the defaults and a list
with all possible arguments is returned. The returned list is used as
the \`control\` argument to the \`nlme\` function.

## Usage

``` r
nlmixr2NlmeControl(
  maxIter = 100,
  pnlsMaxIter = 100,
  msMaxIter = 100,
  minScale = 0.001,
  tolerance = NULL,
  niterEM = 25,
  pnlsTol = NULL,
  msTol = NULL,
  returnObject = FALSE,
  msVerbose = FALSE,
  msWarnNoConv = TRUE,
  gradHess = TRUE,
  apVar = TRUE,
  .relStep = .Machine$double.eps^(1/3),
  minAbsParApVar = 0.05,
  opt = c("nlminb", "nlm"),
  natural = TRUE,
  sigma = NULL,
  optExpression = TRUE,
  literalFix = TRUE,
  sumProd = FALSE,
  rxControl = NULL,
  method = c("ML", "REML"),
  random = NULL,
  fixed = NULL,
  weights = NULL,
  verbose = TRUE,
  returnNlme = FALSE,
  addProp = c("combined2", "combined1"),
  calcTables = TRUE,
  compress = TRUE,
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  muRefCovAlg = TRUE,
  eventSens = c("jump", "fd"),
  print = NULL,
  covMethod = c("nlme", "analytic", "r,s", "r", "s", "sa", "imp", ""),
  ...
)

nlmeControl(
  maxIter = 100,
  pnlsMaxIter = 100,
  msMaxIter = 100,
  minScale = 0.001,
  tolerance = NULL,
  niterEM = 25,
  pnlsTol = NULL,
  msTol = NULL,
  returnObject = FALSE,
  msVerbose = FALSE,
  msWarnNoConv = TRUE,
  gradHess = TRUE,
  apVar = TRUE,
  .relStep = .Machine$double.eps^(1/3),
  minAbsParApVar = 0.05,
  opt = c("nlminb", "nlm"),
  natural = TRUE,
  sigma = NULL,
  optExpression = TRUE,
  literalFix = TRUE,
  sumProd = FALSE,
  rxControl = NULL,
  method = c("ML", "REML"),
  random = NULL,
  fixed = NULL,
  weights = NULL,
  verbose = TRUE,
  returnNlme = FALSE,
  addProp = c("combined2", "combined1"),
  calcTables = TRUE,
  compress = TRUE,
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  muRefCovAlg = TRUE,
  eventSens = c("jump", "fd"),
  print = NULL,
  covMethod = c("nlme", "analytic", "r,s", "r", "s", "sa", "imp", ""),
  ...
)
```

## Arguments

- maxIter:

  maximum number of iterations for the `nlme` optimization algorithm.
  Default is 50.

- pnlsMaxIter:

  maximum number of iterations for the `PNLS` optimization step inside
  the `nlme` optimization. Default is 7.

- msMaxIter:

  maximum number of iterations for
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html) (`iter.max`) or the
  [`nlm`](https://rdrr.io/r/stats/nlm.html) (`iterlim`, from the 10-th
  step) optimization step inside the `nlme` optimization. Default is 50
  (which may be too small for e.g. for overparametrized cases).

- minScale:

  minimum factor by which to shrink the default step size in an attempt
  to decrease the sum of squares in the `PNLS` step. Default `0.001`.

- tolerance:

  tolerance for the convergence criterion in the `nlme` algorithm.
  Default is `1e-6`.

- niterEM:

  number of iterations for the EM algorithm used to refine the initial
  estimates of the random effects variance-covariance coefficients.
  Default is 25.

- pnlsTol:

  tolerance for the convergence criterion in `PNLS` step. Default is
  `1e-3`.

- msTol:

  tolerance for the convergence criterion in `nlm`, passed as the
  `gradtol` argument to the function (see documentation on `nlm`).
  Default is `1e-7`.

- returnObject:

  a logical value indicating whether the fitted object should be
  returned when the maximum number of iterations is reached without
  convergence of the algorithm. Default is `FALSE`.

- msVerbose:

  a logical value passed as the `trace` to
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html)`(.., control= list(trace = *, ..))`
  or as argument `print.level` to
  [`nlm()`](https://rdrr.io/r/stats/nlm.html). Default is `FALSE`.

- msWarnNoConv:

  logical indicating if a
  [`warning`](https://rdrr.io/r/base/warning.html) should be signalled
  whenever the minimization (by `opt`) in the LME step does not
  converge; defaults to `TRUE`.

- gradHess:

  a logical value indicating whether numerical gradient vectors and
  Hessian matrices of the log-likelihood function should be used in the
  `nlm` optimization. This option is only available when the correlation
  structure (`corStruct`) and the variance function structure
  (`varFunc`) have no "varying" parameters and the `pdMat` classes used
  in the random effects structure are `pdSymm` (general
  positive-definite), `pdDiag` (diagonal), `pdIdent` (multiple of the
  identity), or `pdCompSymm` (compound symmetry). Default is `TRUE`.

- apVar:

  a logical value indicating whether the approximate covariance matrix
  of the variance-covariance parameters should be calculated. Default is
  `TRUE`.

- .relStep:

  relative step for numerical derivatives calculations. Default is
  `.Machine$double.eps^(1/3)`.

- minAbsParApVar:

  numeric value - minimum absolute parameter value in the approximate
  variance calculation. The default is `0.05`.

- opt:

  the optimizer to be used, either
  `"`[`nlminb`](https://rdrr.io/r/stats/nlminb.html)`"` (the default) or
  `"`[`nlm`](https://rdrr.io/r/stats/nlm.html)`"`.

- natural:

  a logical value indicating whether the `pdNatural` parametrization
  should be used for general positive-definite matrices (`pdSymm`) in
  `reStruct`, when the approximate covariance matrix of the estimators
  is calculated. Default is `TRUE`.

- sigma:

  optionally a positive number to fix the residual error at. If `NULL`,
  as by default, or `0`, sigma is estimated.

- optExpression:

  Optimize the rxode2 expression to speed up calculation. By default
  this is turned on.

- literalFix:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- sumProd:

  Is a boolean indicating if the model should change multiplication to
  high precision multiplication and sums to high precision sums using
  the PreciseSums package. By default this is `FALSE`.

- rxControl:

  \`rxode2\` ODE solving options during fitting, created with
  \`rxControl()\`

- method:

  a character string. If `"REML"` the model is fit by maximizing the
  restricted log-likelihood. If `"ML"` the log-likelihood is maximized.
  Defaults to `"ML"`.

- random:

  optionally, any of the following: (i) a two-sided formula of the form
  `r1+...+rn~x1+...+xm | g1/.../gQ`, with `r1,...,rn` naming parameters
  included on the right hand side of `model`, `x1+...+xm` specifying the
  random-effects model for these parameters and `g1/.../gQ` the grouping
  structure (`Q` may be equal to 1, in which case no `/` is required).
  The random effects formula will be repeated for all levels of
  grouping, in the case of multiple levels of grouping; (ii) a two-sided
  formula of the form `r1+...+rn~x1+..+xm`, a list of two-sided formulas
  of the form `r1~x1+...+xm`, with possibly different random-effects
  models for different parameters, a `pdMat` object with a two-sided
  formula, or list of two-sided formulas (i.e. a non-`NULL` value for
  `formula(random)`), or a list of pdMat objects with two-sided
  formulas, or lists of two-sided formulas. In this case, the grouping
  structure formula will be given in `groups`, or derived from the data
  used to fit the nonlinear mixed-effects model, which should inherit
  from class `groupedData`,; (iii) a named list of formulas, lists of
  formulas, or `pdMat` objects as in (ii), with the grouping factors as
  names. The order of nesting will be assumed the same as the order of
  the order of the elements in the list; (iv) an `reStruct` object. See
  the documentation on
  [`pdClasses`](https://rdrr.io/pkg/nlme/man/pdClasses.html) for a
  description of the available `pdMat` classes. Defaults to `fixed`,
  resulting in all fixed effects having also random effects.

- fixed:

  a two-sided linear formula of the form `f1+...+fn~x1+...+xm`, or a
  list of two-sided formulas of the form `f1~x1+...+xm`, with possibly
  different models for different parameters. The `f1,...,fn` are the
  names of parameters included on the right hand side of `model` and the
  `x1+...+xm` expressions define linear models for these parameters
  (when the left hand side of the formula contains several parameters,
  they all are assumed to follow the same linear model, described by the
  right hand side expression). A `1` on the right hand side of the
  formula(s) indicates a single fixed effects for the corresponding
  parameter(s).

- weights:

  an optional `varFunc` object or one-sided formula describing the
  within-group heteroscedasticity structure. If given as a formula, it
  is used as the argument to `varFixed`, corresponding to fixed variance
  weights. See the documentation on
  [`varClasses`](https://rdrr.io/pkg/nlme/man/varClasses.html) for a
  description of the available `varFunc` classes. Defaults to `NULL`,
  corresponding to homoscedastic within-group errors.

- verbose:

  an optional logical value. If `TRUE` information on the evolution of
  the iterative algorithm is printed. Default is `FALSE`.

- returnNlme:

  Returns the nlme object instead of the nlmixr object (by default
  FALSE). If any of the nlme specific options of \`random\`, \`fixed\`,
  \`sens\`, the nlme object is returned

- addProp:

  Type of additive-plus-proportional error: \`"combined1"\`, where
  standard deviations add: \$\$y = f + (a + b\times f^c) \times
  \varepsilon\$\$; or \`"combined2"\`, where variances add: \$\$y = f +
  \sqrt{a^2 + b^2\times f^{2\times c}} \times \varepsilon\$\$. Here y =
  observed, f = predicted, a = additive sd, b = proportional/power sd, c
  = power exponent (1 in the proportional case).

- calcTables:

  This boolean is to determine if the foceiFit will calculate tables. By
  default this is `TRUE`

- compress:

  Should the object have compressed items

- adjObf:

  is a boolean to indicate if the objective function should be adjusted
  to be closer to NONMEM's default objective function. By default this
  is `TRUE`

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- sigdig:

  Optimization significant digits; controls the inner/outer optimization
  tolerance (`10^-sigdig`), the boundary check tolerance
  (`5*10^(-sigdig+1)`), and the ODE solver tolerances. The solver
  tolerances are split by solver stiffness and keep `atol` well below
  `rtol`: a stiff solver (`lsoda`/`liblsoda` – and any auto-switching
  method – plus `indLin`, the default) uses `rtol = 10^(-sigdig-3)`,
  `atol = 10^(-sigdig-5)`, while the non-stiff explicit `dop853` uses
  the looser `rtol = 10^-sigdig`, `atol = 10^(-sigdig-3)`. The
  sensitivity (`atolSens`/`rtolSens`) and steady-state
  (`ssAtol`/`ssRtol`) tolerances run one order looser than the
  corresponding main tolerance. At the default `sigdig = 4` a stiff
  solve is `atol = 1e-9`, `rtol = 1e-7`.

- sigdigTable:

  Significant digits in the final output table. If not specified, then
  it matches the significant digits in the \`sigdig\` optimization
  algorithm. If \`sigdig\` is NULL, use 3.

- muRefCovAlg:

  When \`TRUE\` (default), algebraic expressions that can be
  mu-referenced are internally rewritten as mu-referenced covariates and
  restored after optimization. Mirrors
  `saemControl(muRefCovAlg=)`/`nlmeControl(muRefCovAlg=)`; for
  [`foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.md)
  only takes effect when `muModel != "none"`.

- eventSens:

  Controls how dosing/event-parameter (\`alag\`, \`F\`, \`rate\`,
  \`dur\`) sensitivities are computed for THETA/ETA gradients:
  \`"jump"\` (default) uses rxode2's analytic event sensitivities;
  \`"fd"\` uses the legacy finite-difference behavior.

- print:

  Convenience alias for the shared nlmixr \`print\` control. \`nlme\`
  prints progress through its own \`verbose\` option, so \`print\` maps
  to it: \`print=0\` runs quietly (\`verbose=FALSE\`) and any positive
  value is verbose (\`verbose=TRUE\`). When \`print\` is not supplied an
  explicit \`verbose\` is used as given.

- covMethod:

  Covariance method: \`"analytic"\` (default) computes the focei
  observed-information covariance at the converged nlme estimates
  post-fit (falling back to the finite-difference \`"r,s"\` -\>
  \`"r"\`/\`"s"\` chain when out of analytic scope); \`"r,s"\`, \`"r"\`,
  \`"s"\` request the finite-difference forms directly; \`"nlme"\` and
  \`""\` skip the recompute and keep nlme's own standard errors. When
  the recompute fails the \`"nlme"\` covariance is kept.

- ...:

  Further, named control arguments to be passed to
  [`nlminb`](https://rdrr.io/r/stats/nlminb.html) (apart from `trace`
  and `iter.max` mentioned above), where used (`eval.max` and those from
  `abs.tol` down).

## Value

a nlmixr-nlme list

## See also

Other Estimation control:
[`foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.md),
[`saemControl()`](https://nlmixr2.github.io/nlmixr2est/reference/saemControl.md)

## Examples

``` r
nlmeControl()
#> $maxIter
#> [1] 100
#> 
#> $pnlsMaxIter
#> [1] 100
#> 
#> $msMaxIter
#> [1] 100
#> 
#> $minScale
#> [1] 0.001
#> 
#> $tolerance
#> [1] 1e-05
#> 
#> $niterEM
#> [1] 25
#> 
#> $pnlsTol
#> [1] 0.001
#> 
#> $msTol
#> [1] 1e-06
#> 
#> $returnObject
#> [1] FALSE
#> 
#> $msVerbose
#> [1] FALSE
#> 
#> $msWarnNoConv
#> [1] TRUE
#> 
#> $gradHess
#> [1] TRUE
#> 
#> $apVar
#> [1] TRUE
#> 
#> $.relStep
#> [1] 6.055454e-06
#> 
#> $minAbsParApVar
#> [1] 0.05
#> 
#> $opt
#> [1] "nlminb"
#> 
#> $natural
#> [1] TRUE
#> 
#> $sigma
#> [1] 0
#> 
#> $optExpression
#> [1] TRUE
#> 
#> $literalFix
#> [1] TRUE
#> 
#> $sumProd
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
#> [1] 1e-09
#> 
#> $rtol
#> [1] 1e-07
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
#> $method
#> [1] "ML"
#> 
#> $verbose
#> [1] TRUE
#> 
#> $returnNlme
#> [1] FALSE
#> 
#> $addProp
#> [1] "combined2"
#> 
#> $calcTables
#> [1] TRUE
#> 
#> $compress
#> [1] TRUE
#> 
#> $random
#> NULL
#> 
#> $fixed
#> NULL
#> 
#> $weights
#> NULL
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
#> $muRefCovAlg
#> [1] TRUE
#> 
#> $eventSens
#> [1] "jump"
#> 
#> $covMethod
#> [1] "nlme"
#> 
#> $genRxControl
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "nlmeControl"
nlmixr2NlmeControl()
#> $maxIter
#> [1] 100
#> 
#> $pnlsMaxIter
#> [1] 100
#> 
#> $msMaxIter
#> [1] 100
#> 
#> $minScale
#> [1] 0.001
#> 
#> $tolerance
#> [1] 1e-05
#> 
#> $niterEM
#> [1] 25
#> 
#> $pnlsTol
#> [1] 0.001
#> 
#> $msTol
#> [1] 1e-06
#> 
#> $returnObject
#> [1] FALSE
#> 
#> $msVerbose
#> [1] FALSE
#> 
#> $msWarnNoConv
#> [1] TRUE
#> 
#> $gradHess
#> [1] TRUE
#> 
#> $apVar
#> [1] TRUE
#> 
#> $.relStep
#> [1] 6.055454e-06
#> 
#> $minAbsParApVar
#> [1] 0.05
#> 
#> $opt
#> [1] "nlminb"
#> 
#> $natural
#> [1] TRUE
#> 
#> $sigma
#> [1] 0
#> 
#> $optExpression
#> [1] TRUE
#> 
#> $literalFix
#> [1] TRUE
#> 
#> $sumProd
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
#> [1] 1e-09
#> 
#> $rtol
#> [1] 1e-07
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
#> $method
#> [1] "ML"
#> 
#> $verbose
#> [1] TRUE
#> 
#> $returnNlme
#> [1] FALSE
#> 
#> $addProp
#> [1] "combined2"
#> 
#> $calcTables
#> [1] TRUE
#> 
#> $compress
#> [1] TRUE
#> 
#> $random
#> NULL
#> 
#> $fixed
#> NULL
#> 
#> $weights
#> NULL
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
#> $muRefCovAlg
#> [1] TRUE
#> 
#> $eventSens
#> [1] "jump"
#> 
#> $covMethod
#> [1] "nlme"
#> 
#> $genRxControl
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "nlmeControl"
```
