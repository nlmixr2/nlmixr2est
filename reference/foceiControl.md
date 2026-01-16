# Control Options for FOCEi

Control Options for FOCEi

## Usage

``` r
foceiControl(
  sigdig = 3,
  ...,
  epsilon = NULL,
  maxInnerIterations = 1000,
  maxOuterIterations = 5000,
  n1qn1nsim = NULL,
  print = 1L,
  printNcol = floor((getOption("width") - 23)/12),
  scaleTo = 1,
  scaleObjective = 0,
  normType = c("rescale2", "mean", "rescale", "std", "len", "constant"),
  scaleType = c("nlmixr2", "norm", "mult", "multAdd"),
  scaleCmax = 1e+05,
  scaleCmin = 1e-05,
  scaleC = NULL,
  scaleC0 = 1e+05,
  derivEps = rep(20 * sqrt(.Machine$double.eps), 2),
  derivMethod = c("switch", "forward", "central"),
  derivSwitchTol = NULL,
  covDerivMethod = c("central", "forward"),
  covMethod = c("r,s", "r", "s", ""),
  hessEps = (.Machine$double.eps)^(1/3),
  hessEpsLlik = (.Machine$double.eps)^(1/3),
  optimHessType = c("central", "forward"),
  optimHessCovType = c("central", "forward"),
  eventType = c("central", "forward"),
  centralDerivEps = rep(20 * sqrt(.Machine$double.eps), 2),
  lbfgsLmm = 7L,
  lbfgsPgtol = 0,
  lbfgsFactr = NULL,
  eigen = TRUE,
  diagXform = c("sqrt", "log", "identity"),
  iovXform = c("sd", "var", "logsd", "logvar"),
  sumProd = FALSE,
  optExpression = TRUE,
  literalFix = TRUE,
  literalFixRes = TRUE,
  ci = 0.95,
  useColor = crayon::has_color(),
  boundTol = NULL,
  calcTables = TRUE,
  noAbort = TRUE,
  interaction = TRUE,
  cholSEtol = (.Machine$double.eps)^(1/3),
  cholAccept = 0.001,
  resetEtaP = 0.15,
  resetThetaP = 0.05,
  resetThetaFinalP = 0.15,
  diagOmegaBoundUpper = 5,
  diagOmegaBoundLower = 100,
  cholSEOpt = FALSE,
  cholSECov = FALSE,
  fo = FALSE,
  covTryHarder = FALSE,
  outerOpt = c("nlminb", "bobyqa", "lbfgsb3c", "L-BFGS-B", "mma", "lbfgsbLG", "slsqp",
    "Rvmmin"),
  innerOpt = c("n1qn1", "BFGS"),
  rhobeg = 0.2,
  rhoend = NULL,
  npt = NULL,
  rel.tol = NULL,
  x.tol = NULL,
  eval.max = 4000,
  iter.max = 2000,
  abstol = NULL,
  reltol = NULL,
  resetHessianAndEta = FALSE,
  stateTrim = Inf,
  shi21maxOuter = 0L,
  shi21maxInner = 20L,
  shi21maxInnerCov = 20L,
  shi21maxFD = 20L,
  gillK = 10L,
  gillStep = 4,
  gillFtol = 0,
  gillRtol = sqrt(.Machine$double.eps),
  gillKcov = 10L,
  gillKcovLlik = 10L,
  gillStepCovLlik = 4.5,
  gillStepCov = 2,
  gillFtolCov = 0,
  gillFtolCovLlik = 0,
  rmatNorm = TRUE,
  rmatNormLlik = TRUE,
  smatNorm = TRUE,
  smatNormLlik = TRUE,
  covGillF = TRUE,
  optGillF = TRUE,
  covSmall = 1e-05,
  adjLik = TRUE,
  gradTrim = Inf,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  gradCalcCentralSmall = 1e-04,
  gradCalcCentralLarge = 10000,
  etaNudge = qnorm(1 - 0.05/2)/sqrt(3),
  etaNudge2 = qnorm(1 - 0.05/2) * sqrt(3/5),
  nRetries = 3,
  seed = 42,
  resetThetaCheckPer = 0.1,
  etaMat = NULL,
  repeatGillMax = 1,
  stickyRecalcN = 4,
  gradProgressOfvTime = 10,
  addProp = c("combined2", "combined1"),
  badSolveObjfAdj = 100,
  compress = TRUE,
  rxControl = NULL,
  sigdigTable = NULL,
  fallbackFD = FALSE,
  smatPer = 0.6,
  sdLowerFact = 0.001,
  zeroGradFirstReset = TRUE,
  zeroGradRunReset = TRUE,
  zeroGradBobyqa = TRUE,
  mceta = -1L,
  nAGQ = 0,
  agqLow = -Inf,
  agqHi = Inf
)
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

  Ignored parameters

- epsilon:

  Precision of estimate for n1qn1 optimization.

- maxInnerIterations:

  Number of iterations for n1qn1 optimization.

- maxOuterIterations:

  Maximum number of L-BFGS-B optimization for outer problem.

- n1qn1nsim:

  Number of function evaluations for n1qn1 optimization.

- print:

  Integer representing when the outer step is printed. When this is 0 or
  do not print the iterations. 1 is print every function evaluation
  (default), 5 is print every 5 evaluations.

- printNcol:

  Number of columns to printout before wrapping parameter
  estimates/gradient

- scaleTo:

  Scale the initial parameter estimate to this value. By default this
  is 1. When zero or below, no scaling is performed.

- scaleObjective:

  Scale the initial objective function to this value. By default this is
  0 (meaning do not scale)

- normType:

  This is the type of parameter normalization/scaling used to get the
  scaled initial values for nlmixr2. These are used with `scaleType` of.

  With the exception of `rescale2`, these come from [Feature
  Scaling](https://en.wikipedia.org/wiki/Feature_scaling). The
  `rescale2` The rescaling is the same type described in the
  [OptdesX](http://apmonitor.com/me575/uploads/Main/optimization_book.pdf)
  software manual.

  In general, all all scaling formula can be described by:

  \$\$v\_{scaled}\$\$ = (\$\$v\_{unscaled}-C\_{1}\$\$)/\$\$C\_{2}\$\$

  Where

  The other data normalization approaches follow the following formula

  \$\$v\_{scaled}\$\$ = (\$\$v\_{unscaled}-C\_{1}\$\$)/\$\$C\_{2}\$\$

  - `rescale2` This scales all parameters from (-1 to 1). The relative
    differences between the parameters are preserved with this approach
    and the constants are:

    \$\$C\_{1}\$\$ = (max(all unscaled values)+min(all unscaled
    values))/2

    \$\$C\_{2}\$\$ = (max(all unscaled values) - min(all unscaled
    values))/2

  - `rescale` or min-max normalization. This rescales all parameters
    from (0 to 1). As in the `rescale2` the relative differences are
    preserved. In this approach:

    \$\$C\_{1}\$\$ = min(all unscaled values)

    \$\$C\_{2}\$\$ = max(all unscaled values) - min(all unscaled values)

  - `mean` or mean normalization. This rescales to center the parameters
    around the mean but the parameters are from 0 to 1. In this
    approach:

    \$\$C\_{1}\$\$ = mean(all unscaled values)

    \$\$C\_{2}\$\$ = max(all unscaled values) - min(all unscaled values)

  - `std` or standardization. This standardizes by the mean and standard
    deviation. In this approach:

    \$\$C\_{1}\$\$ = mean(all unscaled values)

    \$\$C\_{2}\$\$ = sd(all unscaled values)

  - `len` or unit length scaling. This scales the parameters to the unit
    length. For this approach we use the Euclidean length, that is:

    \$\$C\_{1}\$\$ = 0

    \$\$C\_{2}\$\$ = \$\$\sqrt(v_1^2 + v_2^2 + \cdots + v_n^2)\$\$

  - `constant` which does not perform data normalization. That is

    \$\$C\_{1}\$\$ = 0

    \$\$C\_{2}\$\$ = 1

- scaleType:

  The scaling scheme for nlmixr2. The supported types are:

  - `nlmixr2` In this approach the scaling is performed by the following
    equation:

    \$\$v\_{scaled}\$\$ = (\$\$v\_{current} -
    v\_{init}\$\$)\*scaleC\[i\] + scaleTo

    The `scaleTo` parameter is specified by the `normType`, and the
    scales are specified by `scaleC`.

  - `norm` This approach uses the simple scaling provided by the
    `normType` argument.

  - `mult` This approach does not use the data normalization provided by
    `normType`, but rather uses multiplicative scaling to a constant
    provided by the `scaleTo` argument.

    In this case:

    \$\$v\_{scaled}\$\$ =
    \$\$v\_{current}\$\$/\$\$v\_{init}\$\$\*scaleTo

  - `multAdd` This approach changes the scaling based on the parameter
    being specified. If a parameter is defined in an exponential block
    (ie exp(theta)), then it is scaled on a linearly, that is:

    \$\$v\_{scaled}\$\$ = (\$\$v\_{current}-v\_{init}\$\$) + scaleTo

    Otherwise the parameter is scaled multiplicatively.

    \$\$v\_{scaled}\$\$ =
    \$\$v\_{current}\$\$/\$\$v\_{init}\$\$\*scaleTo

- scaleCmax:

  Maximum value of the scaleC to prevent overflow.

- scaleCmin:

  Minimum value of the scaleC to prevent underflow.

- scaleC:

  The scaling constant used with `scaleType=nlmixr2`. When not
  specified, it is based on the type of parameter that is estimated. The
  idea is to keep the derivatives similar on a log scale to have similar
  gradient sizes. Hence parameters like log(exp(theta)) would have a
  scaling factor of 1 and log(theta) would have a scaling factor of
  ini_value (to scale by 1/value; ie d/dt(log(ini_value)) = 1/ini_value
  or scaleC=ini_value)

  - For parameters in an exponential (ie exp(theta)) or parameters
    specifying powers, boxCox or yeoJohnson transformations , this is 1.

  - For additive, proportional, lognormal error structures, these are
    given by 0.5\*abs(initial_estimate)

  - Factorials are scaled by abs(1/digamma(initial_estimate+1))

  - parameters in a log scale (ie log(theta)) are transformed by
    log(abs(initial_estimate))\*abs(initial_estimate)

  These parameter scaling coefficients are chose to try to keep similar
  slopes among parameters. That is they all follow the slopes
  approximately on a log-scale.

  While these are chosen in a logical manner, they may not always apply.
  You can specify each parameters scaling factor by this parameter if
  you wish.

- scaleC0:

  Number to adjust the scaling factor by if the initial gradient is
  zero.

- derivEps:

  Forward difference tolerances, which is a vector of relative
  difference and absolute difference. The central/forward difference
  step size h is calculated as:

  `h = abs(x)*derivEps[1] + derivEps[2]`

- derivMethod:

  indicates the method for calculating derivatives of the outer problem.
  Currently supports "switch", "central" and "forward" difference
  methods. Switch starts with forward differences. This will switch to
  central differences when abs(delta(OFV)) \<= derivSwitchTol and switch
  back to forward differences when abs(delta(OFV)) \> derivSwitchTol.

- derivSwitchTol:

  The tolerance to switch forward to central differences.

- covDerivMethod:

  indicates the method for calculating the derivatives while calculating
  the covariance components (Hessian and S).

- covMethod:

  Method for calculating covariance. In this discussion, R is the
  Hessian matrix of the objective function. The S matrix is the sum of
  individual gradient cross-product (evaluated at the individual
  empirical Bayes estimates).

  - "`r,s`" Uses the sandwich matrix to calculate the covariance, that
    is: `solve(R) %*% S %*% solve(R)`

  - "`r`" Uses the Hessian matrix to calculate the covariance as
    `2 %*% solve(R)`

  - "`s`" Uses the cross-product matrix to calculate the covariance as
    `4 %*% solve(S)`

  - "" Does not calculate the covariance step.

- hessEps:

  is a double value representing the epsilon for the Hessian
  calculation. This is used for the R matrix calculation.

- hessEpsLlik:

  is a double value representing the epsilon for the Hessian calculation
  when doing focei generalized log-likelihood estimation. This is used
  for the R matrix calculation.

- optimHessType:

  The hessian type for when calculating the individual hessian by
  numeric differences (in generalized log-likelihood estimation). The
  options are "central", and "forward". The central differences is what
  R's \`optimHess()\` uses and is the default for this method. (Though
  the "forward" is faster and still reasonable for most cases). The
  Shi21 cannot be changed for the Gill83 algorithm with the optimHess in
  a generalized likelihood problem.

- optimHessCovType:

  The hessian type for when calculating the individual hessian by
  numeric differences (in generalized log-likelihood estimation). The
  options are "central", and "forward". The central differences is what
  R's \`optimHess()\` uses. While this takes longer in optimization, it
  is more accurate, so for calculating the covariance and final
  likelihood, the central differences are used. This also uses the
  modified Shi21 method

- eventType:

  Event gradient type for dosing events; Can be "central" or "forward"

- centralDerivEps:

  Central difference tolerances. This is a numeric vector of relative
  difference and absolute difference. The central/forward difference
  step size h is calculated as:

  `h = abs(x)*derivEps[1] + derivEps[2]`

- lbfgsLmm:

  An integer giving the number of BFGS updates retained in the
  "L-BFGS-B" method, It defaults to 7.

- lbfgsPgtol:

  is a double precision variable.

  On entry pgtol \>= 0 is specified by the user. The iteration will stop
  when:

  `max(\| proj g_i \| i = 1, ..., n) <= lbfgsPgtol`

  where pg_i is the ith component of the projected gradient.

  On exit pgtol is unchanged. This defaults to zero, when the check is
  suppressed.

- lbfgsFactr:

  Controls the convergence of the "L-BFGS-B" method. Convergence occurs
  when the reduction in the objective is within this factor of the
  machine tolerance. Default is 1e10, which gives a tolerance of about
  `2e-6`, approximately 4 sigdigs. You can check your exact tolerance by
  multiplying this value by `.Machine$double.eps`

- eigen:

  A boolean indicating if eigenvectors are calculated to include a
  condition number calculation.

- diagXform:

  This is the transformation used on the diagonal of the
  `chol(solve(omega))`. This matrix and values are the parameters
  estimated in FOCEi. The possibilities are:

  - `sqrt` Estimates the sqrt of the diagonal elements of
    `chol(solve(omega))`. This is the default method.

  - `log` Estimates the log of the diagonal elements of
    `chol(solve(omega))`

  - `identity` Estimates the diagonal elements without any
    transformations

- iovXform:

  This is the transformation used on the diagonal of the \`iov\`. The
  possibilities are:

  - `sd` Estimate the IOV as the standard deviation for IOV

  - `var` Estimate the IOV as the variance for IOV.

  - `logsd` Estimate the IOV as the log(sd) instead of sd.

  - `logvar` Estimate the IOV as the log(var) instead of variance.

- sumProd:

  Is a boolean indicating if the model should change multiplication to
  high precision multiplication and sums to high precision sums using
  the PreciseSums package. By default this is `FALSE`.

- optExpression:

  Optimize the rxode2 expression to speed up calculation. By default
  this is turned on.

- literalFix:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- literalFixRes:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- useColor:

  Boolean indicating if focei can use ASCII color codes

- boundTol:

  Tolerance for boundary issues.

- calcTables:

  This boolean is to determine if the foceiFit will calculate tables. By
  default this is `TRUE`

- noAbort:

  Boolean to indicate if you should abort the FOCEi evaluation if it
  runs into troubles. (default TRUE)

- interaction:

  Boolean indicate FOCEi should be used (TRUE) instead of FOCE (FALSE)

- cholSEtol:

  tolerance for Generalized Cholesky Decomposition. Defaults to
  suggested (.Machine\$double.eps)^(1/3)

- cholAccept:

  Tolerance to accept a Generalized Cholesky Decomposition for a R or S
  matrix.

- resetEtaP:

  represents the p-value for reseting the individual ETA to 0 during
  optimization (instead of the saved value). The two test statistics
  used in the z-test are either chol(omega^-1) %\*% eta or
  eta/sd(allEtas). A p-value of 0 indicates the ETAs never reset. A
  p-value of 1 indicates the ETAs always reset.

- resetThetaP:

  represents the p-value for reseting the population mu-referenced THETA
  parameters based on ETA drift during optimization, and resetting the
  optimization. A p-value of 0 indicates the THETAs never reset. A
  p-value of 1 indicates the THETAs always reset and is not allowed. The
  theta reset is checked at the beginning and when nearing a local
  minima. The percent change in objective function where a theta reset
  check is initiated is controlled in `resetThetaCheckPer`.

- resetThetaFinalP:

  represents the p-value for reseting the population mu-referenced THETA
  parameters based on ETA drift during optimization, and resetting the
  optimization one final time.

- diagOmegaBoundUpper:

  This represents the upper bound of the diagonal omega matrix. The
  upper bound is given by diag(omega)\*diagOmegaBoundUpper. If
  `diagOmegaBoundUpper` is 1, there is no upper bound on Omega.

- diagOmegaBoundLower:

  This represents the lower bound of the diagonal omega matrix. The
  lower bound is given by diag(omega)/diagOmegaBoundUpper. If
  `diagOmegaBoundLower` is 1, there is no lower bound on Omega.

- cholSEOpt:

  Boolean indicating if the generalized Cholesky should be used while
  optimizing.

- cholSECov:

  Boolean indicating if the generalized Cholesky should be used while
  calculating the Covariance Matrix.

- fo:

  is a boolean indicating if this is a FO approximation routine.

- covTryHarder:

  If the R matrix is non-positive definite and cannot be corrected to be
  non-positive definite try estimating the Hessian on the unscaled
  parameter space.

- outerOpt:

  optimization method for the outer problem

- innerOpt:

  optimization method for the inner problem (not implemented yet.)

- rhobeg:

  Beginning change in parameters for bobyqa algorithm (trust region). By
  default this is 0.2 or 20 parameters when the parameters are scaled
  to 1. rhobeg and rhoend must be set to the initial and final values of
  a trust region radius, so both must be positive with 0 \< rhoend \<
  rhobeg. Typically rhobeg should be about one tenth of the greatest
  expected change to a variable. Note also that smallest difference
  abs(upper-lower) should be greater than or equal to rhobeg\*2. If this
  is not the case then rhobeg will be adjusted. (bobyqa)

- rhoend:

  The smallest value of the trust region radius that is allowed. If not
  defined, then 10^(-sigdig-1) will be used. (bobyqa)

- npt:

  The number of points used to approximate the objective function via a
  quadratic approximation for bobyqa. The value of npt must be in the
  interval \[n+2,(n+1)(n+2)/2\] where n is the number of parameters in
  par. Choices that exceed 2\*n+1 are not recommended. If not defined,
  it will be set to 2\*n + 1. (bobyqa)

- rel.tol:

  Relative tolerance before nlminb stops (nlmimb).

- x.tol:

  X tolerance for nlmixr2 optimizer

- eval.max:

  Number of maximum evaluations of the objective function (nlmimb)

- iter.max:

  Maximum number of iterations allowed (nlmimb)

- abstol:

  Absolute tolerance for nlmixr2 optimizer (BFGS)

- reltol:

  tolerance for nlmixr2 (BFGS)

- resetHessianAndEta:

  is a boolean representing if the individual Hessian is reset when ETAs
  are reset using the option `resetEtaP`.

- stateTrim:

  Trim state amounts/concentrations to this value.

- shi21maxOuter:

  The maximum number of steps for the optimization of the
  forward-difference step size. When not zero, use this instead of Gill
  differences.

- shi21maxInner:

  The maximum number of steps for the optimization of the individual
  Hessian matrices in the generalized likelihood problem. When 0,
  un-optimized finite differences are used.

- shi21maxInnerCov:

  The maximum number of steps for the optimization of the individual
  Hessian matrices in the generalized likelihood problem for the
  covariance step. When 0, un-optimized finite differences are used.

- shi21maxFD:

  The maximum number of steps for the optimization of the forward
  difference step size when using dosing events (lag time, modeled
  duration/rate and bioavailability)

- gillK:

  The total number of possible steps to determine the optimal
  forward/central difference step size per parameter (by the Gill 1983
  method). If 0, no optimal step size is determined. Otherwise this is
  the optimal step size determined.

- gillStep:

  When looking for the optimal forward difference step size, this is
  This is the step size to increase the initial estimate by. So each
  iteration the new step size = (prior step size)\*gillStep

- gillFtol:

  The gillFtol is the gradient error tolerance that is acceptable before
  issuing a warning/error about the gradient estimates.

- gillRtol:

  The relative tolerance used for Gill 1983 determination of optimal
  step size.

- gillKcov:

  The total number of possible steps to determine the optimal
  forward/central difference step size per parameter (by the Gill 1983
  method) during the covariance step. If 0, no optimal step size is
  determined. Otherwise this is the optimal step size determined.

- gillKcovLlik:

  The total number of possible steps to determine the optimal
  forward/central difference step per parameter when using the
  generalized focei log-likelihood method (by the Gill 1986 method). If
  0, no optimal step size is determined. Otherwise this is the optimal
  step size is determined

- gillStepCovLlik:

  Same as above but during generalized focei log-likelihood

- gillStepCov:

  When looking for the optimal forward difference step size, this is
  This is the step size to increase the initial estimate by. So each
  iteration during the covariance step is equal to the new step size =
  (prior step size)\*gillStepCov

- gillFtolCov:

  The gillFtol is the gradient error tolerance that is acceptable before
  issuing a warning/error about the gradient estimates during the
  covariance step.

- gillFtolCovLlik:

  Same as above but applied during generalized log-likelihood
  estimation.

- rmatNorm:

  A parameter to normalize gradient step size by the parameter value
  during the calculation of the R matrix

- rmatNormLlik:

  A parameter to normalize gradient step size by the parameter value
  during the calculation of the R matrix if you are using generalized
  log-likelihood Hessian matrix.

- smatNorm:

  A parameter to normalize gradient step size by the parameter value
  during the calculation of the S matrix

- smatNormLlik:

  A parameter to normalize gradient step size by the parameter value
  during the calculation of the S matrix if you are using the
  generalized log-likelihood.

- covGillF:

  Use the Gill calculated optimal Forward difference step size for the
  instead of the central difference step size during the central
  difference gradient calculation.

- optGillF:

  Use the Gill calculated optimal Forward difference step size for the
  instead of the central difference step size during the central
  differences for optimization.

- covSmall:

  The covSmall is the small number to compare covariance numbers before
  rejecting an estimate of the covariance as the final estimate (when
  comparing sandwich vs R/S matrix estimates of the covariance). This
  number controls how small the variance is before the covariance matrix
  is rejected.

- adjLik:

  In nlmixr2, the objective function matches NONMEM's objective
  function, which removes a 2\*pi constant from the likelihood
  calculation. If this is TRUE, the likelihood function is adjusted by
  this 2\*pi factor. When adjusted this number more closely matches the
  likelihood approximations of nlme, and SAS approximations. Regardless
  of if this is turned on or off the objective function matches NONMEM's
  objective function.

- gradTrim:

  The parameter to adjust the gradient to if the \|gradient\| is very
  large.

- maxOdeRecalc:

  Maximum number of times to reduce the ODE tolerances and try to
  resolve the system if there was a bad ODE solve.

- odeRecalcFactor:

  The ODE recalculation factor when ODE solving goes bad, this is the
  factor the rtol/atol is reduced

- gradCalcCentralSmall:

  A small number that represents the value where \|grad\| \<
  gradCalcCentralSmall where forward differences switch to central
  differences.

- gradCalcCentralLarge:

  A large number that represents the value where \|grad\| \>
  gradCalcCentralLarge where forward differences switch to central
  differences.

- etaNudge:

  By default initial ETA estimates start at zero; Sometimes this doesn't
  optimize appropriately. If this value is non-zero, when the n1qn1
  optimization didn't perform appropriately, reset the Hessian, and
  nudge the ETA up by this value; If the ETA still doesn't move, nudge
  the ETA down by this value. By default this value is
  qnorm(1-0.05/2)\*1/sqrt(3), the first of the Gauss Quadrature numbers
  times by the 0.95% normal region. If this is not successful try the
  second eta nudge number (below). If +-etaNudge2 is not successful,
  then assign to zero and do not optimize any longer

- etaNudge2:

  This is the second eta nudge. By default it is
  qnorm(1-0.05/2)\*sqrt(3/5), which is the n=3 quadrature point
  (excluding zero) times by the 0.95% normal region

- nRetries:

  If FOCEi doesn't fit with the current parameter estimates, randomly
  sample new parameter estimates and restart the problem. This is
  similar to 'PsN' resampling.

- seed:

  an object specifying if and how the random number generator should be
  initialized

- resetThetaCheckPer:

  represents objective function % percentage below which resetThetaP is
  checked.

- etaMat:

  Eta matrix for initial estimates or final estimates of the ETAs.

  This can also be a fit to take use the final estimation estimates and
  use them as the initial eta value of the next fit.

  By default, it will be the initial values of the etas from the last
  fit (if supplied) or missing, meaning all ETAs start at zero
  (\`NULL\`)

  When this value is \`NA\`, the initial ETA estimates are not taken
  from the last fit.

- repeatGillMax:

  If the tolerances were reduced when calculating the initial Gill
  differences, the Gill difference is repeated up to a maximum number of
  times defined by this parameter.

- stickyRecalcN:

  The number of bad ODE solves before reducing the atol/rtol for the
  rest of the problem.

- gradProgressOfvTime:

  This is the time for a single objective function evaluation (in
  seconds) to start progress bars on gradient evaluations

- addProp:

  specifies the type of additive plus proportional errors, the one where
  standard deviations add (combined1) or the type where the variances
  add (combined2).

  The combined1 error type can be described by the following equation:

  \$\$y = f + (a + b\times f^c) \times \varepsilon\$\$

  The combined2 error model can be described by the following equation:

  \$\$y = f + \sqrt{a^2 + b^2\times f^{2\times c}} \times
  \varepsilon\$\$

  Where:

  \- y represents the observed value

  \- f represents the predicted value

  \- a is the additive standard deviation

  \- b is the proportional/power standard deviation

  \- c is the power exponent (in the proportional case c=1)

- badSolveObjfAdj:

  The objective function adjustment when the ODE system cannot be
  solved. It is based on each individual bad solve.

- compress:

  Should the object have compressed items

- rxControl:

  \`rxode2\` ODE solving options during fitting, created with
  \`rxControl()\`

- sigdigTable:

  Significant digits in the final output table. If not specified, then
  it matches the significant digits in the \`sigdig\` optimization
  algorithm. If \`sigdig\` is NULL, use 3.

- fallbackFD:

  Fallback to the finite differences if the sensitivity equations do not
  solve.

- smatPer:

  A percentage representing the number of failed parameter gradients for
  each individual (which are replaced with the overall gradient for the
  parameter) out of the total number of gradients parameters (ie
  \`ntheta\*nsub\`) before the S matrix is considered to be a bad
  matrix.

- sdLowerFact:

  A factor for multiplying the estimate by when the lower estimate is
  zero and the error is known to represent a standard deviation of a
  parameter (like add.sd, prop.sd, pow.sd, lnorm.sd, etc). When zero, no
  factor is applied. If your initial estimate is 0.15 and your lower
  bound is zero, then the lower bound would be assumed to be 0.00015.

- zeroGradFirstReset:

  boolean, when \`TRUE\` if the first gradient is zero, reset the zero
  gradient to \`sqrt(.Machine\$double.eps)\` to get past the bad initial
  estimate, otherwise error (and possibly reset), when \`FALSE\` error
  when the first gradient is zero. When \`NA\` on the last reset, have
  the zero gradient ignored, otherwise error and look for another value.
  Default is \`TRUE\`

- zeroGradRunReset:

  boolean, when \`TRUE\` if a gradient is zero, reset the zero gradient
  to \`sqrt(.Machine\$double.eps)\` to get past the bad estimate while
  running. Otherwise error (and possibly reset). Default is \`TRUE\`

- zeroGradBobyqa:

  boolean, when \`TRUE\` if a gradient is zero, the reset will change
  the method to the gradient free bobyqa method. When \`NA\`, the zero
  gradient will change to bobyqa only when the first gradient is zero.
  Default is \`TRUE\`

- mceta:

  Integer indicating the type of Monte Carlo sampling to perform for the
  best initial ETA estimate (based on \`omega\`). When:

  \- \`-1\` the last eta is used for the optimization (default)

  \- \`0\` eta=0 is used for each inner optimization

  For the rest of the \`mceta\`, each parameter's inner objective
  function is calculated and the eta set with the best objective
  function is used. With these further options:

  \- \`1\` the last eta and eta=0 are used

  \- \`2\` the last eta and eta=0 are used, as well as 1 sampled eta
  from the omega matrix

  \- \`n\` the last eta and eta=0 are used, as well as n-1 sampled etas
  from the omega matrix

- nAGQ:

  Number of Gauss-Hermite Adaptive Quadrature points to take. When
  \`nAGQ=0\`, the AGQ is not used. With \`nAGQ=1\`, this is equivalent
  to the Laplace method. The adaptive quadrature expands every node for
  each of the ETAs, so it can be quite expensive with a large amount of
  ETAs. Once the EBE is obtained for a subject, you will have nAGQ^neta
  additional function evaluations for even nAGQ numbers and
  (nAGQ^neta)-1 additional function evaluations for odd nAGQ numbers.

- agqLow:

  The lower bound for adaptive quadrature log-likelihood. By default
  this is -Inf; in the original nlmixr's gnlmm it was -700.

- agqHi:

  The upper bound for adaptive quadrature log-likelihood. By default
  this is Inf; in the original nlmixr's gnlmm was 400.

## Value

The control object that changes the options for the FOCEi family of
estimation methods

## Details

Note this uses the R's L-BFGS-B in
[`optim`](https://rdrr.io/r/stats/optim.html) for the outer problem and
the BFGS [`n1qn1`](https://rdrr.io/pkg/n1qn1/man/n1qn1.html) with that
allows restoring the prior individual Hessian (for faster optimization
speed).

However the inner problem is not scaled. Since most eta estimates start
near zero, scaling for these parameters do not make sense.

This process of scaling can fix some ill conditioning for the unscaled
problem. The covariance step is performed on the unscaled problem, so
the condition number of that matrix may not be reflective of the scaled
problem's condition-number.

## References

Gill, P.E., Murray, W., Saunders, M.A., & Wright, M.H. (1983). Computing
Forward-Difference Intervals for Numerical Optimization. Siam Journal on
Scientific and Statistical Computing, 4, 310-321.

Shi, H.M., Xie, Y., Xuan, M.Q., & Nocedal, J. (2021). Adaptive
Finite-Difference Interval Estimation for Noisy Derivative-Free
Optimization.

## See also

[`optim`](https://rdrr.io/r/stats/optim.html)

[`n1qn1`](https://rdrr.io/pkg/n1qn1/man/n1qn1.html)

[`rxSolve`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)

Other Estimation control:
[`nlmixr2NlmeControl()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2NlmeControl.md),
[`saemControl()`](https://nlmixr2.github.io/nlmixr2est/reference/saemControl.md)

## Author

Matthew L. Fidler
