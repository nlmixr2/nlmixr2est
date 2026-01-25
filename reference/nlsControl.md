# nlmixr2 defaults controls for nls

nlmixr2 defaults controls for nls

## Usage

``` r
nlsControl(
  maxiter = 10000,
  tol = 1e-05,
  minFactor = 1/1024,
  printEval = FALSE,
  warnOnly = FALSE,
  scaleOffset = 0,
  nDcentral = FALSE,
  algorithm = c("LM", "default", "plinear", "port"),
  ftol = sqrt(.Machine$double.eps),
  ptol = sqrt(.Machine$double.eps),
  gtol = 0,
  diag = list(),
  epsfcn = 0,
  factor = 100,
  maxfev = integer(),
  nprint = 0,
  solveType = c("grad", "fun"),
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  eventType = c("central", "forward"),
  shiErr = (.Machine$double.eps)^(1/3),
  shi21maxFD = 20L,
  useColor = crayon::has_color(),
  printNcol = floor((getOption("width") - 23)/12),
  print = 1L,
  normType = c("rescale2", "mean", "rescale", "std", "len", "constant"),
  scaleType = c("nlmixr2", "norm", "mult", "multAdd"),
  scaleCmax = 1e+05,
  scaleCmin = 1e-05,
  scaleC = NULL,
  scaleTo = 1,
  gradTo = 1,
  trace = FALSE,
  rxControl = NULL,
  optExpression = TRUE,
  sumProd = FALSE,
  literalFix = TRUE,
  returnNls = FALSE,
  addProp = c("combined2", "combined1"),
  calcTables = TRUE,
  compress = TRUE,
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  ...
)
```

## Arguments

- maxiter:

  A positive integer specifying the maximum number of iterations
  allowed.

- tol:

  A positive numeric value specifying the tolerance level for the
  relative offset convergence criterion.

- minFactor:

  A positive numeric value specifying the minimum step-size factor
  allowed on any step in the iteration. The increment is calculated with
  a Gauss-Newton algorithm and successively halved until the residual
  sum of squares has been decreased or until the step-size factor has
  been reduced below this limit.

- printEval:

  a logical specifying whether the number of evaluations (steps in the
  gradient direction taken each iteration) is printed.

- warnOnly:

  a logical specifying whether
  [`nls()`](https://rdrr.io/r/stats/nls.html) should return instead of
  signalling an error in the case of termination before convergence.
  Termination before convergence happens upon completion of `maxiter`
  iterations, in the case of a singular gradient, and in the case that
  the step-size factor is reduced below `minFactor`.

- scaleOffset:

  a constant to be added to the denominator of the relative offset
  convergence criterion calculation to avoid a zero divide in the case
  where the fit of a model to data is very close. The default value of
  `0` keeps the legacy behaviour of
  [`nls()`](https://rdrr.io/r/stats/nls.html). A value such as `1` seems
  to work for problems of reasonable scale with very small residuals.

- nDcentral:

  only when *numerical* derivatives are used:
  [`logical`](https://rdrr.io/r/base/logical.html) indicating if
  *central* differences should be employed, i.e.,
  [`numericDeriv`](https://rdrr.io/r/stats/numericDeriv.html)`(*, central=TRUE)`
  be used.

- algorithm:

  character string specifying the algorithm to use. The default
  algorithm is a Gauss-Newton algorithm. Other possible values are
  `"plinear"` for the Golub-Pereyra algorithm for partially linear
  least-squares models and `"port"` for the ‘nl2sol’ algorithm from the
  Port library – see the references. Can be abbreviated.

- ftol:

  non-negative numeric. Termination occurs when both the actual and
  predicted relative reductions in the sum of squares are at most
  `ftol`. Therefore, `ftol` measures the relative error desired in the
  sum of squares.

- ptol:

  non-negative numeric. Termination occurs when the relative error
  between two consecutive iterates is at most `ptol`. Therefore, `ptol`
  measures the relative error desired in the approximate solution.

- gtol:

  non-negative numeric. Termination occurs when the cosine of the angle
  between result of `fn` evaluation \\fvec\\ and any column of the
  Jacobian is at most `gtol` in absolute value. Therefore, `gtol`
  measures the orthogonality desired between the function vector and the
  columns of the Jacobian.

- diag:

  a list or numeric vector containing positive entries that serve as
  multiplicative scale factors for the parameters. Length of `diag`
  should be equal to that of `par`. If not, user-provided `diag` is
  ignored and `diag` is internally set.

- epsfcn:

  (used if `jac` is not provided) is a numeric used in determining a
  suitable step for the forward-difference approximation. This
  approximation assumes that the relative errors in the functions are of
  the order of `epsfcn`. If `epsfcn` is less than the machine precision,
  it is assumed that the relative errors in the functions are of the
  order of the machine precision.

- factor:

  positive numeric, used in determining the initial step bound. This
  bound is set to the product of `factor` and the
  \\\|\code{diag}\*\code{par}\|\\ if nonzero, or else to `factor`
  itself. In most cases `factor` should lie in the interval (0.1,100).
  100 is a generally recommended value.

- maxfev:

  integer; termination occurs when the number of calls to `fn` has
  reached `maxfev`. Note that `nls.lm` sets the value of `maxfev` to
  `100*(length(par) + 1)` if `maxfev = integer()`, where `par` is the
  list or vector of parameters to be optimized.

- nprint:

  is an integer; set `nprint` to be positive to enable printing of
  iterates

- solveType:

  tells if \`nlm\` will use nlmixr2's analytical gradients when
  available (finite differences will be used for event-related
  parameters like parameters controlling lag time, duration/rate of
  infusion, and modeled bioavailability). This can be:

  \- \`"hessian"\` which will use the analytical gradients to create a
  Hessian with finite differences.

  \- \`"gradient"\` which will use the gradient and let \`nlm\`
  calculate the finite difference hessian

  \- \`"fun"\` where nlm will calculate both the finite difference
  gradient and the finite difference Hessian

  When using nlmixr2's finite differences, the "ideal" step size for
  either central or forward differences are optimized for with the
  Shi2021 method which may give more accurate derivatives

- stickyRecalcN:

  The number of bad ODE solves before reducing the atol/rtol for the
  rest of the problem.

- maxOdeRecalc:

  Maximum number of times to reduce the ODE tolerances and try to
  resolve the system if there was a bad ODE solve.

- odeRecalcFactor:

  The ODE recalculation factor when ODE solving goes bad, this is the
  factor the rtol/atol is reduced

- eventType:

  Event gradient type for dosing events; Can be "central" or "forward"

- shiErr:

  This represents the epsilon when optimizing the ideal step size for
  numeric differentiation using the Shi2021 method

- shi21maxFD:

  The maximum number of steps for the optimization of the forward
  difference step size when using dosing events (lag time, modeled
  duration/rate and bioavailability)

- useColor:

  Boolean indicating if focei can use ASCII color codes

- printNcol:

  Number of columns to printout before wrapping parameter
  estimates/gradient

- print:

  Integer representing when the outer step is printed. When this is 0 or
  do not print the iterations. 1 is print every function evaluation
  (default), 5 is print every 5 evaluations.

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

- scaleTo:

  Scale the initial parameter estimate to this value. By default this
  is 1. When zero or below, no scaling is performed.

- gradTo:

  this is the factor that the gradient is scaled to before optimizing.
  This only works with scaleType="nlmixr2".

- trace:

  logical value indicating if a trace of the iteration progress should
  be printed. Default is `FALSE`. If `TRUE` the residual (weighted)
  sum-of-squares, the convergence criterion and the parameter values are
  printed at the conclusion of each iteration. Note that
  [`format()`](https://rdrr.io/r/base/format.html) is used, so these
  mostly depend on
  [`getOption`](https://rdrr.io/r/base/options.html)`("digits")`. When
  the `"plinear"` algorithm is used, the conditional estimates of the
  linear parameters are printed after the nonlinear parameters. When the
  `"port"` algorithm is used the objective function value printed is
  half the residual (weighted) sum-of-squares.

- rxControl:

  \`rxode2\` ODE solving options during fitting, created with
  \`rxControl()\`

- optExpression:

  Optimize the rxode2 expression to speed up calculation. By default
  this is turned on.

- sumProd:

  Is a boolean indicating if the model should change multiplication to
  high precision multiplication and sums to high precision sums using
  the PreciseSums package. By default this is `FALSE`.

- literalFix:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- returnNls:

  logical; when TRUE, will return the nls object instead of the nlmixr
  object

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

  Optimization significant digits. This controls:

  - The tolerance of the inner and outer optimization is `10^-sigdig`

  - The tolerance of the ODE solvers is `0.5*10^(-sigdig-2)`; For the
    sensitivity equations and steady-state solutions the default is
    `0.5*10^(-sigdig-1.5)` (sensitivity changes only applicable for
    liblsoda)

  - The tolerance of the boundary check is `5 * 10 ^ (-sigdig + 1)`

- sigdigTable:

  Significant digits in the final output table. If not specified, then
  it matches the significant digits in the \`sigdig\` optimization
  algorithm. If \`sigdig\` is NULL, use 3.

- ...:

  Additional optional arguments. None are used at present.

## Value

nls control object

## Author

Matthew L. Fidler

## Examples

``` r
# \donttest{

one.cmt <- function() {
  ini({
   tka <- 0.45
   tcl <- log(c(0, 2.7, 100))
   tv <- 3.45
   add.sd <- 0.7
 })
 model({
   ka <- exp(tka)
   cl <- exp(tcl)
   v <- exp(tv)
   linCmt() ~ add(add.sd)
 })
}

# Uses nlsLM from minpack.lm if available

fit1 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="nls", nlsControl(algorithm="LM"))
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of nls model...
#> ✔ done
#> → calculate jacobian
#> → calculate ∂(f)/∂(θ)
#> → finding duplicate expressions in nls gradient...
#> → optimizing duplicate expressions in nls gradient...
#> → finding duplicate expressions in nls pred-only...
#>  
#>  
#>  
#>  
#> → calculating covariance
#> ✔ done
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → finding duplicate expressions in EBE model...
#> → optimizing duplicate expressions in EBE model...
#> → compiling EBE model...
#>  
#>  
#> ✔ done
#> → Calculating residuals/tables
#> ✔ done
#> → compress origData in nlmixr2 object, save 6584
#> → compress parHistData in nlmixr2 object, save 2232

# Uses port and respect parameter boundaries
fit2 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="nls", nlsControl(algorithm="port"))
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of nls model...
#> ✔ done
#> → calculate jacobian
#> → calculate ∂(f)/∂(θ)
#> → finding duplicate expressions in nls gradient...
#> → optimizing duplicate expressions in nls gradient...
#> → finding duplicate expressions in nls pred-only...
#>  
#>  
#>  
#>  
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → finding duplicate expressions in EBE model...
#> → optimizing duplicate expressions in EBE model...
#> → compiling EBE model...
#>  
#>  
#> ✔ done
#> → Calculating residuals/tables
#> ✔ done
#> → compress origData in nlmixr2 object, save 6584
#> → compress parHistData in nlmixr2 object, save 2184

# You can access the underlying nls object with `$nls`
fit2$nls
#> Nonlinear regression model
#>   model: 0 ~ nlmixr2est::.nlmixrNlsFunValGrad(DV, tka, tcl, tv)
#>    data: nlmixr2est::.nlmixrNlsData()
#>     tka     tcl      tv 
#> -1.0097 -0.6696  1.0423 
#>  residual sum-of-squares: 249.7
#> 
#> Algorithm "port", convergence message: relative convergence (4)

# }
```
