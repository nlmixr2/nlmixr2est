# Control for n1qn1 estimation method in nlmixr2

Control for n1qn1 estimation method in nlmixr2

## Usage

``` r
n1qn1Control(
  epsilon = (.Machine$double.eps)^0.25,
  max_iterations = 10000,
  nsim = 10000,
  imp = 0,
  print.functions = FALSE,
  returnN1qn1 = FALSE,
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
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
  rxControl = NULL,
  optExpression = TRUE,
  sumProd = FALSE,
  literalFix = TRUE,
  literalFixRes = TRUE,
  addProp = c("combined2", "combined1"),
  calcTables = TRUE,
  compress = FALSE,
  covMethod = c("r", "n1qn1", ""),
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  ...
)
```

## Arguments

- epsilon:

  Precision of estimate for n1qn1 optimization.

- max_iterations:

  Number of iterations

- nsim:

  Number of function evaluations

- imp:

  Verbosity of messages.

- print.functions:

  Boolean to control if the function value and parameter estimates are
  echoed every time a function is called.

- returnN1qn1:

  return the n1qn1 output instead of the nlmixr2 fit

- stickyRecalcN:

  The number of bad ODE solves before reducing the atol/rtol for the
  rest of the problem.

- maxOdeRecalc:

  Maximum number of times to reduce the ODE tolerances and try to
  resolve the system if there was a bad ODE solve.

- odeRecalcFactor:

  The ODE recalculation factor when ODE solving goes bad, this is the
  factor the rtol/atol is reduced

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

- literalFixRes:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

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

  Ignored parameters

## Value

bobqya control structure

## Author

Matthew L. Fidler

## Examples

``` r
# \donttest{
# A logit regression example with emax model

dsn <- data.frame(i=1:1000)
dsn$time <- exp(rnorm(1000))
dsn$DV=rbinom(1000,1,exp(-1+dsn$time)/(1+exp(-1+dsn$time)))

mod <- function() {
 ini({
   E0 <- 0.5
   Em <- 0.5
   E50 <- 2
   g <- fix(2)
 })
 model({
   v <- E0+Em*time^g/(E50^g+time^g)
   ll(bin) ~ DV * v - log(1 + exp(v))
 })
}

fit2 <- nlmixr(mod, dsn, est="n1qn1")
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of population log-likelihood model...
#> ✔ done
#> → calculate jacobian
#> → calculate ∂(f)/∂(θ)
#> → finding duplicate expressions in nlm llik gradient...
#> → optimizing duplicate expressions in nlm llik gradient...
#> → finding duplicate expressions in nlm pred-only...
#> → optimizing duplicate expressions in nlm pred-only...
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

print(fit2)
#> ── nlmixr² log-likelihood n1qn1 ──
#> 
#>           OBJF     AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -716.4269 1127.45 1142.173      -560.7251        536.7183        68.24358
#> 
#> ── Time (sec $time): ──
#> 
#>            setup table compress    other
#> elapsed 0.002988 0.033    0.001 0.896012
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>        Est.     SE  %RSE   Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> E0  -0.6426 0.2388 37.16 -0.6426 (-1.111, -0.1746)                    
#> Em    6.071  2.854    47     6.071 (0.4785, 11.66)                    
#> E50   3.011  1.308 43.43     3.011 (0.4482, 5.575)                    
#> g         2  FIXED FIXED                         2                    
#>  
#>   Covariance Type ($covMethod): r
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED      v
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1     0.0394     0 -0.423 -0.642
#> 2 1     0.0557     1 -1.06  -0.641
#> 3 1     0.0604     1 -1.06  -0.640
#> # ℹ 997 more rows

# you can also get the nlm output with fit2$n1qn1

fit2$n1qn1
#> $value
#> [1] 560.7251
#> 
#> $par
#>         E0         Em        E50 
#> -0.6426484  6.0713089  3.0114656 
#> 
#> $H
#>              [,1]         [,2]         [,3]
#> [1,]  0.001772187  0.002754645 -0.007253828
#> [2,]  0.002754645  0.009784363 -0.021157083
#> [3,] -0.007253828 -0.021157083  0.050480308
#> 
#> $c.hess
#>  [1]  0.001772187  0.002754645 -0.007253828  0.009784363 -0.021157083
#>  [6]  0.050480308  0.000000000  0.000000000  0.000000000  0.000000000
#> [11]  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000
#> [16]  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000
#> [21]  0.000000000  0.000000000  0.000000000  0.000000000
#> 
#> $n.fn
#> [1] 40
#> 
#> $n.gr
#> [1] 40
#> 
#> $scaleC
#> [1] 0.003022222 0.040289311 0.035870774
#> 
#> $par.scaled
#>         E0         Em        E50 
#> -379.08228  137.28256   29.19749 
#> 
#> $hessian
#>               E0           Em          E50
#> E0   0.001774938  0.002754683 -0.007254274
#> Em   0.002754683  0.009784327 -0.021155968
#> E50 -0.007254274 -0.021155968  0.050439198
#> 
#> $covMethod
#> [1] "r"
#> 
#> $cov.scaled
#>           E0       Em      E50
#> E0  6244.783 1974.786 1726.433
#> Em  1974.786 5016.360 2388.055
#> E50 1726.433 2388.055 1329.237
#> 
#> $cov
#>             E0        Em       E50
#> E0  0.05703874 0.2404563 0.1871616
#> Em  0.24045630 8.1426982 3.4512380
#> E50 0.18716156 3.4512380 1.7103453
#> 
#> $r
#>                E0           Em          E50
#> E0   0.0008874692  0.001377342 -0.003627137
#> Em   0.0013773416  0.004892164 -0.010577984
#> E50 -0.0036271370 -0.010577984  0.025219599
#> 

# The nlm control has been modified slightly to include
# extra components and name the parameters
# }
```
