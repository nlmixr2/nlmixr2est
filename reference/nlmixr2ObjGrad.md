# Objective function and analytic gradient of a FOCEi-family fit at new parameters

Evaluates the fit's own objective function and its analytic outer
gradient (the \`foceiControl(fast=TRUE)\` C++ gradient) at user-supplied
population parameters, without running an outer optimization. The fit's
estimation method (\`"focei"\`, \`"foce"\`, \`"laplace"\`, \`"agq"\`,
...) is reused, since it determines the objective being differentiated.

## Usage

``` r
nlmixr2ObjGrad(
  fit,
  theta = NULL,
  omega = NULL,
  scale = c("untransformed", "transformed"),
  individual = FALSE
)
```

## Arguments

- fit:

  A fit from a FOCEi-family estimation method.

- theta:

  Population fixed effects and residual parameters to evaluate at, on
  the scale they are written in \`ini()\`. Either named (a subset of
  \`names(fit\$theta)\`; the rest stay at the fit's values) or unnamed
  and as long as \`fit\$theta\`. \`NULL\` uses \`fit\$theta\`.

- omega:

  Between-subject covariance matrix to evaluate at, with the same
  dimensions as \`fit\$omega\`. \`NULL\` uses \`fit\$omega\`.

- scale:

  The scale of the returned gradient for the omega elements.
  \`"untransformed"\` differentiates with respect to the variances and
  covariances themselves, named as in \`ini()\` (\`eta.cl\`,
  \`(eta.cl,eta.v)\`). \`"transformed"\` differentiates with respect to
  the Cholesky-based estimation-scale omega parameters the optimizer
  uses (\`om.chol.\<k\>\`, with \`foceiControl(diagXform=)\` applied to
  the diagonal). Thetas are always on their \`ini()\` scale.

- individual:

  Also return each subject's contribution to the objective and gradient.

## Value

A list with \`objective\` (the fit's objective at the point),
\`gradient\` (the derivative of \`objective\` for every estimated
parameter, named), \`par\` (the parameter values on the same scale as
\`gradient\`) and \`scale\`. With \`individual = TRUE\` it also has
\`individual\`, a data frame of each subject's \`ID\`, objective
contribution (\`OBJI\`) and gradient (one column per parameter), and
\`population\`, the part of the objective and gradient that belongs to
no subject (such as an \`ini()\` prior), so the subjects plus
\`population\` sum to the totals. The attribute \`nAnalyticGradDirect\`
counts the analytic gradient evaluations behind the result.

## Details

At the fit's own estimates the fit's empirical Bayes estimates are used
as-is, so the objective equals \`fit\$objective\`. At any other point
the etas are re-optimized (warm-started from the fit's) before the
objective and gradient are evaluated.

## Examples

``` r
# \donttest{
one.cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}
fit <- nlmixr2(one.cmt, theo_sd, "focei", control = foceiControl(print = 0))
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → calculate sensitivities
#> → calculate ∂(f)/∂(η)
#> → calculate ∂(R²)/∂(η)
#> → finding duplicate expressions in inner model...
#> → finding duplicate expressions in EBE model...
#> → compiling inner model...
#>  
#>  
#> ✔ done
#> → finding duplicate expressions in FD model...
#> → compiling EBE model...
#>  
#>  
#> ✔ done
#> → compiling events FD model...
#>  
#>  
#> ✔ done
#> → Calculating residuals/tables
#> ✔ done
nlmixr2ObjGrad(fit)
#> $objective
#> [1] 130.2909
#> 
#> $gradient
#>         tka         tcl          tv      add.sd      eta.ka      eta.cl 
#>  0.09317026  0.11469669  0.22822036  0.18756629 -0.17602287 -1.25862230 
#> 
#> $par
#>       tka       tcl        tv    add.sd    eta.ka    eta.cl 
#> 0.4034116 1.0280965 3.4295877 0.7823108 0.3271547 0.1192014 
#> 
#> $scale
#> [1] "untransformed"
#> 
#> attr(,"nAnalyticGradDirect")
#> [1] 1
nlmixr2ObjGrad(fit, theta = c(tka = 0.5))
#> $objective
#> [1] 130.621
#> 
#> $gradient
#>         tka         tcl          tv      add.sd      eta.ka      eta.cl 
#>   6.7356585   0.1335134 -11.0762779   0.3176872  -1.3964653  -1.2753121 
#> 
#> $par
#>       tka       tcl        tv    add.sd    eta.ka    eta.cl 
#> 0.5000000 1.0280965 3.4295877 0.7823108 0.3271547 0.1192014 
#> 
#> $scale
#> [1] "untransformed"
#> 
#> attr(,"nAnalyticGradDirect")
#> [1] 1
# }
```
