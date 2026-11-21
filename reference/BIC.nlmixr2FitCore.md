# Bayesian information criterion for nlmixr2 fits

\`BIC(fit)\` charges every estimated parameter \`log(nobs)\`. The
\`type\` variants charge the parameters tied to random effects
\`log(nsub)\` instead (Delattre et al. 2014), which keeps BIC from
favouring over-parameterised random-effect structures.

## Usage

``` r
# S3 method for class 'nlmixr2FitCore'
BIC(object, ..., type = c("fixed", "mixed", "iiv", "random"))

# S3 method for class 'nlmixr2FitCoreSilent'
BIC(object, ..., type = c("fixed", "mixed", "iiv", "random"))
```

## Arguments

- object:

  nlmixr2 fit

- ...:

  other fitted models to compare, giving a data frame with columns
  \`df\` and \`BIC\`, one row per model; a character value here is taken
  as \`type\`, so \`BIC(fit, "mixed")\` works

- type:

  which parameters are charged \`log(nsub)\`:

  \- \`"fixed"\`: none; equal to \`BIC(fit)\`.

  \- \`"mixed"\`: omega elements and the thetas of eta-bearing
  structural parameters.

  \- \`"iiv"\`: omega (and kappa) elements only.

  \- \`"random"\`: every parameter tied to a random effect, which adds
  to \`"mixed"\` the residual-error thetas when an eta sits on the
  residual error.

## Value

the BIC value, or a data frame when several models are given

## Details

Parameters are classified from the fit's model:

\- Omega elements (including inter-occasion kappas) are always
random-effect parameters.

\- A theta is a structural random-effect parameter when it enters the
expression that defines a parameter carrying an eta (directly or through
intermediate assignments, but not through ODE states), as with \`tcl\`
and the covariate effect \`cwt\` in \`cl \<- exp(tcl + cwt \* log(WT /
70) + eta.cl)\`.

\- A residual-error theta is a random-effect parameter only when an eta
sits on the residual error, as with \`prop.sd\` in \`rsd \<- prop.sd \*
exp(eta.rsd); cp ~ prop(rsd)\`.

The linkage is read from the expression around each eta: the whole sum
it is added to (\`tcl + eta.cl\`), or the product it multiplies, where
the other factors are its typical value (\`tvcl \* exp(eta.cl)\`) unless
the eta's own factor already has one; then only covariate effects such
as \`(WT / 70)^cwt\` join it (in \`exp(tcl + eta.cl) / v\`, \`v\` stays
its own parameter). This is exact for mu-referenced parameters; for
other forms it is an approximation, so write each eta on its own
parameter line for an exact count. Fixed parameters are not counted. A
theta whose linkage cannot be determined is charged \`log(nobs)\`, as in
\`BIC(fit)\`.

## References

Delattre M, Lavielle M, Poursat MA (2014). A note on BIC in
mixed-effects models. Electronic Journal of Statistics, 8(1), 456-475.

## Author

Matthew L. Fidler

## Examples

``` r
# \donttest{
one.cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- log(2.7)
    tv <- 3.45
    eta.cl ~ 0.3
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv)
    linCmt() ~ add(add.sd)
  })
}
fit <- nlmixr2(one.cmt, theo_sd, est = "focei",
               control = foceiControl(print = 0))
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
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
#> rxode2 5.1.8.9000 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
#> 
#> Attaching package: ‘rxode2’
#> The following objects are masked from ‘package:nlmixr2est’:
#> 
#>     boxCox, yeoJohnson
#> → Calculating residuals/tables
#> ✔ done
BIC(fit)
#> [1] 460.6031
BIC(fit, type = "mixed")
#> [1] 455.8073
# }
```
