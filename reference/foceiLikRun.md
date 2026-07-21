# Evaluate a loaded general FOCE-family likelihood at supplied etas

Writes the population parameter vector into the loaded system and
returns the per-subject log-likelihood at the supplied etas, computed in
parallel over subjects. Requires a system loaded by \[foceiLikLoad()\].

## Usage

``` r
foceiLikRun(
  theta,
  eta,
  type = c("joint", "cond"),
  cores = rxode2::getRxThreads()
)
```

## Arguments

- theta:

  The estimation-scale parameter vector (length \`handle\$npars\`),
  matching the FOCEi optimizer parameterization: population thetas
  followed by the estimated Omega elements. \`handle\$initPar\` from
  \[foceiLikLoad()\] is a ready starting value.

- eta:

  A \`nid\` by \`neta\` matrix of random effects (one row per subject,
  in the loaded system's subject order).

- type:

  \`"joint"\` (default) returns the individual joint log density \`log
  p(y_i, eta_i)\`; \`"cond"\` returns the conditional data
  log-likelihood \`log p(y_i \| eta_i)\` alone. See Details.

- cores:

  Number of threads for the parallel per-subject evaluation.

## Value

A named numeric vector (length \`nid\`, named by subject id) of
per-subject log-likelihoods.

## Details

Both types are evaluated at the etas you supply, so both use each
subject's individual predictions; neither is a population (eta = 0)
quantity. They differ only by the random-effect prior term:

\- \`"cond"\` is the conditional data log-likelihood \`log p(y_i \|
eta_i)\`, the observation contribution alone. - \`"joint"\` is \`log
p(y_i, eta_i) = log p(y_i \| eta_i) + log p(eta_i)\`, which adds the
Gaussian random-effect prior \`log p(eta_i) = -0.5 eta_i' Omega^-1
eta_i + 0.5 log\|Omega^-1\| - neta/2 log(2 pi)\`.

So \`"joint"\` minus \`"cond"\` is exactly \`log p(eta_i)\`. \`"joint"\`
is the default because it is the usual target for MCMC/SAMBA-style
samplers: as a function of \`eta_i\` it is the individual's posterior
kernel, and it is the quantity the FOCEi inner problem optimizes over
the etas. Use \`"cond"\` when you supply the random-effect density
yourself, or when you need the observation contribution separately.

The prior is built from the loaded system's own \`Omega^-1\` and its log
determinant – the same Omega the inner likelihood uses – so \`"joint"\`
stays internally consistent with the engine rather than with the nominal
\`ini()\` values (the two differ by a small amount through Omega's
internal \`rxSymInv\` representation).

For Gaussian endpoints the observation contribution follows nlmixr2's
internal residual-likelihood convention, \`-0.5 err^2/r - 0.5 log(r)\`,
which omits the additive \`-0.5 log(2 pi)\` per observation; general
log-likelihood (\`ll()\`) endpoints contribute the user's log density as
written. The eta prior above is fully normalized. Both types are
therefore proper log densities up to a fixed per-observation constant
that does not depend on \`theta\` or \`eta\`, so likelihood ratios, and
any sampler that uses them, are unaffected.

## See also

\[foceiLikLoad()\], \[foceiLikUnload()\]

## Author

Matthew L. Fidler

## Examples

``` r

# \donttest{

one.cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

h <- foceiLikLoad(one.cmt, theo_sd, "focei")
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments

eta <- matrix(0, h$nid, h$neta)

# The individual joint log density log p(y_i, eta_i) (the default)
foceiLikRun(h$initPar, eta)
#>         1         2         3         4         5         6         7         8 
#> -94.15985 -62.28137 -63.71703 -63.26625 -85.38786 -36.32018 -43.59111 -48.51995 
#>         9        10        11        12 
#> -60.63507 -83.32305 -52.96870 -77.66747 

# The conditional data log-likelihood log p(y_i | eta_i) alone; the two
# differ by the Gaussian eta prior
foceiLikRun(h$initPar, eta, type = "cond")
#>         1         2         3         4         5         6         7         8 
#> -93.41180 -61.53332 -62.96899 -62.51820 -84.63981 -35.57213 -42.84306 -47.77190 
#>         9        10        11        12 
#> -59.88703 -82.57500 -52.22065 -76.91942 

# Non-zero etas
set.seed(42)
foceiLikRun(h$initPar, matrix(stats::rnorm(h$nid * h$neta, 0, 0.1), h$nid, h$neta))
#>         1         2         3         4         5         6         7         8 
#> -94.39492 -62.29327 -63.72068 -63.42463 -85.40293 -36.45668 -43.72089 -48.57601 
#>         9        10        11        12 
#> -60.72634 -83.39178 -52.99722 -77.87510 

# A new population parameter vector needs no reload
theta <- h$initPar
theta[1] <- theta[1] + 0.1
foceiLikRun(theta, eta)
#>         1         2         3         4         5         6         7         8 
#> -94.15971 -62.28123 -63.71688 -63.26615 -85.38773 -36.32011 -43.59105 -48.51985 
#>         9        10        11        12 
#> -60.63491 -83.32295 -52.96853 -77.66737 

foceiLikUnload()
# }
```
