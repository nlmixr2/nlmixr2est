# Load a general FOCE-family likelihood into memory

Compiles the inner (FOCEi sensitivity) model from an rxode2 UI model,
preprocesses the data, and sets up the FOCEi inner problem in memory so
that individual log-likelihoods can be evaluated repeatedly (in parallel
per subject) at supplied etas without recompiling – the setup used
internally by \`est="advi"\`, \`est="vae"\` and the f-SAEM fast kernel,
exposed here for MCMC/SAMBA-style callers (issue \#414).

## Usage

``` r
foceiLikLoad(
  object,
  data,
  likelihood = c("focei", "focep", "foce"),
  rxControl = rxode2::rxControl(),
  ...
)
```

## Arguments

- object:

  An \`rxode2\`/\`nlmixr2\` UI model (a model function or its compiled
  UI).

- data:

  The estimation data (a data frame with the usual nlmixr2 columns).

- likelihood:

  The individual likelihood type: \`"focei"\` (FOCE with interaction),
  \`"focep"\` (FOCE+, interaction off with the residual variance at the
  conditional eta) or \`"foce"\` (NONMEM-style FOCE, residual variance
  frozen at eta=0).

- rxControl:

  An \[rxode2::rxControl()\] object for the ODE solving options.

- ...:

  Additional solving/model options passed to \`.foceiLikControl\` (e.g.
  \`optExpression\`, \`addProp\`, \`eventSens\`).

## Value

Invisibly, a handle list with the loaded system's dimensions:
\`initPar\` (the estimation-scale parameter vector at the model's
initial estimates, a ready \`theta\` for \[foceiLikRun()\]), \`npars\`,
\`ntheta\`, \`neta\`, \`nid\`, \`thetaNames\`, \`etaNames\`, \`idLvl\`
and \`likelihood\`.

## Details

Only one likelihood system may be loaded at a time; loading errors if
one is already loaded. Use \[foceiLikRun()\] to evaluate and
\[foceiLikUnload()\] to free.

## See also

\[foceiLikRun()\], \[foceiLikUnload()\]

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

# Set the likelihood up in memory once; only one may be loaded at a time
h <- foceiLikLoad(one.cmt, theo_sd, "focei")
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments

# The handle carries the dimensions and a ready starting parameter vector
h$nid
#> [1] 12
h$neta
#> [1] 3
h$initPar
#> [1] 0.450000 1.000000 3.450000 0.700000 1.136219 1.351200 1.778279

# Individual joint log-likelihood at eta = 0, one value per subject
eta <- matrix(0, h$nid, h$neta)
foceiLikRun(h$initPar, eta)
#>         1         2         3         4         5         6         7         8 
#> -94.15985 -62.28137 -63.71703 -63.26625 -85.38786 -36.32018 -43.59111 -48.51995 
#>         9        10        11        12 
#> -60.63507 -83.32305 -52.96870 -77.66747 

# Free it when done (loading again before this errors)
foceiLikUnload()
# }
```
