# Load a general FOCE-family likelihood into memory

Compiles the inner (FOCEi sensitivity) model from an rxode2 UI model,
preprocesses the data, and sets up the FOCEi inner problem in memory so
that individual log-likelihoods can be evaluated repeatedly (in parallel
per subject) at supplied etas without recompiling – the setup used
internally by \`est="emvi"\`/\`est="fbvi"\`, \`est="vae"\` and the
f-SAEM fast kernel, exposed here for MCMC/SAMBA-style callers (issue
\#414).

## Usage

``` r
foceiLikLoad(
  object,
  data,
  likelihood = c("focei", "focep", "foce"),
  rxControl = rxode2::rxControl(),
  scale = c("focei", "natural"),
  thetaSens = FALSE,
  est = "focei",
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

- scale:

  The parameter scale \[foceiLikRun()\]'s \`theta\` is on: \`"focei"\`
  (default) is the FOCEi estimation scale (the historical behavior);
  \`"natural"\` pins the scaling to the identity (\`scaleType="mult"\`,
  \`scaleTo=0\`), so the population thetas in \`theta\` (and in the
  handle's \`initPar\`) are the natural-scale values directly comparable
  with \`ui\$iniDf\$est\` – no unscaling step for external callers
  (#939). With either choice the trailing omega-block entries of the
  parameter vector remain in the internal \`diagXform\` parameterization
  of \`chol(Omega^-1)\`; \`"natural"\` leaves them unscaled but does not
  change that parameterization.

- thetaSens:

  When \`TRUE\`, also build and wire the theta-sensitivity model
  (\`d(f)/d(theta)\`, \`d(V)/d(theta)\` forward sensitivities for the
  estimated non-mu structural and residual-error thetas), the model the
  imp/advi engines use for their outer population gradient. Off by
  default: building it costs compile and solve time a value-only caller
  should not pay (#939). The handle's \`thetaSensIdx\` reports which
  \`ntheta\` indices carry sensitivities (\`integer(0)\` when none do,
  e.g. when every theta is mu-referenced).

- est:

  Estimation-method name the standard pre-process hooks run against
  (default \`"focei"\`). The hooks read the named method's capability
  attributes (for example whether a bounded-transform reparameterization
  is wanted), so a caller preparing the problem for a different consumer
  can name it here and get that consumer's preprocessing guarantees
  instead of focei's (#939). The inner engine itself is unaffected.

- ...:

  Additional solving/model options passed to \`.foceiLikControl\` (e.g.
  \`optExpression\`, \`addProp\`, \`eventSens\`).

## Value

Invisibly, a handle list with the loaded system's dimensions:
\`initPar\` (the parameter vector at the model's initial estimates on
the requested \`scale\`, a ready \`theta\` for \[foceiLikRun()\]),
\`npars\`, \`ntheta\`, \`neta\`, \`nid\`, \`thetaNames\`, \`etaNames\`,
\`idLvl\`, \`likelihood\`, \`scale\`, \`thetaSens\` and
\`thetaSensIdx\`.

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
#>           1           2           3           4           5           6 
#> -24.4017587  -9.3365473  -0.6718843  -5.5542461 -19.0895200 -33.8105271 
#>           7           8           9          10          11          12 
#> -24.9543558  -8.8147973 -32.7217521 -30.0008367 -18.5848785 -13.1123082 

# Free it when done (loading again before this errors)
foceiLikUnload()
# }
```
