# Likelihood penalty specification

Describes one compiled likelihood penalty (see
\`vignette("likelihood-penalties", package = "nlmixr2est")\`) to apply
during a single fit. A model extension attaches specs to a model with
\[nlmixr2LikPenaltyAdd()\]; they can also be passed to \`nlmixr2()\` as
\`likPenalty=\`.

## Usage

``` r
nlmixr2LikPenaltySpec(
  ptr,
  thetas = character(0),
  par = numeric(0),
  data = NULL,
  inObjf = NA,
  name = NULL
)
```

## Arguments

- ptr:

  the penalty: an external pointer to a C \`nlmixrLikPenalty\`, a
  function returning one (resolved when the fit starts, so the spec
  survives \`saveRDS()\`), or the name of a built-in penalty
  (\`"thetaL2"\`, \`"etaL2"\`, \`"omegaL2"\`).

- thetas:

  names of the population parameters the penalty acts on; handed to the
  C callbacks as 0-based theta indices (\`ctx-\>sel\`).

- par:

  numeric settings handed to the C callbacks (\`ctx-\>par\`).

- data:

  optional external pointer handed to the C callbacks (\`ctx-\>data\`);
  defaults to the bundle's own \`data\`.

- inObjf:

  should the penalty be counted in the reported objective function?
  \`NA\` uses the bundle's \`NLMIXR_PEN_IN_OBJF\` flag.

- name:

  label used in \`fit\$likPenalty\`; defaults to the bundle's name.

## Value

a \`nlmixr2LikPenaltySpec\` object

## Author

Matthew L. Fidler

## Examples

``` r
# ridge penalty of 10*tka^2 (on the -2LL scale), shaping the fit only
nlmixr2LikPenaltySpec("thetaL2", thetas = "tka", par = 10)
#> likelihood penalty thetaL2 on tka par: 10 inObjf: NA 
```
