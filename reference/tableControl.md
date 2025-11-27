# Output table/data.frame options

Output table/data.frame options

## Usage

``` r
tableControl(
  npde = NULL,
  cwres = NULL,
  nsim = 300,
  ties = TRUE,
  censMethod = c("truncated-normal", "cdf", "ipred", "pred", "epred", "omit"),
  seed = 1009,
  cholSEtol = (.Machine$double.eps)^(1/3),
  state = TRUE,
  lhs = TRUE,
  eta = TRUE,
  covariates = TRUE,
  addDosing = FALSE,
  subsetNonmem = TRUE,
  cores = NULL,
  keep = NULL,
  drop = NULL
)
```

## Arguments

- npde:

  When TRUE, request npde regardless of the algorithm used.

- cwres:

  When TRUE, request CWRES and FOCEi likelihood regardless of the
  algorithm used.

- nsim:

  represents the number of simulations. For rxode2, if you supply single
  subject event tables (created with `[eventTable()]`)

- ties:

  When \`TRUE\` jitter prediction-discrepancy points to discourage ties
  in cdf.

- censMethod:

  Handle censoring method:

  \- \`"truncated-normal"\` Simulates from a truncated normal
  distribution under the assumption of the model and censoring.

  \- \`"cdf"\` Use the cdf-method for censoring with npde and use this
  for any other residuals (\`cwres\` etc)

  \- \`"omit"\` omit the residuals for censoring

- seed:

  an object specifying if and how the random number generator should be
  initialized

- cholSEtol:

  The tolerance for the \`rxode2::choleSE\` function

- state:

  is a Boolean indicating if \`state\` values will be included (default
  \`TRUE\`)

- lhs:

  is a Boolean indicating if remaining \`lhs\` values will be included
  (default \`TRUE\`)

- eta:

  is a Boolean indicating if \`eta\` values will be included (default
  \`TRUE\`)

- covariates:

  is a Boolean indicating if covariates will be included (default
  \`TRUE\`)

- addDosing:

  Boolean indicating if the solve should add rxode2 EVID and related
  columns. This will also include dosing information and estimates at
  the doses. Be default, rxode2 only includes estimates at the
  observations. (default `FALSE`). When `addDosing` is `NULL`, only
  include `EVID=0` on solve and exclude any model-times or `EVID=2`. If
  `addDosing` is `NA` the classic `rxode2` EVID events are returned.
  When `addDosing` is `TRUE` add the event information in NONMEM-style
  format; If `subsetNonmem=FALSE` rxode2 will also include extra event
  types (`EVID`) for ending infusion and modeled times:

  - `EVID=-1` when the modeled rate infusions are turned off (matches
    `rate=-1`)

  - `EVID=-2` When the modeled duration infusions are turned off
    (matches `rate=-2`)

  - `EVID=-10` When the specified `rate` infusions are turned off
    (matches `rate>0`)

  - `EVID=-20` When the specified `dur` infusions are turned off
    (matches `dur>0`)

  - `EVID=101,102,103,...` Modeled time where 101 is the first model
    time, 102 is the second etc.

- subsetNonmem:

  subset to NONMEM compatible EVIDs only. By default `TRUE`.

- cores:

  Number of cores used in parallel ODE solving. This is equivalent to
  calling
  [`setRxThreads()`](https://nlmixr2.github.io/rxode2/reference/getRxThreads.html)

- keep:

  is the keep sent to the table

- drop:

  is the dropped variables sent to the table

## Value

A list of table options for nlmixr2

## Details

If you ever want to add CWRES/FOCEi objective function you can use the
[`addCwres`](https://nlmixr2.github.io/nlmixr2est/reference/addCwres.md)

If you ever want to add NPDE/EPRED columns you can use the
[`addNpde`](https://nlmixr2.github.io/nlmixr2est/reference/addNpde.md)

## Author

Matthew L. Fidler
