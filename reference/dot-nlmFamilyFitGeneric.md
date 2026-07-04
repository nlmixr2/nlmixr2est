# Shared fit driver for the nlm-family estimation methods

Shared fit driver for the nlm-family estimation methods

## Usage

``` r
.nlmFamilyFitGeneric(
  env,
  method,
  fitModel,
  getTheta,
  controlToFocei,
  returnFlag,
  objective = NULL,
  message = function(fit) fit$message,
  emitFitWarnings = FALSE,
  extra = "",
  adjustOutput = TRUE,
  postSetup = NULL
)
```

## Arguments

- env:

  dispatch environment (provides \`ui\`, \`control\`, \`data\`,
  \`table\`)

- method:

  estimation-method string; also the slot the raw fit is stored under
  (e.g. \`"nlm"\` -\> \`.ret\[\["nlm"\]\]\`)

- fitModel:

  \`function(ui, dataSav)\` running the optimizer

- getTheta:

  \`function(fit, ui)\` returning the full theta vector

- controlToFocei:

  \`function(env)\` translating the control to a focei-style control for
  output assembly

- returnFlag:

  rxode2 control flag name that short-circuits and returns the raw
  optimizer result (e.g. \`"returnNlm"\`)

- objective:

  optional \`function(fit)\` returning the raw objective; when \`NULL\`
  the driver does not set \`\$objective\` (a \`postSetup\` closure did)

- message:

  \`function(fit)\` returning the \`\$message\` (default
  \`fit\$message\`)

- emitFitWarnings:

  when TRUE, re-emit the warnings collected from \`fitModel\` via
  \`warning()\` (nlm does this; the others do not)

- extra:

  \`\$extra\` print string, or a \`function(control)\` returning it

- adjustOutput:

  when TRUE, run \`.nlmFamilyAdjustOutput()\`

- postSetup:

  optional \`function(ret, ui, fitList)\` returning a modified \`ret\`,
  run right after the raw fit is stored and before
  \`.nlmFamilyAdjustOutput()\` (for methods that set
  cov/covMethod/objective with custom values)

## Value

the assembled nlmixr2 fit (or the raw optimizer result if
\`returnFlag\`)

## Author

Matthew L. Fidler
