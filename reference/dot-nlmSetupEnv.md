# Setup a nonlinear system for optimization

Setup a nonlinear system for optimization

## Usage

``` r
.nlmSetupEnv(par, ui, data, modelInfo, control, lower = NULL, upper = NULL)
```

## Arguments

- par:

  A named vector of initial estimates to setup the nonlinear model
  solving environment. The names of the parameter should match the names
  of the model to run (not \`THETA\[#\]\` as required in the
  \`modelInfo\` argument)

- ui:

  rxode2 ui model

- data:

  rxode2 compatible data for solving/setting up

- modelInfo:

  A list containing the following elements:

  \- \`predOnly\` – A model with only predictions calculated. These
  predictions should be in terms of \`THETA\[#\]\` and \`DV\`. The

  \- \`eventTheta\` is an indicator if the \`THETA\[#\]\` is related to
  an event (like \`dur(x)\` \`f(x)\`). These variables will use Shi2021
  finite differences and need to be indicated when setting up the
  solving environment. When finite differences are required, this is
  \`1L\` when they are not it should be \`0L\`. This should match the
  length of \`par\`

  \- \`thetaGrad\` – needed when solveType != 1; a model that gives the
  value and gradient of each \`THETA\[#\]\`

  An example can be found with \`ui\$nlmSensModel\` or
  \`ui\$nlmRxModel\`

- control:

  is a control structure with a few required elements:

  \- \`rxControl\` represents the rxode2 solving options - \`solveType\`
  integer indicating the solveType (optional) - \`stickyRecalcN\` -
  \`maxOdeRecalc\` - \`odeRecalcFactor\` - \`eventType\` (optional) -
  \`shi21maxFD\` (optional) - \`shiErr\` (optional) - \`optimHessType\`
  (optional) - \`shi21maxHess\` (optional) - \`hessErr\` (optional) -
  \`useColor\` - \`printNcol\` - \`print\` - \`normType\` -
  \`scaleType\` - \`scaleCmin\` - \`scaleCmax\` - \`scaleTo\` -
  \`scaleC\` - \`gradTo\` (optional); if missing assumed gradTo=0

- lower:

  lower bounds, will be scaled if present

- upper:

  upper bounds, will be scaled if present

## Value

nlm solve environment; of interest

\`\$par.ini\` – scaled parameter initial value

\`\$lower\` – scaled parameter lower value

\`\$upper\` – scaled parameter upper value

\`\$.ctl\` – control structure

## Details

In between using this, rxode2 solving should not be called.

This will also print the header for solving (if print != 0)

## Author

Matthew Fidler
