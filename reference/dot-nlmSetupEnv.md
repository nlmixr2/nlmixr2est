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

  A list with \`predOnly\` (predictions-only model in terms of
  \`THETA\[#\]\`/\`DV\`), \`eventTheta\` (0/1 per THETA flagging
  event-related parameters that need Shi2021 finite differences, same
  length as \`par\`), and \`thetaGrad\` (needed when solveType != 1;
  gives value/gradient per THETA). See \`ui\$nlmSensModel\` or
  \`ui\$nlmRxModel\` for examples.

- control:

  control structure; required: \`rxControl\`, \`stickyRecalcN\`,
  \`maxOdeRecalc\`, \`odeRecalcFactor\`. Optional: \`solveType\`,
  \`eventType\`, \`shi21maxFD\`, \`shiErr\`, \`optimHessType\`,
  \`shi21maxHess\`, \`hessErr\`, \`useColor\`, \`printNcol\`, \`print\`,
  \`normType\`, \`scaleType\`, \`scaleCmin\`, \`scaleCmax\`,
  \`scaleTo\`, \`scaleC\`, \`gradTo\` (default 0 if missing).

- lower:

  lower bounds, will be scaled if present

- upper:

  upper bounds, will be scaled if present

## Value

nlm solve environment; key fields: \`\$par.ini\`, \`\$lower\`,
\`\$upper\` (all scaled), and \`\$.ctl\` (control structure).

## Details

No rxode2 solving should occur between setup calls; prints the solving
header if \`print != 0\`.

## Author

Matthew Fidler
