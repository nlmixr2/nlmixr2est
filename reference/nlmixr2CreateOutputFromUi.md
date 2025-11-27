# Create nlmixr output from the UI

Create nlmixr output from the UI

## Usage

``` r
nlmixr2CreateOutputFromUi(
  ui,
  data = NULL,
  control = NULL,
  table = NULL,
  env = NULL,
  est = "none"
)
```

## Arguments

- ui:

  This is the UI that will be used for the translation

- data:

  This has the data

- control:

  focei control for data creation

- table:

  Table options

- env:

  Environment setup which needs the following: - \`\$table\` for table
  options - \`\$origData\` – Original Data - \`\$dataSav\` – Processed
  data from .foceiPreProcessData - \`\$idLvl\` – Level information for
  ID factor added - \`\$covLvl\` – Level information for items to
  convert to factor - \`\$ui\` for ui object - \`\$fullTheta\` Full
  theta information - \`\$etaObf\` data frame with ID, etas and OBJI -
  \`\$cov\` For covariance - \`\$covMethod\` for the method of
  calculating the covariance - \`\$adjObf\` Should the objective
  function value be adjusted - \`\$objective\` objective function
  value - \`\$extra\` Extra print information - \`\$method\` Estimation
  method (for printing) - \`\$omega\` Omega matrix - \`\$theta\` Is a
  theta data frame - \`\$model\` a list of model information for table
  generation. Needs a \`predOnly\` model - \`\$message\` Message for
  display - \`\$est\` estimation method - \`\$ofvType\` (optional) tells
  the type of ofv is currently being use

  There are some more details that need to be described here

- est:

  Estimation method

## Value

nlmixr fit object

## Author

Matthew L. Fidler
