# Covariates explored by the VAE covariate search

Returns the candidate columns that \`nlmixr2(..., est = "vae")\` would
explore during automated covariate selection, using the same discovery
rules as the fit: every non-reserved data column that is constant within
each subject is a candidate. A numeric candidate with more than two
unique values is continuous and contributes one column per eligible
shape; anything else is categorical and contributes an indicator per
testable level. Columns sharing a \`group\` are alternate shapes of one
covariate, so at most one of them can enter a given parameter.
Time-varying columns cannot be searched and are excluded with a warning.

## Usage

``` r
vaeCovariates(
  data,
  warn = TRUE,
  shapes = c("power", "lin", "log", "identity", "center"),
  covCenterType = c("median", "mean"),
  covCenter = NULL,
  catCutoff = 0.05
)
```

## Arguments

- data:

  estimation dataset containing at least an \`ID\` column; column names
  are matched case-insensitively, as in the VAE fit

- warn:

  when \`TRUE\` (default) warn about time-varying columns excluded from
  the search; when \`FALSE\` exclude them silently

- shapes, covCenterType, covCenter, catCutoff:

  as in \[vaeControl()\]; control which shapes are explored and how
  covariates are centered

## Value

a data frame with one row per candidate search column and columns
\`covariate\` (the column name), \`raw\` (upper-cased data column it
comes from), \`shape\`, \`level\` (for categorical indicators),
\`group\` (mutual exclusion group), \`type\` and \`center\`; zero rows
when nothing qualifies

## Author

Matthew L. Fidler

## Examples

``` r
d <- data.frame(id = rep(1:3, each = 2), time = rep(0:1, 3), dv = rnorm(6),
                wt = rep(c(70, 80, 60), each = 2),
                sex = rep(c(0, 1, 0), each = 2))
vaeCovariates(d)
#>   covariate raw shape level group        type center
#> 1  WT_power  WT power  <NA>     1  continuous     70
#> 2    WT_lin  WT   lin  <NA>     1  continuous     70
#> 3       SEX SEX   cat  <NA>     2 categorical      0

# restrict the explored shapes
vaeCovariates(d, shapes = "power")
#>   covariate raw shape level group        type center
#> 1  WT_power  WT power  <NA>     1  continuous     70
#> 2       SEX SEX   cat  <NA>     2 categorical      0
```
