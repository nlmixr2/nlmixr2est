# Covariates explored by the VAE covariate search

Returns the subject-level covariates that \`nlmixr2(..., est = "vae")\`
would explore during automated covariate selection, using the same
discovery rules as the fit: every non-reserved numeric data column that
is constant within each subject is a candidate; a candidate with more
than two unique values (all positive) is treated as continuous (encoded
\`log(value/mean)\`), anything else as categorical (mean-centered).
Time-varying numeric columns cannot be searched and are excluded with a
warning.

## Usage

``` r
vaeCovariates(data, warn = TRUE)
```

## Arguments

- data:

  estimation dataset containing at least an \`ID\` column; column names
  are matched case-insensitively, as in the VAE fit

- warn:

  when \`TRUE\` (default) warn about time-varying numeric columns
  excluded from the search; when \`FALSE\` exclude them silently

## Value

a data frame with one row per explored covariate and columns
\`covariate\` (upper-cased column name), \`type\` (\`"continuous"\` or
\`"categorical"\`) and \`center\` (the population value the covariate is
centered at); zero rows when no covariates qualify

## Author

Matthew L. Fidler

## Examples

``` r
d <- data.frame(id = rep(1:3, each = 2), time = rep(0:1, 3), dv = rnorm(6),
                wt = rep(c(70, 80, 60), each = 2),
                sex = rep(c(0, 1, 0), each = 2))
vaeCovariates(d)
#>   covariate        type     center
#> 1        WT  continuous 70.0000000
#> 2       SEX categorical  0.3333333
```
