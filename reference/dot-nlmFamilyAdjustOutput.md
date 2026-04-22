# Adjust nlm and family output environment

Will take information like \`\$censInformation\`, \`\$parHistData\`,
\`\$cov\` and \`\$covMethod\` from the ret\[\[str\]\] and put it
directly in the environment \`ret\`

## Usage

``` r
.nlmFamilyAdjustOutput(ret, str)
```

## Arguments

- ret:

  environment for fit output that needs to be adjusted

- str:

  string for the fit output

## Value

updated environment

## Author

Matthew L. Fidler
