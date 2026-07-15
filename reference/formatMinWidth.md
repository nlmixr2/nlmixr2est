# Format numeric values to minimize the printing width

Special values (\`NaN\`, \`Inf\`, \`-Inf\`, \`0\`) are shown as-is;
\`NA\` uses \`naValue\`.

## Usage

``` r
formatMinWidth(x, digits = 3, naValue = "NA")
```

## Arguments

- x:

  The numeric vector to convert

- digits:

  The number of significant digits to show

- naValue:

  The value to return if \`is.na(x)\`

## Value

A character vector converting the numbers with minimum width
