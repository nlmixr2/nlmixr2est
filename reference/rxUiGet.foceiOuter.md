# Build the augmented outer-gradient sensitivity model (compiled model + \`dirs\` + \`P2\`) for a UI. This is the persistent \`..outer\` sibling of the inner model: it depends only on the model + direction set (NOT theta/eta/omega), so it is built once during model setup (via \`rxUiGet.foceiModel\`/\`foceModel\`, which qs2-cache the whole model list) and reused across every outer-gradient call. Callable independently as \`ui\$foceiOuter\`. \`NULL\` when out of analytic scope (the gradient then falls back to finite differences).

Build the augmented outer-gradient sensitivity model (compiled model +
\`dirs\` + \`P2\`) for a UI. This is the persistent \`..outer\` sibling
of the inner model: it depends only on the model + direction set (NOT
theta/eta/omega), so it is built once during model setup (via
\`rxUiGet.foceiModel\`/\`foceModel\`, which qs2-cache the whole model
list) and reused across every outer-gradient call. Callable
independently as \`ui\$foceiOuter\`. \`NULL\` when out of analytic scope
(the gradient then falls back to finite differences).

## Usage

``` r
# S3 method for class 'foceiOuter'
rxUiGet(x, ...)
```
