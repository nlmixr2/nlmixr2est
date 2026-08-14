# Register / remove an estimation interceptor

An interceptor is a \`function(env)\` consulted in \[nlmixr2Est()\]
before the standard estimation method is dispatched. It receives the
fully set-up estimation environment (\`env\$ui\`, \`env\$data\`,
\`env\$control\`, the est name in \`class(env)\[1\]\`, and any extra
\`nlmixr2()\` arguments in \`env\$.nlmixr2Dots\`). If it returns a
non-\`NULL\` value that value becomes the fit – the interceptor has
\*claimed\* the estimation; if it returns \`NULL\` it declines and the
next interceptor (or the standard method) runs. While a claimed
interceptor runs, interceptors are suppressed, so it may call
\`nlmixr2()\` internally and get the ordinary methods.

## Usage

``` r
registerEstInterceptor(name, fun)

removeEstInterceptor(name)
```

## Arguments

- name:

  unique interceptor name.

- fun:

  \`function(env)\` returning a finalized fit to claim, or \`NULL\`.

## Value

invisibly the function (register) or \`TRUE\`/\`FALSE\` (remove).

## Author

Matthew L. Fidler
