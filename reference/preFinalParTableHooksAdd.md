# Register a pre-final parameter table hook run before final table assembly

\`fun\` must take one argument (\`env\`) and update it in place (e.g.
\`cov\`, \`theta\`, \`thetaNames\`, \`thetaDf\`) before the final tables
are built.

## Usage

``` r
preFinalParTableHooksAdd(name, fun)
```

## Arguments

- name:

  Character vector representing the name of the hook

- fun:

  The function to run

## Value

The function that was added (invisibly)

## See also

Other preProcessHooks:
[`postFinalObjectHooks()`](https://nlmixr2.github.io/nlmixr2est/reference/postFinalObjectHooks.md),
[`postFinalObjectHooksAdd()`](https://nlmixr2.github.io/nlmixr2est/reference/postFinalObjectHooksAdd.md),
[`postFinalObjectHooksRm()`](https://nlmixr2.github.io/nlmixr2est/reference/postFinalObjectHooksRm.md),
[`preFinalParTableHooks()`](https://nlmixr2.github.io/nlmixr2est/reference/preFinalParTableHooks.md),
[`preFinalParTableHooksRm()`](https://nlmixr2.github.io/nlmixr2est/reference/preFinalParTableHooksRm.md),
[`preProcessHooks()`](https://nlmixr2.github.io/nlmixr2est/reference/preProcessHooks.md),
[`preProcessHooksAdd()`](https://nlmixr2.github.io/nlmixr2est/reference/preProcessHooksAdd.md),
[`preProcessHooksRm()`](https://nlmixr2.github.io/nlmixr2est/reference/preProcessHooksRm.md)

## Author

Matthew L. Fidler
