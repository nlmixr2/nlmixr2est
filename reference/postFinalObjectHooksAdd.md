# Register a post-final-object hook run after estimation completes

\`fun\` must take one argument (\`ret\`) and return the finalized return
object (or \`NULL\`/nothing to leave it unchanged).

## Usage

``` r
postFinalObjectHooksAdd(name, fun)
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
[`postFinalObjectHooksRm()`](https://nlmixr2.github.io/nlmixr2est/reference/postFinalObjectHooksRm.md),
[`preFinalParTableHooks()`](https://nlmixr2.github.io/nlmixr2est/reference/preFinalParTableHooks.md),
[`preFinalParTableHooksAdd()`](https://nlmixr2.github.io/nlmixr2est/reference/preFinalParTableHooksAdd.md),
[`preFinalParTableHooksRm()`](https://nlmixr2.github.io/nlmixr2est/reference/preFinalParTableHooksRm.md),
[`preProcessHooks()`](https://nlmixr2.github.io/nlmixr2est/reference/preProcessHooks.md),
[`preProcessHooksAdd()`](https://nlmixr2.github.io/nlmixr2est/reference/preProcessHooksAdd.md),
[`preProcessHooksRm()`](https://nlmixr2.github.io/nlmixr2est/reference/preProcessHooksRm.md)

## Author

Matthew L. Fidler
