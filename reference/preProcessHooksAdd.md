# Register a pre-processing hook run before estimation begins

\`fun\` must take four arguments (\`ui\`, \`est\`, \`data\`,
\`control\`) and return a list with any of \`'ui'\`, \`'est'\`,
\`'data'\`, \`'control'\` to override; a non-returned element is left
unchanged.

## Usage

``` r
preProcessHooksAdd(name, fun)
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
[`preFinalParTableHooksAdd()`](https://nlmixr2.github.io/nlmixr2est/reference/preFinalParTableHooksAdd.md),
[`preFinalParTableHooksRm()`](https://nlmixr2.github.io/nlmixr2est/reference/preFinalParTableHooksRm.md),
[`preProcessHooks()`](https://nlmixr2.github.io/nlmixr2est/reference/preProcessHooks.md),
[`preProcessHooksRm()`](https://nlmixr2.github.io/nlmixr2est/reference/preProcessHooksRm.md)

## Author

Matthew L. Fidler
