# This adds a pre-processing hook to nlmixr2est

This pre-processing hook is run before the estimation process begins. It
is useful for modifying the user interface, the estimation object, the
data, or the control object before the estimation process begins. The
function must take four arguments: ui, est, data, and control. The
function must return a list with elements 'ui', 'est', 'data', and/or
'control'. If the element is not returned, the original object is used.
If the element is returned, the original object is replaced with the new
object.

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
[`preProcessHooks()`](https://nlmixr2.github.io/nlmixr2est/reference/preProcessHooks.md),
[`preProcessHooksRm()`](https://nlmixr2.github.io/nlmixr2est/reference/preProcessHooksRm.md)

## Author

Matthew L. Fidler
