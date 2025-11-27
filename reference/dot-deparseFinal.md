# Deparse finalize a control or related object into a language object

This function deparses an object into a language expression, optionally
using a custom function for specific elements.

## Usage

``` r
.deparseFinal(default, object, w, var, fun = NULL)
```

## Arguments

- default:

  A default object used for comparison; This is the estimation control
  procedure. It should have a class matching the function that created
  it.

- object:

  The object to be deparsed into a language exression

- w:

  A vector of indices indicating which elements are different and need
  to be deparsed. This likely comes from \`.deparseDifferent()\`

- var:

  A string representing the variable name to be assigned in the deparsed
  expression.

- fun:

  An optional custom function to handle specific elements during
  deparsing. Default is NULL. This handles things that are specific to
  an estimation control and is used by functions like
  \`rxUiDeparse.saemControl()\`

## Value

A language object representing the deparsed expression.

## Author

Matthew L. Fidler
