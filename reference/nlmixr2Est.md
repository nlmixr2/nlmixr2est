# Generic for nlmixr2 estimation methods

Generic for nlmixr2 estimation methods

## Usage

``` r
# S3 method for class 'agq'
nlmixr2Est(env, ...)

# S3 method for class 'bobyqa'
nlmixr2Est(env, ...)

# S3 method for class 'fo'
nlmixr2Est(env, ...)

# S3 method for class 'foce'
nlmixr2Est(env, ...)

# S3 method for class 'focei'
nlmixr2Est(env, ...)

# S3 method for class 'output'
nlmixr2Est(env, ...)

# S3 method for class 'foi'
nlmixr2Est(env, ...)

# S3 method for class 'laplace'
nlmixr2Est(env, ...)

# S3 method for class 'lbfgsb3c'
nlmixr2Est(env, ...)

# S3 method for class 'n1qn1'
nlmixr2Est(env, ...)

# S3 method for class 'newuoa'
nlmixr2Est(env, ...)

# S3 method for class 'nlm'
nlmixr2Est(env, ...)

# S3 method for class 'nlme'
nlmixr2Est(env, ...)

# S3 method for class 'nlminb'
nlmixr2Est(env, ...)

nlmixr2Est(env, ...)

# Default S3 method
nlmixr2Est(env, ...)

# S3 method for class 'nls'
nlmixr2Est(env, ...)

# S3 method for class 'optim'
nlmixr2Est(env, ...)

# S3 method for class 'posthoc'
nlmixr2Est(env, ...)

# S3 method for class 'rxSolve'
nlmixr2Est(env, ...)

# S3 method for class 'simulate'
nlmixr2Est(env, ...)

# S3 method for class 'simulation'
nlmixr2Est(env, ...)

# S3 method for class 'predict'
nlmixr2Est(env, ...)

# S3 method for class 'saem'
nlmixr2Est(env, ...)

# S3 method for class 'uobyqa'
nlmixr2Est(env, ...)
```

## Arguments

- env:

  Environment for the nlmixr2 estimation routines.

  This needs to have:

  \- rxode2 ui object in \`\$ui\`

  \- data to fit in the estimation routine in \`\$data\`

  \- control for the estimation routine's control options in \`\$ui\`

- ...:

  Other arguments provided to \`nlmixr2Est()\` provided for flexibility
  but not currently used inside nlmixr

## Value

nlmixr2 fit object

## Details

This is a S3 generic that allows others to use the nlmixr2 environment
to do their own estimation routines

## Author

Matthew Fidler
