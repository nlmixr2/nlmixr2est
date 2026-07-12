# Validate the RPEM control (est="rpem")

Validate the RPEM control (est="rpem")

Get valid nlmixr control object

## Usage

``` r
# S3 method for class 'agq'
getValidNlmixrCtl(control)

# S3 method for class 'bobyqa'
getValidNlmixrCtl(control)

# S3 method for class 'fo'
getValidNlmixrCtl(control)

# S3 method for class 'foce'
getValidNlmixrCtl(control)

# S3 method for class 'foceif'
getValidNlmixrCtl(control)

# S3 method for class 'focef'
getValidNlmixrCtl(control)

# S3 method for class 'focepf'
getValidNlmixrCtl(control)

# S3 method for class 'mufoceif'
getValidNlmixrCtl(control)

# S3 method for class 'mufocef'
getValidNlmixrCtl(control)

# S3 method for class 'mufocepf'
getValidNlmixrCtl(control)

# S3 method for class 'irlsfoceif'
getValidNlmixrCtl(control)

# S3 method for class 'irlsfocef'
getValidNlmixrCtl(control)

# S3 method for class 'irlsfocepf'
getValidNlmixrCtl(control)

# S3 method for class 'focep'
getValidNlmixrCtl(control)

# S3 method for class 'foi'
getValidNlmixrCtl(control)

# S3 method for class 'fsaem'
getValidNlmixrCtl(control)

# S3 method for class 'imp'
getValidNlmixrCtl(control)

# S3 method for class 'impmap'
getValidNlmixrCtl(control)

# S3 method for class 'irlsfocep'
getValidNlmixrCtl(control)

# S3 method for class 'laplace'
getValidNlmixrCtl(control)

# S3 method for class 'lbfgsb3c'
getValidNlmixrCtl(control)

# S3 method for class 'mufocei'
getValidNlmixrCtl(control)

# S3 method for class 'irlsfocei'
getValidNlmixrCtl(control)

# S3 method for class 'mufoce'
getValidNlmixrCtl(control)

# S3 method for class 'irlsfoce'
getValidNlmixrCtl(control)

# S3 method for class 'muagq'
getValidNlmixrCtl(control)

# S3 method for class 'irlsagq'
getValidNlmixrCtl(control)

# S3 method for class 'mulaplace'
getValidNlmixrCtl(control)

# S3 method for class 'irlslaplace'
getValidNlmixrCtl(control)

# S3 method for class 'mufocep'
getValidNlmixrCtl(control)

# S3 method for class 'n1qn1'
getValidNlmixrCtl(control)

# S3 method for class 'newuoa'
getValidNlmixrCtl(control)

# S3 method for class 'nlm'
getValidNlmixrCtl(control)

# S3 method for class 'nlminb'
getValidNlmixrCtl(control)

# S3 method for class 'nls'
getValidNlmixrCtl(control)

# S3 method for class 'optim'
getValidNlmixrCtl(control)

# S3 method for class 'posthoc'
getValidNlmixrCtl(control)

# S3 method for class 'rpem'
getValidNlmixrCtl(control)

getValidNlmixrControl(control, est)

getValidNlmixrCtl(control)

# S3 method for class 'focei'
getValidNlmixrCtl(control)

# S3 method for class 'nlme'
getValidNlmixrCtl(control)

# S3 method for class 'saem'
getValidNlmixrCtl(control)

# S3 method for class 'rxSolve'
getValidNlmixrCtl(control)

# S3 method for class 'simulate'
getValidNlmixrCtl(control)

# S3 method for class 'simulation'
getValidNlmixrCtl(control)

# S3 method for class 'predict'
getValidNlmixrCtl(control)

# S3 method for class 'tableControl'
getValidNlmixrCtl(control)

# Default S3 method
getValidNlmixrCtl(control)

# S3 method for class 'uobyqa'
getValidNlmixrCtl(control)

# S3 method for class 'vae'
getValidNlmixrCtl(control)
```

## Arguments

- control:

  nlmixr control object

- est:

  Estimation routine

## Value

Valid control object based on estimation method run.

## Details

This is based on running the S3 method \`getValidNlmixrCtl()\` the
\`control\` object is put into a list and the class of this new list is
\`c(est, "getValidNlmixrControl")\`
