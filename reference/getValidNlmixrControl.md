# Get valid nlmixr control object

Get valid nlmixr control object

## Usage

``` r
# S3 method for class 'advi'
getValidNlmixrCtl(control)

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

# S3 method for class 'mfoceif'
getValidNlmixrCtl(control)

# S3 method for class 'mfocef'
getValidNlmixrCtl(control)

# S3 method for class 'mfocepf'
getValidNlmixrCtl(control)

# S3 method for class 'ifoceif'
getValidNlmixrCtl(control)

# S3 method for class 'ifocef'
getValidNlmixrCtl(control)

# S3 method for class 'ifocepf'
getValidNlmixrCtl(control)

# S3 method for class 'focep'
getValidNlmixrCtl(control)

# S3 method for class 'foi'
getValidNlmixrCtl(control)

# S3 method for class 'fsaem'
getValidNlmixrCtl(control)

# S3 method for class 'ifocep'
getValidNlmixrCtl(control)

# S3 method for class 'imp'
getValidNlmixrCtl(control)

# S3 method for class 'impmap'
getValidNlmixrCtl(control)

# S3 method for class 'laplace'
getValidNlmixrCtl(control)

# S3 method for class 'lbfgsb3c'
getValidNlmixrCtl(control)

# S3 method for class 'mfocep'
getValidNlmixrCtl(control)

# S3 method for class 'mfocei'
getValidNlmixrCtl(control)

# S3 method for class 'ifocei'
getValidNlmixrCtl(control)

# S3 method for class 'mfoce'
getValidNlmixrCtl(control)

# S3 method for class 'ifoce'
getValidNlmixrCtl(control)

# S3 method for class 'magq'
getValidNlmixrCtl(control)

# S3 method for class 'iagq'
getValidNlmixrCtl(control)

# S3 method for class 'mlaplace'
getValidNlmixrCtl(control)

# S3 method for class 'ilaplace'
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

# S3 method for class 'mnpag'
getValidNlmixrCtl(control)

# S3 method for class 'inpag'
getValidNlmixrCtl(control)

# S3 method for class 'mnpb'
getValidNlmixrCtl(control)

# S3 method for class 'inpb'
getValidNlmixrCtl(control)

# S3 method for class 'npag'
getValidNlmixrCtl(control)

# S3 method for class 'npb'
getValidNlmixrCtl(control)

# S3 method for class 'optim'
getValidNlmixrCtl(control)

# S3 method for class 'posthoc'
getValidNlmixrCtl(control)

# S3 method for class 'qrpem'
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
