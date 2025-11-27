# Handle Model Object

Handle Model Object

## Usage

``` r
nmObjHandleModelObject(model, env)

# S3 method for class 'saemModelList'
nmObjHandleModelObject(model, env)

# S3 method for class 'foceiModelList'
nmObjHandleModelObject(model, env)

# Default S3 method
nmObjHandleModelObject(model, env)
```

## Arguments

- model:

  model list should have at least:

  \- \`predOnly\` – this is the prediction model with all the left
  handed equations added so they will be added the table. The model
  should have \`rx_pred\_\`, the model based prediction, as the first
  defined lhs component. The second component should be \`rx_r\_\`, the
  variance of the prediction. These variables may change based on
  distribution type. In additional all interesting calculated variables
  should be included.

  \- \`predNoLhs\` – This is the prediction model. It only has the
  prediction and no left handed equations.

- env:

  Environment for the fit information

## Value

This returns the \`\$model\` object for a fit. It is a s3 method because
it may be different between different model types
