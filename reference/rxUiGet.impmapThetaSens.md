# Build the symengine env carrying the impmap sensitivity model (`..thetaSens`). For each estimated non-mu theta j it outputs rx\_\_sens_rx_pred\_\_BY_THETA_j\_\_\_ = d(f)/d(theta_j) and rx\_\_sens_rx_r\_\_BY_THETA_j\_\_\_ = d(V)/d(theta_j).

Build the symengine env carrying the impmap sensitivity model
(`..thetaSens`). For each estimated non-mu theta j it outputs
rx\_\_sens_rx_pred\_\_BY_THETA_j\_\_\_ = d(f)/d(theta_j) and
rx\_\_sens_rx_r\_\_BY_THETA_j\_\_\_ = d(V)/d(theta_j).

## Usage

``` r
# S3 method for class 'impmapThetaSens'
rxUiGet(x, ...)
```
