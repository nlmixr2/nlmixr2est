# Tagged list of the available nlmixr2 estimation methods

Returns the built-in (and any attribute-tagged third-party) \`est=\`
methods grouped by their estimation category, as used when an
unsupported method is requested.

## Usage

``` r
nlmixr2AllEstType()
```

## Value

data.frame with columns \`est\`, \`type\` and \`description\`

## Examples

``` r
nlmixr2AllEstType()
#>         est                   type                         description
#> 1        fo             Linearized                         First-Order
#> 2      foce             Linearized  First-Order Conditional Estimation
#> 3     focei             Linearized               FOCE with Interaction
#> 4     focep             Linearized FOCE+ (residual at conditional eta)
#> 5       foi             Linearized        First-Order with Interaction
#> 6      nlme             Linearized  Lindstrom-Bates alternating (nlme)
#> 7       agq Integral approximation        Adaptive Gaussian Quadrature
#> 8       imp Integral approximation Importance sampling (no MAP search)
#> 9    impmap Integral approximation           Importance sampling (MAP)
#> 10  laplace Integral approximation               Laplace approximation
#> 11    qrpem          Stochastic EM          Quasi-Random Parametric EM
#> 12     saem          Stochastic EM         Stochastic Approximation EM
#> 13     npag          Nonparametric         NonParametric Adaptive Grid
#> 14      npb          Nonparametric                 Nonparametric Bayes
#> 15     emvi       Machine learning                      Variational EM
#> 16     fbvi       Machine learning    Full Bayes variational inference
#> 17      vae       Machine learning        Variational autoencoder NLME
#> 18   bobyqa Optimizer (NLM family)            BOBYQA (derivative-free)
#> 19 lbfgsb3c Optimizer (NLM family)                            L-BFGS-B
#> 20    n1qn1 Optimizer (NLM family)                        n1qn1 (BFGS)
#> 21   newuoa Optimizer (NLM family)            NEWUOA (derivative-free)
#> 22      nlm Optimizer (NLM family)                    nlm quasi-Newton
#> 23   nlminb Optimizer (NLM family)                         PORT nlminb
#> 24      nls Optimizer (NLM family)             nonlinear least squares
#> 25    optim Optimizer (NLM family)          Nelder-Mead / BFGS (optim)
#> 26   uobyqa Optimizer (NLM family)            UOBYQA (derivative-free)
```
