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
#>         est                   type
#> 1        fo             Linearized
#> 2      foce             Linearized
#> 3     focei             Linearized
#> 4     focep             Linearized
#> 5       foi             Linearized
#> 6      nlme             Linearized
#> 7       agq Integral approximation
#> 8       imp Integral approximation
#> 9    impmap Integral approximation
#> 10  laplace Integral approximation
#> 11    qrpem          Stochastic EM
#> 12     saem          Stochastic EM
#> 13     npag          Nonparametric
#> 14      npb          Nonparametric
#> 15     advi       Machine learning
#> 16      vae       Machine learning
#> 17   bobyqa Optimizer (NLM family)
#> 18 lbfgsb3c Optimizer (NLM family)
#> 19    n1qn1 Optimizer (NLM family)
#> 20   newuoa Optimizer (NLM family)
#> 21      nlm Optimizer (NLM family)
#> 22   nlminb Optimizer (NLM family)
#> 23      nls Optimizer (NLM family)
#> 24    optim Optimizer (NLM family)
#> 25   uobyqa Optimizer (NLM family)
#>                                        description
#> 1                                      First-Order
#> 2               First-Order Conditional Estimation
#> 3                            FOCE with Interaction
#> 4              FOCE+ (residual at conditional eta)
#> 5                     First-Order with Interaction
#> 6               Lindstrom-Bates alternating (nlme)
#> 7                     Adaptive Gaussian Quadrature
#> 8              Importance sampling (no MAP search)
#> 9                        Importance sampling (MAP)
#> 10                           Laplace approximation
#> 11                      Quasi-Random Parametric EM
#> 12                     Stochastic Approximation EM
#> 13                     NonParametric Adaptive Grid
#> 14                             Nonparametric Bayes
#> 15 Automatic Differentiation Variational Inference
#> 16                    Variational autoencoder NLME
#> 17                        BOBYQA (derivative-free)
#> 18                                        L-BFGS-B
#> 19                                    n1qn1 (BFGS)
#> 20                        NEWUOA (derivative-free)
#> 21                                nlm quasi-Newton
#> 22                                     PORT nlminb
#> 23                         nonlinear least squares
#> 24                      Nelder-Mead / BFGS (optim)
#> 25                        UOBYQA (derivative-free)
```
