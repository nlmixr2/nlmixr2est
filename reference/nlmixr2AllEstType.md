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
#> 11    fsaem          Stochastic EM
#> 12    qrpem          Stochastic EM
#> 13     saem          Stochastic EM
#> 14     npag          Nonparametric
#> 15      npb          Nonparametric
#> 16     advi       Machine learning
#> 17      vae       Machine learning
#> 18   bobyqa Optimizer (NLM family)
#> 19 lbfgsb3c Optimizer (NLM family)
#> 20    n1qn1 Optimizer (NLM family)
#> 21   newuoa Optimizer (NLM family)
#> 22      nlm Optimizer (NLM family)
#> 23   nlminb Optimizer (NLM family)
#> 24      nls Optimizer (NLM family)
#> 25    optim Optimizer (NLM family)
#> 26   uobyqa Optimizer (NLM family)
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
#> 11                     Fast SAEM (Karimi-Lavielle)
#> 12                      Quasi-Random Parametric EM
#> 13                     Stochastic Approximation EM
#> 14                     NonParametric Adaptive Grid
#> 15                             Nonparametric Bayes
#> 16 Automatic Differentiation Variational Inference
#> 17                    Variational autoencoder NLME
#> 18                        BOBYQA (derivative-free)
#> 19                                        L-BFGS-B
#> 20                                    n1qn1 (BFGS)
#> 21                        NEWUOA (derivative-free)
#> 22                                nlm quasi-Newton
#> 23                                     PORT nlminb
#> 24                         nonlinear least squares
#> 25                      Nelder-Mead / BFGS (optim)
#> 26                        UOBYQA (derivative-free)
```
