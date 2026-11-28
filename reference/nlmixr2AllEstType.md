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
#>          est                   type
#> 1         fo             Linearized
#> 2       foce             Linearized
#> 3      focei             Linearized
#> 4      focep             Linearized
#> 5        foi             Linearized
#> 6        its             Linearized
#> 7       nlme             Linearized
#> 8        agq Integral approximation
#> 9        imp Integral approximation
#> 10    impmap Integral approximation
#> 11   laplace Integral approximation
#> 12 stieltjes Integral approximation
#> 13     qrpem          Stochastic EM
#> 14      saem          Stochastic EM
#> 15      npag          Nonparametric
#> 16       npb          Nonparametric
#> 17      emvi       Machine learning
#> 18      fbvi       Machine learning
#> 19       vae       Machine learning
#> 20    bobyqa Optimizer (NLM family)
#> 21  lbfgsb3c Optimizer (NLM family)
#> 22     n1qn1 Optimizer (NLM family)
#> 23    newuoa Optimizer (NLM family)
#> 24       nlm Optimizer (NLM family)
#> 25    nlminb Optimizer (NLM family)
#> 26       nls Optimizer (NLM family)
#> 27     optim Optimizer (NLM family)
#> 28     trust Optimizer (NLM family)
#> 29    uobyqa Optimizer (NLM family)
#>                                      description
#> 1                                    First-Order
#> 2             First-Order Conditional Estimation
#> 3                          FOCE with Interaction
#> 4            FOCE+ (residual at conditional eta)
#> 5                   First-Order with Interaction
#> 6                            Iterative Two-Stage
#> 7             Lindstrom-Bates alternating (nlme)
#> 8                   Adaptive Gaussian Quadrature
#> 9            Importance sampling (no MAP search)
#> 10                     Importance sampling (MAP)
#> 11                         Laplace approximation
#> 12                Two-stage Stieltjes quadrature
#> 13                    Quasi-Random Parametric EM
#> 14                   Stochastic Approximation EM
#> 15                   NonParametric Adaptive Grid
#> 16                           Nonparametric Bayes
#> 17                                Variational EM
#> 18              Full Bayes variational inference
#> 19                  Variational autoencoder NLME
#> 20                      BOBYQA (derivative-free)
#> 21                                      L-BFGS-B
#> 22                                  n1qn1 (BFGS)
#> 23                      NEWUOA (derivative-free)
#> 24                              nlm quasi-Newton
#> 25                                   PORT nlminb
#> 26                       nonlinear least squares
#> 27                    Nelder-Mead / BFGS (optim)
#> 28 Trust-region Newton (C++-resident, RcppTrust)
#> 29                      UOBYQA (derivative-free)
```
