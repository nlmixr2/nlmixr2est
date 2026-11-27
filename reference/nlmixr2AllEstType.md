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
#> 6       nlme             Linearized
#> 7        agq Integral approximation
#> 8        imp Integral approximation
#> 9     impmap Integral approximation
#> 10   laplace Integral approximation
#> 11 stieltjes Integral approximation
#> 12     qrpem          Stochastic EM
#> 13      saem          Stochastic EM
#> 14      npag          Nonparametric
#> 15       npb          Nonparametric
#> 16      emvi       Machine learning
#> 17      fbvi       Machine learning
#> 18       vae       Machine learning
#> 19    bobyqa Optimizer (NLM family)
#> 20  lbfgsb3c Optimizer (NLM family)
#> 21     n1qn1 Optimizer (NLM family)
#> 22    newuoa Optimizer (NLM family)
#> 23       nlm Optimizer (NLM family)
#> 24    nlminb Optimizer (NLM family)
#> 25       nls Optimizer (NLM family)
#> 26     optim Optimizer (NLM family)
#> 27     trust Optimizer (NLM family)
#> 28    uobyqa Optimizer (NLM family)
#>                                      description
#> 1                                    First-Order
#> 2             First-Order Conditional Estimation
#> 3                          FOCE with Interaction
#> 4            FOCE+ (residual at conditional eta)
#> 5                   First-Order with Interaction
#> 6             Lindstrom-Bates alternating (nlme)
#> 7                   Adaptive Gaussian Quadrature
#> 8            Importance sampling (no MAP search)
#> 9                      Importance sampling (MAP)
#> 10                         Laplace approximation
#> 11                Two-stage Stieltjes quadrature
#> 12                    Quasi-Random Parametric EM
#> 13                   Stochastic Approximation EM
#> 14                   NonParametric Adaptive Grid
#> 15                           Nonparametric Bayes
#> 16                                Variational EM
#> 17              Full Bayes variational inference
#> 18                  Variational autoencoder NLME
#> 19                      BOBYQA (derivative-free)
#> 20                                      L-BFGS-B
#> 21                                  n1qn1 (BFGS)
#> 22                      NEWUOA (derivative-free)
#> 23                              nlm quasi-Newton
#> 24                                   PORT nlminb
#> 25                       nonlinear least squares
#> 26                    Nelder-Mead / BFGS (optim)
#> 27 Trust-region Newton (C++-resident, RcppTrust)
#> 28                      UOBYQA (derivative-free)
```
