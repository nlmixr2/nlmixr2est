
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nlmixr: an R package for population PKPD modeling

<!-- badges: start -->

[![R build
status](https://github.com/nlmixr2/nlmixr2/workflows/R-CMD-check/badge.svg)](https://github.com/nlmixr2/nlmixr2/actions)
[![Codecov test
coverage](https://codecov.io/gh/nlmixr2/nlmixr2/branch/master/graph/badge.svg)](https://codecov.io/gh/nlmixr2/nlmixr2?branch=master)
[![CRAN
checks](https://cranchecks.info/badges/summary/nlmixr)](https://cran.r-project.org/web/checks/check_results_nlmixr.html)
[![CRAN total
downloads](https://cranlogs.r-pkg.org/badges/grand-total/nlmixr)](https://cran.r-project.org/package=nlmixr)
[![CRAN monthly
downloads](https://cranlogs.r-pkg.org/badges/nlmixr)](https://cran.r-project.org/package=nlmixr)
<!-- badges: end -->

------------------------------------------------------------------------

![nlmixr](logo.png)

`nlmixr` is an R package for fitting general dynamic models,
pharmacokinetic (PK) models and pharmacokinetic-pharmacodynamic (PKPD)
models in particular, with either individual data or population data.
The nlme and SAEM estimation routines can be accessed using a universal
user interface (UUI), that provides universal model and parameter
definition syntax and results in a fit object that can be used as input
into the `Xpose` package. Running nlmixr using the UUI is described in
[this
vignette](https://nlmixr2.github.io/nlmixr2/articles/running_nlmixr.html).

Under the hood `nlmixr` has five main modules:

1.  `dynmodel()` and its mcmc cousin `dynmodel.mcmc()` for nonlinear
    dynamic models of individual data;
2.  `nlme_lin_cmpt()`for one to three linear compartment models of
    population data with first order absorption, or i.v. bolus, or i.v.
    infusion using the nlme algorithm;
3.  `nlme_ode()` for general dynamic models defined by ordinary
    differential equations (ODEs) of population data using the nlme
    algorithm;
4.  `saem_fit` for general dynamic models defined by ordinary
    differential equations (ODEs) of population data by the Stochastic
    Approximation Expectation-Maximization (SAEM) algorithm;  
5.  `gnlmm` for generalized non-linear mixed-models (possibly defined by
    ordinary differential equations) of population data by the adaptive
    Gaussian quadrature algorithm.

A few utilities to facilitate population model building are also
included in `nlmixr`.

Documentation can be found at <https://nlmixr2.github.io/nlmixr2/>, and
we maintain a comprehensive and ever-growing guide to using `nlmixr` at
our [bookdown
site](https://nlmixr2.github.io/nlmixr_bookdown/index.html).

More examples and the associated data files are available at
<https://github.com/nlmixr2/nlmixr2/tree/master/vignettes>.

We recommend you have a look at
[`rxode2`](https://nlmixr2.github.io/rxode2/articles/RxODE-intro.html),
the engine upon which `nlmixr` depends, as well as
[`xpose.nlmixr`](https://github.com/nlmixr2/xpose.nlmixr), which
provides a link to the seminal nonlinear mixed-effects model diagnostics
package [`xpose`](https://uupharmacometrics.github.io/xpose/), and
[`shinyMixR`](https://github.com/RichardHooijmaijers/shinyMixR), which
provides a means to build a project-centric workflow around nlmixr from
the R command line and from a streamlined
[`shiny`](https://shiny.rstudio.com/) front-end application. Members of
the nlmixr team also contribute to the
[`ggPMX`](https://github.com/ggPMXdevelopment/ggPMX),
[`xgxr`](https://github.com/Novartis/xgxr) and
[`pmxTools`](https://github.com/kestrel99/pmxTools) packages. For PKPD
modeling (with ODE and dosing history) with [Stan](http://mc-stan.org/),
check out Yuan Xiong’s package
[`PMXStan`](https://github.com/yxiong1/pmxstan).

## Installation

When on CRAN, you can install the released version of nlmixr from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("nlmixr2")
```

### Installing the Development Version

First, install [rxode2](https://github.com/nlmixr2/rxode2#installation).
If not on Windows, you will also need to follow the instructions to
install
[symengine](https://github.com/symengine/symengine.R#installation).

And then install the development version of nlmixr2 from
[GitHub](https://github.com/) with:

``` r
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("nlmixr2/nlmixr2")
```
