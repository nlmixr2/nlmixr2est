# nlmixr2 fits population PK and PKPD non-linear mixed effects models.

nlmixr2 is an R package for fitting population pharmacokinetic (PK) and
pharmacokinetic-pharmacodynamic (PKPD) models.

## Usage

``` r
nlmixr2(
  object,
  data,
  est = NULL,
  control = list(),
  table = tableControl(),
  ...,
  save = NULL,
  envir = parent.frame()
)

nlmixr(
  object,
  data,
  est = NULL,
  control = list(),
  table = tableControl(),
  ...,
  save = NULL,
  envir = parent.frame()
)

# S3 method for class '`function`'
nlmixr2(
  object,
  data = NULL,
  est = NULL,
  control = NULL,
  table = tableControl(),
  ...,
  save = NULL,
  envir = parent.frame()
)

# S3 method for class 'rxUi'
nlmixr2(
  object,
  data = NULL,
  est = NULL,
  control = NULL,
  table = tableControl(),
  ...,
  save = NULL,
  envir = parent.frame()
)

# S3 method for class 'nlmixr2FitCore'
nlmixr2(
  object,
  data = NULL,
  est = NULL,
  control = NULL,
  table = tableControl(),
  ...,
  save = NULL,
  envir = parent.frame()
)

# S3 method for class 'nlmixr2FitData'
nlmixr2(
  object,
  data = NULL,
  est = NULL,
  control = NULL,
  table = tableControl(),
  ...,
  save = NULL,
  envir = parent.frame()
)
```

## Arguments

- object:

  Fitted object or function specifying the model.

- data:

  nlmixr data

- est:

  estimation method (all methods are shown by \`nlmixr2AllEst()\`).
  Methods can be added for other tools

- control:

  The estimation control object. These are expected to be different for
  each type of estimation method

- table:

  The output table control object (like \`tableControl()\`)

- ...:

  Other parameters

- save:

  Boolean to save a nlmixr2 object in a rds file in the working
  directory. If `NULL`, uses option "nlmixr2.save"

- envir:

  Environment where the nlmixr object/function is evaluated before
  running the estimation routine.

## Value

Either a nlmixr2 model or a nlmixr2 fit object

## Details

The nlmixr2 generalized function allows common access to the nlmixr2
estimation routines.

The nlmixr object has the following fields:

|                 |      |                                                                                                                                                                                                                                                           |
|-----------------|------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| Field           | Note | Description                                                                                                                                                                                                                                               |
| censInfo        |      | Gives the censorng information abot the fit (the type of censoring that was seend and handled in the dataset)                                                                                                                                             |
| conditionNumber |      | Condition number, that is the highest divided by the lowest eigenvalue in the population covariance matrix                                                                                                                                                |
| cor             |      | Correlation matrix                                                                                                                                                                                                                                        |
| cov             |      | Variance-covariance matrix                                                                                                                                                                                                                                |
| covMethod       |      | Method used to calculate covariance of the fixed effects                                                                                                                                                                                                  |
| dataLloq        |      | Gives the lloq from the dataset (average) when cesoring has occured; Requires the fit to have a table step                                                                                                                                                |
| dataMergeFull   |      | Full data merge with the fit output and the original dataset; Also includes nlmixrLlikObs which includes the individual observation contribution to the likelihood                                                                                        |
| dataMergeInner  |      | Inner data merge with the fit output and the original dataset; Also includes nlmixrLlikObs which includes the individual observation contribution to the likelihood                                                                                       |
| dataMergeLeft   |      | Left data merge with the fit output and the original dataset; Also includes nlmixrLlikObs which includes the individual observation contribution to the likelihood                                                                                        |
| dataMergeRight  |      | Right data merge with the fit output and the original dataset; Also includes nlmixrLlikObs which includes the individual observation contribution to the likelihood                                                                                       |
| dataUloq        |      | Gives the uloq from the dataset (average) when censoring has occured; requires the fit to have a table step                                                                                                                                               |
| env             |      | This is the environment where all the information for the fit is stored outside of the data-frame. It is an R environment hence \$env                                                                                                                     |
| runInfo         |      | This returns a list of all the warnings or fit information                                                                                                                                                                                                |
| rxControl       |      | Integration options used to control rxode2                                                                                                                                                                                                                |
| scaleInfo       |      | The scaling factors used for nlmixr2 estimation in focei; The can be changed by foceiControl(scaleC=…) if you think these are unreasonable. It also tells the Gill83 outcome of trying to find the best step size (High gradient error, bad gradient etc) |
| seed            |      | This is the initial seed used for saem                                                                                                                                                                                                                    |
| shrink          |      | This is a table of shrinkages for all the individual ETAs as well as the variance shrinkage as well as summary statistics for the ETAs and Residual Error                                                                                                 |
| simInfo         |      | This returns a list of all the fit information used for a traditional rxode2 simulation, which you can tweak yourself if you wish                                                                                                                         |
| table           |      | These are the table options that were used when generating the table output (were CWRES included, etc                                                                                                                                                     |
| theta           |      | Estimates for eta for each individual                                                                                                                                                                                                                     |
| time            |      | Duration of different parts of the analysis (e.g. setup, optimization, calculation of covariance, etc.)                                                                                                                                                   |
| ui              |      | Final estimates for the model                                                                                                                                                                                                                             |
| atol            | n2r  | Absolute tolerance that NONMEM specified; will be used when solving                                                                                                                                                                                       |
| dfObs           | n2r  | Degrees of freedom by observation                                                                                                                                                                                                                         |
| dfSub           | n2r  | Degrees of freedom by subject                                                                                                                                                                                                                             |
| etaData         | n2r  | Subject level IIV values                                                                                                                                                                                                                                  |
| ipredAtol       | n2r  | Absolute tolerance difference between NONMEM and rxode2 individual predictions                                                                                                                                                                            |
| ipredCompare    | n2r  | Data frame with ipred values                                                                                                                                                                                                                              |
| ipredRtol       | n2r  | Relative tolerance difference between NONMEM and rxode2 individual predictions                                                                                                                                                                            |
| nonmemData      | n2r  | Original dataset used for NONMEM analysis                                                                                                                                                                                                                 |
| predAtol        | n2r  | Absolute tolerance difference between NONMEM and rxode2 population predictions                                                                                                                                                                            |
| predCompare     | n2r  | Data frame with pred values                                                                                                                                                                                                                               |
| predRtol        | n2r  | Relative tolerance difference between NONMEM and rxode2 population predictions                                                                                                                                                                            |
| rtol            | n2r  | Relative tolerance that NONMEM specified; will be used when solving                                                                                                                                                                                       |
| sigma           | n2r  | Error model matrix                                                                                                                                                                                                                                        |
| ssAtol          | n2r  | Steady state absolute tolerance that NONMEM specified; will be used for solving.                                                                                                                                                                          |
| ssRtol          | n2r  | Steady state relative tolerance that NONMEM specified will be used for solving                                                                                                                                                                            |
| thetaMat        | n2r  | Covariance Matrix (matches rxSolve(thetaMat=)                                                                                                                                                                                                             |

n2r - These fields are added when a NONMEM model is imported using
`nonmem2rx()`

## nlmixr modeling mini-language

**Rationale**

nlmixr estimation routines each have their own way of specifying models.
Often the models are specified in ways that are most intuitive for one
estimation routine, but do not make sense for another estimation
routine. Sometimes, legacy estimation routines like
[`nlme`](https://rdrr.io/pkg/nlme/man/nlme.html) have their own syntax
that is outside of the control of the nlmixr package.

The unique syntax of each routine makes the routines themselves easier
to maintain and expand, and allows interfacing with existing packages
that are outside of nlmixr (like
[`nlme`](https://rdrr.io/pkg/nlme/man/nlme.html)). However, a model
definition language that is common between estimation methods, and an
output object that is uniform, will make it easier to switch between
estimation routines and will facilitate interfacing output with external
packages like Xpose.

The nlmixr mini-modeling language, attempts to address this issue by
incorporating a common language. This language is inspired by both R and
NONMEM, since these languages are familiar to many pharmacometricians.

**Initial Estimates and boundaries for population parameters**

nlmixr models are contained in a R function with two blocks: `ini` and
`model`. This R function can be named anything, but is not meant to be
called directly from R. In fact if you try you will likely get an error
such as `Error: could not find function "ini"`.

The `ini` model block is meant to hold the initial estimates for the
model, and the boundaries of the parameters for estimation routines that
support boundaries (note nlmixr's `saem` and `nlme` do not currently
support parameter boundaries).

To explain how these initial estimates are specified we will start with
an annotated example:

    f <- function(){ ## Note the arguments to the function are currently
                     ## ignored by nlmixr
        ini({
            ## Initial conditions for population parameters (sometimes
            ## called theta parameters) are defined by either `<-` or '='
            lCl <- 1.6      #log Cl (L/hr)
            ## Note that simple expressions that evaluate to a number are
            ## OK for defining initial conditions (like in R)
            lVc = log(90)  #log V (L)
            ## Also a comment on a parameter is captured as a parameter label
            lKa <- 1 #log Ka (1/hr)
            ## Bounds may be specified by c(lower, est, upper), like NONMEM:
            ## Residuals errors are assumed to be population parameters
            prop.err <- c(0, 0.2, 1)
        })
        ## The model block will be discussed later
        model({})
    }

As shown in the above examples:

- Simple parameter values are specified as a R-compatible assignment

- Boundaries my be specified by `c(lower, est, upper)`.

- Like NONMEM, `c(lower,est)` is equivalent to `c(lower,est,Inf)`

- Also like NONMEM, `c(est)` does not specify a lower bound, and is
  equivalent to specifying the parameter without R's \`c\` function.

- The initial estimates are specified on the variance scale, and in
  analogy with NONMEM, the square roots of the diagonal elements
  correspond to coefficients of variation when used in the exponential
  IIV implementation

These parameters can be named almost any R compatible name. Please note
that:

- Residual error estimates should be coded as population estimates (i.e.
  using an '=' or '\<-' statement, not a '~').

- Naming variables that start with "`_`" are not supported. Note that R
  does not allow variable starting with "`_`" to be assigned without
  quoting them.

- Naming variables that start with "`rx_`" or "`nlmixr_`" is not
  supported since
  [rxode2](https://nlmixr2.github.io/rxode2/reference/rxode2.html) and
  nlmixr2 use these prefixes internally for certain estimation routines
  and calculating residuals.

- Variable names are case sensitive, just like they are in R. "`CL`" is
  not the same as "`Cl`".

**Initial Estimates for between subject error distribution (NONMEM's
\$OMEGA)**

In mixture models, multivariate normal individual deviations from the
population parameters are estimated (in NONMEM these are called `eta`
parameters). Additionally the variance/covariance matrix of these
deviations is also estimated (in NONMEM this is the OMEGA matrix). These
also have initial estimates. In nlmixr these are specified by the \`~\`
operator that is typically used in R for "modeled by", and was chosen to
distinguish these estimates from the population and residual error
parameters.

Continuing the prior example, we can annotate the estimates for the
between subject error distribution

    f <- function(){
        ini({
            lCl <- 1.6      #log Cl (L/hr)
            lVc = log(90)  #log V (L)
            lKa <- 1 #log Ka (1/hr)
            prop.err <- c(0, 0.2, 1)
            ## Initial estimate for ka IIV variance
            ## Labels work for single parameters
            eta.ka ~ 0.1 # BSV Ka

            ## For correlated parameters, you specify the names of each
            ## correlated parameter separated by a addition operator `+`
            ## and the left handed side specifies the lower triangular
            ## matrix initial of the covariance matrix.
            eta.cl + eta.vc ~ c(0.1,
                                0.005, 0.1)
            ## Note that labels do not currently work for correlated
            ## parameters.  Also do not put comments inside the lower
            ## triangular matrix as this will currently break the model.
        })
        ## The model block will be discussed later
        model({})
    }

As shown in the above examples:

- Simple variances are specified by the variable name and the estimate
  separated by \`~\`.

- Correlated parameters are specified by the sum of the variable labels
  and then the lower triangular matrix of the covariance is specified on
  the left handed side of the equation. This is also separated by \`~\`.

Currently the model syntax does not allow comments inside the lower
triangular matrix.

**Model Syntax for ODE based models (NONMEM's \$PK, \$PRED, \$DES and
\$ERROR)**

Once the initialization block has been defined, you can define a model
in terms of the defined variables in the `ini` block. You can also mix
in RxODE blocks into the model.

The current method of defining a nlmixr model is to specify the
parameters, and then possibly the RxODE lines:

Continuing describing the syntax with an annotated example:

    f <- function(){
        ini({
            lCl <- 1.6      #log Cl (L/hr)
            lVc <- log(90)   #log Vc (L)
            lKA <- 0.1      #log Ka (1/hr)
            prop.err <- c(0, 0.2, 1)
            eta.Cl ~ 0.1 ## BSV Cl
            eta.Vc ~ 0.1 ## BSV Vc
            eta.KA ~ 0.1 ## BSV Ka
        })
        model({
            ## First parameters are defined in terms of the initial estimates
            ## parameter names.
            Cl <- exp(lCl + eta.Cl)
            Vc = exp(lVc + eta.Vc)
            KA <- exp(lKA + eta.KA)
            ## After the differential equations are defined
            kel <- Cl / Vc;
            d/dt(depot)    = -KA*depot;
            d/dt(centr)  =  KA*depot-kel*centr;
            ## And the concentration is then calculated
            cp = centr / Vc;
            ## Last, nlmixr is told that the plasma concentration follows
            ## a proportional error (estimated by the parameter prop.err)
            cp ~ prop(prop.err)
        })
    }

A few points to note:

- Parameters are often defined before the differential equations.

- The differential equations, parameters and error terms are in a single
  block, instead of multiple sections.

- State names, calculated variables cannot start with either "`rx_`" or
  "`nlmixr_`" since these are used internally in some estimation
  routines.

- Errors are specified using the \`~\`. Currently you can use either
  `add(parameter)` for additive error, prop(parameter) for proportional
  error or `add(parameter1) + prop(parameter2)` for additive plus
  proportional error. You can also specify `norm(parameter)` for the
  additive error, since it follows a normal distribution.

- Some routines, like `saem` require parameters in terms of
  `Pop.Parameter + Individual.Deviation.Parameter + Covariate*Covariate.Parameter`.
  The order of these parameters do not matter. This is similar to
  NONMEM's mu-referencing, though not quite so restrictive.

- The type of parameter in the model is determined by the initial block;
  Covariates used in the model are missing in the `ini` block. These
  variables need to be present in the modeling dataset for the model to
  run.

**Model Syntax for solved PK systems**

Solved PK systems are also currently supported by nlmixr with the
\`linCmt()\` pseudo-function. An annotated example of a solved system is
below:

\##'

    f <- function(){
        ini({
            lCl <- 1.6      #log Cl (L/hr)
            lVc <- log(90)   #log Vc (L)
            lKA <- 0.1      #log Ka (1/hr)
            prop.err <- c(0, 0.2, 1)
            eta.Cl ~ 0.1 ## BSV Cl
            eta.Vc ~ 0.1 ## BSV Vc
            eta.KA ~ 0.1 ## BSV Ka
        })
        model({
            Cl <- exp(lCl + eta.Cl)
            Vc = exp(lVc + eta.Vc)
            KA <- exp(lKA + eta.KA)
            ## Instead of specifying the ODEs, you can use
            ## the linCmt() function to use the solved system.
            ##
            ## This function determines the type of PK solved system
            ## to use by the parameters that are defined.  In this case
            ## it knows that this is a one-compartment model with first-order
            ## absorption.
            linCmt() ~ prop(prop.err)
        })
    }

A few things to keep in mind:

- While RxODE allows mixing of solved systems and ODEs, this has not
  been implemented in nlmixr yet.

- The solved systems implemented are the one, two and three compartment
  models with or without first-order absorption. Each of the models
  support a lag time with a tlag parameter.

- In general the linear compartment model figures out the model by the
  parameter names. nlmixr currently knows about numbered volumes, Vc/Vp,
  Clearances in terms of both Cl and Q/CLD. Additionally nlmixr knows
  about elimination micro-constants (ie K12). Mixing of these parameters
  for these models is currently not supported.

**Checking model syntax**

After specifying the model syntax you can check that nlmixr is
interpreting it correctly by using the `nlmixr` function on it.

Using the above function we can get:

    > nlmixr(f)
    ## 1-compartment model with first-order absorption in terms of Cl
    ## Initialization:
    ################################################################################
    Fixed Effects ($theta):
        lCl     lVc     lKA
    1.60000 4.49981 0.10000

    Omega ($omega):
         [,1] [,2] [,3]
    [1,]  0.1  0.0  0.0
    [2,]  0.0  0.1  0.0
    [3,]  0.0  0.0  0.1

    ## Model:
    ################################################################################
    Cl <- exp(lCl + eta.Cl)
    Vc = exp(lVc + eta.Vc)
    KA <- exp(lKA + eta.KA)
    ## Instead of specifying the ODEs, you can use
    ## the linCmt() function to use the solved system.
    ##
    ## This function determines the type of PK solved system
    ## to use by the parameters that are defined.  In this case
    ## it knows that this is a one-compartment model with first-order
    ## absorption.
    linCmt() ~ prop(prop.err)

In general this gives you information about the model (what type of
solved system/RxODE), initial estimates as well as the code for the
model block.

**Using the model syntax for estimating a model**

Once the model function has been created, you can use it and a dataset
to estimate the parameters for a model given a dataset.

This dataset has to have RxODE compatible events IDs. Both Monolix and
NONMEM use a a very similar standard to what nlmixr can support.

Once the data has been converted to the appropriate format, you can use
the `nlmixr` function to run the appropriate code.

The method to estimate the model is:

    fit <- nlmixr(model.function, dataset, est="est", control=estControl(options))

Currently `nlme` and `saem` are implemented. For example, to run the
above model with `saem`, we could have the following:

    > f <- function(){
        ini({
            lCl <- 1.6      #log Cl (L/hr)
            lVc <- log(90)   #log Vc (L)
            lKA <- 0.1      #log Ka (1/hr)
            prop.err <- c(0, 0.2, 1)
            eta.Cl ~ 0.1 ## BSV Cl
            eta.Vc ~ 0.1 ## BSV Vc
            eta.KA ~ 0.1 ## BSV Ka
        })
        model({
            ## First parameters are defined in terms of the initial estimates
            ## parameter names.
            Cl <- exp(lCl + eta.Cl)
            Vc = exp(lVc + eta.Vc)
            KA <- exp(lKA + eta.KA)
            ## After the differential equations are defined
            kel <- Cl / Vc;
            d/dt(depot)    = -KA*depot;
            d/dt(centr)  =  KA*depot-kel*centr;
            ## And the concentration is then calculated
            cp = centr / Vc;
            ## Last, nlmixr is told that the plasma concentration follows
            ## a proportional error (estimated by the parameter prop.err)
            cp ~ prop(prop.err)
        })
    }
    > fit.s <- nlmixr(f,d,est="saem",control=saemControl(n.burn=50,n.em=100,print=50));
    Compiling RxODE differential equations...done.
    c:/Rtools/mingw_64/bin/g++  -I"c:/R/R-34~1.1/include" -DNDEBUG     -I"d:/Compiler/gcc-4.9.3/local330/include"  -Ic:/nlmixr/inst/include -Ic:/R/R-34~1.1/library/STANHE~1/include -Ic:/R/R-34~1.1/library/Rcpp/include -Ic:/R/R-34~1.1/library/RCPPAR~1/include -Ic:/R/R-34~1.1/library/RCPPEI~1/include -Ic:/R/R-34~1.1/library/BH/include   -O2 -Wall  -mtune=core2 -c saem3090757b4bd1x64.cpp -o saem3090757b4bd1x64.o
    In file included from c:/R/R-34~1.1/library/RCPPAR~1/include/armadillo:52:0,
                     from c:/R/R-34~1.1/library/RCPPAR~1/include/RcppArmadilloForward.h:46,
                     from c:/R/R-34~1.1/library/RCPPAR~1/include/RcppArmadillo.h:31,
                     from saem3090757b4bd1x64.cpp:1:
    c:/R/R-34~1.1/library/RCPPAR~1/include/armadillo_bits/compiler_setup.hpp:474:96: note: #pragma message: WARNING: use of OpenMP disabled; this compiler doesn't support OpenMP 3.0+
       #pragma message ("WARNING: use of OpenMP disabled; this compiler doesn't support OpenMP 3.0+")
                                                                                                    ^
    c:/Rtools/mingw_64/bin/g++ -shared -s -static-libgcc -o saem3090757b4bd1x64.dll tmp.def saem3090757b4bd1x64.o c:/nlmixr/R/rx_855815def56a50f0e7a80e48811d947c_x64.dll -Lc:/R/R-34~1.1/bin/x64 -lRblas -Lc:/R/R-34~1.1/bin/x64 -lRlapack -lgfortran -lm -lquadmath -Ld:/Compiler/gcc-4.9.3/local330/lib/x64 -Ld:/Compiler/gcc-4.9.3/local330/lib -Lc:/R/R-34~1.1/bin/x64 -lR
    done.
    1:    1.8174   4.6328   0.0553   0.0950   0.0950   0.0950   0.6357
    50:    1.3900   4.2039   0.0001   0.0679   0.0784   0.1082   0.1992
    100:    1.3894   4.2054   0.0107   0.0686   0.0777   0.1111   0.1981
    150:    1.3885   4.2041   0.0089   0.0683   0.0778   0.1117   0.1980
    Using sympy via SnakeCharmR
    ## Calculate ETA-based prediction and error derivatives:
    Calculate Jacobian...................done.
    Calculate sensitivities.......
    done.
    ## Calculate d(f)/d(eta)
    ## ...
    ## done
    ## ...
    ## done
    The model-based sensitivities have been calculated
    Calculating Table Variables...
    done

The options for `saem` are controlled by
[`saemControl`](https://nlmixr2.github.io/nlmixr2est/reference/saemControl.md).
You may wish to make sure the minimization is complete in the case of
`saem`. You can do that with `traceplot` which shows the iteration
history with the divided by burn-in and EM phases. In this case, the
burn in seems reasonable; you may wish to increase the number of
iterations in the EM phase of the estimation. Overall it is probably a
semi-reasonable solution.

**nlmixr output objects**

In addition to unifying the modeling language sent to each of the
estimation routines, the outputs currently have a unified structure.

You can see the fit object by typing the object name:

    > fit.s
     -- nlmixr SAEM fit (ODE); OBJF calculated from FOCEi approximation -------------
          OBJF      AIC      BIC Log-likelihood Condition Number
      62337.09 62351.09 62399.01      -31168.55          82.6086

     -- Time (sec; fit.s$time): -----------------------------------------------------
               saem setup Likelihood Calculation covariance table
     elapsed 430.25 31.64                   1.19          0  3.44

     -- Parameters (fit.s$par.fixed): -----------------------------------------------
                  Parameter Estimate     SE  
     lCl      log Cl (L/hr)     1.39 0.0240  1.73       4.01 (3.83, 4.20)    26.6
     lVc         log Vc (L)     4.20 0.0256 0.608       67.0 (63.7, 70.4)    28.5
     lKA      log Ka (1/hr)  0.00924 0.0323  349.      1.01 (0.947, 1.08)    34.3
     prop.err      prop.err    0.198                             19.8
              Shrink(SD)
     lCl          0.248
     lVc           1.09
     lKA           4.19
     prop.err      1.81

       No correlations in between subject variability (BSV) matrix
       Full BSV covariance (fit.s$omega) or correlation (fit.s$omega.R; diagonals=SDs)
       Distribution stats (mean/skewness/kurtosis/p-value) available in fit.s$shrink

     -- Fit Data (object fit.s is a modified data.frame): ---------------------------
     # A tibble: 6,947 x 22
       ID     TIME    DV  PRED    RES    WRES IPRED  IRES  IWRES CPRED   CRES
     * <fct> <dbl> <dbl> <dbl>  <dbl>   <dbl> <dbl> <dbl>  <dbl> <dbl>  <dbl>
     1 1      0.25  205.  198.   6.60  0.0741  189.  16.2  0.434  198.   6.78
     2 1      0.5   311.  349. -38.7  -0.261   330. -19.0 -0.291  349. -38.3
     3 1      0.75  389.  464. -74.5  -0.398   434. -45.2 -0.526  463. -73.9
     # ... with 6,944 more rows, and 11 more variables: CWRES <dbl>, eta.Cl <dbl>,
     #   eta.Vc <dbl>, eta.KA <dbl>, depot <dbl>, centr <dbl>, Cl <dbl>, Vc <dbl>,
     #   KA <dbl>, kel <dbl>, cp <dbl>

This example shows what is typical printout of a nlmixr fit object. The
elements of the fit are:

- The type of fit ([`nlme`](https://rdrr.io/pkg/nlme/man/nlme.html),
  `saem`, etc)

- Metrics of goodness of fit ([`AIC`](https://rdrr.io/r/stats/AIC.html),
  [`BIC`](https://rdrr.io/r/stats/AIC.html), and
  [`logLik`](https://rdrr.io/r/stats/logLik.html)).

  - To align the comparison between methods, the FOCEi likelihood
    objective is calculated regardless of the method used and used for
    goodness of fit metrics.

  - This FOCEi likelihood has been compared to NONMEM's objective
    function and gives the same values (based on the data in Wang 2007)

  - Also note that `saem` does not calculate an objective function, and
    the FOCEi is used as the only objective function for the fit.

  - Even though the objective functions are calculated in the same
    manner, caution should be used when comparing fits from various
    estimation routines.

- The next item is the timing of each of the steps of the fit.

  - These can be also accessed by (`fit.s$time`).

  - As a mnemonic, the access for this item is shown in the printout.
    This is true for almost all of the other items in the printout.

- After the timing of the fit, the parameter estimates are displayed
  (can be accessed by `fit.s$par.fixed`)

  - While the items are rounded for R printing, each estimate without
    rounding is still accessible by the \`\$\` syntax. For example, the
    \`\$Untransformed\` gives the untransformed parameter values.

  - The Untransformed parameter takes log-space parameters and
    back-transforms them to normal parameters. Not the CIs are listed on
    the back-transformed parameter space.

  - Proportional Errors are converted to

- Omega block (accessed by `fit.s$omega`)

- The table of fit data. Please note:

  - A nlmixr fit object is actually a data frame. Saving it as a Rdata
    object and then loading it without nlmixr will just show the data by
    itself. Don't worry; the fit information has not vanished, you can
    bring it back by simply loading nlmixr, and then accessing the data.

  - Special access to fit information (like the `$omega`) needs nlmixr
    to extract the information.

  - If you use the `$` to access information, the order of precedence
    is:

    - Fit data from the overall data.frame

    - Information about the parsed nlmixr model (via `$uif`)

    - Parameter history if available (via `$par.hist` and
      `$par.hist.stacked`)

    - Fixed effects table (via `$par.fixed`)

    - Individual differences from the typical population parameters (via
      `$eta`)

    - Fit information from the list of information generated during the
      post-hoc residual calculation.

    - Fit information from the environment where the post-hoc residual
      were calculated

    - Fit information about how the data and options interacted with the
      specified model (such as estimation options or if the solved
      system is for an infusion or an IV bolus).

  - While the printout may displays the data as a `data.table` object or
    `tbl` object, the data is NOT any of these objects, but rather a
    derived data frame.

  - Since the object *is* a data.frame, you can treat it like one.

In addition to the above properties of the fit object, there are a few
additional that may be helpful for the modeler:

- `$theta` gives the fixed effects parameter estimates (in NONMEM the
  `theta`s). This can also be accessed in
  [`fixed.effects`](https://rdrr.io/pkg/nlme/man/fixed.effects.html)
  function. Note that the residual variability is treated as a fixed
  effect parameter and is included in this list.

- `$eta` gives the random effects parameter estimates, or in NONMEM the
  `eta`s. This can also be accessed in using the
  [`random.effects`](https://rdrr.io/pkg/nlme/man/random.effects.html)
  function.

## Author

Matthew L. Fidler

## Examples

``` r
# \donttest{

one.cmt <- function() {
 ini({
   ## You may label each parameter with a comment
   tka <- 0.45 # Ka
   tcl <- log(c(0, 2.7, 100)) # Log Cl
   ## This works with interactive models
   ## You may also label the preceding line with label("label text")
   tv <- 3.45; label("log V")
   ## the label("Label name") works with all models
   eta.ka ~ 0.6
   eta.cl ~ 0.3
   eta.v ~ 0.1
   add.sd <- 0.7
   prop.sd <- 0.01
 })
 model({
   ka <- exp(tka + eta.ka)
   cl <- exp(tcl + eta.cl)
   v <- exp(tv + eta.v)
   linCmt() ~ add(add.sd) + prop(prop.sd)
 })
}

# fitF <- nlmixr(one.cmt, theo_sd, "focei")

fitS <- nlmixr(one.cmt, theo_sd, "saem")
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#>  
#>  
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of saem model...
#> ✔ done
#> → finding duplicate expressions in saem model...
#> ✔ done
#> ℹ calculate uninformed etas
#> ℹ done
#> params:  tka tcl tv  V(eta.ka)   V(eta.cl)   V(eta.v)    add.sd  prop.sd
#> Calculating covariance matrix
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of saem model...
#> ✔ done
#> → finding duplicate expressions in saem predOnly model 0...
#> → finding duplicate expressions in saem predOnly model 1...
#> → finding duplicate expressions in saem predOnly model 2...
#> → optimizing duplicate expressions in saem predOnly model 2...
#> ✔ done
#>  
#>  
#> → Calculating residuals/tables
#> ✔ done
#> → compress origData in nlmixr2 object, save 6584
#> → compress parHistData in nlmixr2 object, save 8376
#> → compress phiM in nlmixr2 object, save 313376

# }
```
