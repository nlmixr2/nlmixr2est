.regFloat1 <- rex::rex(
  or(
    group(some_of("0":"9"), ".", any_of("0":"9")),
    group(any_of("0":"9"), ".", some_of("0":"9"))
  ),
  maybe(group(one_of("E", "e"), maybe(one_of("+", "-")), some_of("0":"9")))
)
.regFloat2 <- rex::rex(some_of("0":"9"), one_of("E", "e"), maybe(one_of("-", "+")), some_of("0":"9"))
.regDecimalint <- rex::rex(or("0", group("1":"9", any_of("0":"9"))))
.regNum <- rex::rex(maybe("-"), or(.regDecimalint, .regFloat1, .regFloat2))

use.utf <- function() {
  opt <- getOption("cli.unicode", NULL)
  if (!is.null(opt)) {
    isTRUE(opt)
  } else {
    l10n_info()$`UTF-8` && !is.latex()
  }
}

is.latex <- function() {
  if (!("knitr" %in% loadedNamespaces())) {
    return(FALSE)
  }
  get("is_latex_output", asNamespace("knitr"))()
}

#' Control Options for FOCEi
#'
#' @param sigdig Optimization significant digits. This controls:
#'
#' \itemize{
#'
#'  \item The tolerance of the inner and outer optimization is \code{10^-sigdig}
#'
#'  \item The tolerance of the ODE solvers is
#'  \code{0.5*10^(-sigdig-2)}; For the sensitivity equations and
#'  steady-state solutions the default is \code{0.5*10^(-sigdig-1.5)}
#'  (sensitivity changes only applicable for liblsoda)
#'
#'  \item The tolerance of the boundary check is \code{5 * 10 ^ (-sigdig + 1)}
#'
#'  \item The significant figures that some tables are rounded to.
#' }
#'
#' @param atolSens Sensitivity atol, can be different than atol with
#'     liblsoda.  This allows a less accurate solve for gradients (if desired)
#'
#' @param rtolSens Sensitivity rtol, can be different than rtol with
#'     liblsoda.  This allows a less accurate solve for gradients (if desired)
#'
#' @param ssAtol Steady state absolute tolerance (atol) for calculating if steady-state
#'     has been archived.
#'
#' @param ssRtol Steady state relative tolerance (rtol) for
#'     calculating if steady-state has been achieved.
#'
#' @param ssAtolSens Sensitivity absolute tolerance (atol) for
#'     calculating if steady state has been achieved for sensitivity compartments.
#'
#' @param ssRtolSens Sensitivity relative tolerance (rtol) for
#'     calculating if steady state has been achieved for sensitivity compartments.
#'
#' @param epsilon Precision of estimate for n1qn1 optimization.
#'
#' @param maxstepsOde Maximum number of steps for ODE solver.
#'
#' @param print Integer representing when the outer step is
#'     printed. When this is 0 or do not print the iterations.  1 is
#'     print every function evaluation (default), 5 is print every 5
#'     evaluations.
#'
#' @param scaleTo Scale the initial parameter estimate to this value.
#'     By default this is 1.  When zero or below, no scaling is performed.
#'
#' @param scaleObjective Scale the initial objective function to this
#'     value.  By default this is 1.
#'
#' @param derivEps Forward difference tolerances, which is a
#'     vector of relative difference and absolute difference.  The
#'     central/forward difference step size h is calculated as:
#'
#'         \code{h = abs(x)*derivEps[1] + derivEps[2]}
#'
#' @param derivMethod indicates the method for calculating
#'     derivatives of the outer problem.  Currently supports
#'     "switch", "central" and "forward" difference methods.  Switch
#'     starts with forward differences.  This will switch to central
#'     differences when abs(delta(OFV)) <= derivSwitchTol and switch
#'     back to forward differences when abs(delta(OFV)) >
#'     derivSwitchTol.
#'
#' @param derivSwitchTol The tolerance to switch forward to central
#'     differences.
#'
#' @param covDerivMethod indicates the method for calculating the
#'     derivatives while calculating the covariance components
#'     (Hessian and S).
#'
#' @param covMethod Method for calculating covariance.  In this
#'     discussion, R is the Hessian matrix of the objective
#'     function. The S matrix is the sum of individual
#'     gradient cross-product (evaluated at the individual empirical
#'     Bayes estimates).
#'
#' \itemize{
#'
#'  \item "\code{r,s}" Uses the sandwich matrix to calculate the
#'  covariance, that is: \code{solve(R) \%*\% S \%*\% solve(R)}
#'
#'  \item "\code{r}" Uses the Hessian matrix to calculate the
#'  covariance as \code{2 \%*\% solve(R)}
#'
#'  \item "\code{s}" Uses the cross-product matrix to calculate the
#'  covariance as \code{4 \%*\% solve(S)}
#'
#'  \item "" Does not calculate the covariance step.
#' }
#'
#' @param covTryHarder If the R matrix is non-positive definite and
#'     cannot be corrected to be non-positive definite try estimating
#'     the Hessian on the unscaled parameter space.
#'
#' @param hessEps is a double value representing the epsilon for the Hessian calculation.
#'
#' @param centralDerivEps Central difference tolerances.  This is a
#'     numeric vector of relative difference and absolute difference.
#'     The central/forward difference step size h is calculated as:
#'
#'         \code{h = abs(x)*derivEps[1] + derivEps[2]}
#'
#' @param lbfgsLmm An integer giving the number of BFGS updates
#'     retained in the "L-BFGS-B" method, It defaults to 7.
#'
#' @param lbfgsPgtol is a double precision variable.
#'
#'     On entry pgtol >= 0 is specified by the user.  The iteration
#'     will stop when:
#'
#'        \code{max(\| proj g_i \| i = 1, ..., n) <= lbfgsPgtol}
#'
#'     where pg_i is the ith component of the projected gradient.
#'
#'     On exit pgtol is unchanged.  This defaults to zero, when the
#'     check is suppressed.
#'
#' @param lbfgsFactr Controls the convergence of the "L-BFGS-B"
#'     method.  Convergence occurs when the reduction in the
#'     objective is within this factor of the machine
#'     tolerance. Default is 1e10, which gives a tolerance of about
#'     \code{2e-6}, approximately 4 sigdigs.  You can check your
#'     exact tolerance by multiplying this value by
#'     \code{.Machine$double.eps}
#'
#' @param diagXform This is the transformation used on the diagonal
#'     of the \code{chol(solve(omega))}. This matrix and values are the
#'     parameters estimated in FOCEi. The possibilities are:
#'
#' \itemize{
#'  \item \code{sqrt} Estimates the sqrt of the diagonal elements of \code{chol(solve(omega))}.  This is the default method.
#'
#'  \item \code{log} Estimates the log of the diagonal elements of \code{chol(solve(omega))}
#'
#'  \item \code{identity} Estimates the diagonal elements without any transformations
#' }
#' @param sumProd Is a boolean indicating if the model should change
#'     multiplication to high precision multiplication and sums to
#'     high precision sums using the PreciseSums package.  By default
#'     this is \code{FALSE}.
#'
#'
#' @param optExpression Optimize the rxode2 expression to speed up
#'     calculation. By default this is turned on.
#'
#' @param ci Confidence level for some tables.  By default this is
#'     0.95 or 95\% confidence.
#'
#' @param useColor Boolean indicating if focei can use ASCII color codes
#'
#' @param boundTol Tolerance for boundary issues.
#'
#' @param calcTables This boolean is to determine if the foceiFit
#'     will calculate tables. By default this is \code{TRUE}
#'
#' @param ... Ignored parameters
#'
#' @param maxInnerIterations Number of iterations for n1qn1
#'     optimization.
#'
#' @param maxOuterIterations Maximum number of L-BFGS-B optimization
#'     for outer problem.
#'
#' @param n1qn1nsim Number of function evaluations for n1qn1
#'     optimization.
#'
#' @param eigen A boolean indicating if eigenvectors are calculated
#'     to include a condition number calculation.
#'
#' @param addPosthoc Boolean indicating if posthoc parameters are
#'     added to the table output.
#'
#' @param printNcol Number of columns to printout before wrapping
#'     parameter estimates/gradient
#'
#' @param noAbort Boolean to indicate if you should abort the FOCEi
#'     evaluation if it runs into troubles.  (default TRUE)
#'
#' @param interaction Boolean indicate FOCEi should be used (TRUE)
#'     instead of FOCE (FALSE)
#'
#' @param cholSEOpt Boolean indicating if the generalized Cholesky
#'     should be used while optimizing.
#'
#' @param cholSECov Boolean indicating if the generalized Cholesky
#'     should be used while calculating the Covariance Matrix.
#'
#' @param fo is a boolean indicating if this is a FO approximation routine.
#'
#' @param cholSEtol tolerance for Generalized Cholesky
#'     Decomposition.  Defaults to suggested (.Machine$double.eps)^(1/3)
#'
#' @param cholAccept Tolerance to accept a Generalized Cholesky
#'     Decomposition for a R or S matrix.
#'
#' @param outerOpt optimization method for the outer problem
#'
#' @param innerOpt optimization method for the inner problem (not
#'     implemented yet.)
#'
#' @param stateTrim Trim state amounts/concentrations to this value.
#'
#' @param resetEtaP represents the p-value for reseting the
#'     individual ETA to 0 during optimization (instead of the saved
#'     value).  The two test statistics used in the z-test are either
#'     chol(omega^-1) \%*\% eta or eta/sd(allEtas).  A p-value of 0
#'     indicates the ETAs never reset.  A p-value of 1 indicates the
#'     ETAs always reset.
#'
#' @param resetThetaP represents the p-value for reseting the
#'     population mu-referenced THETA parameters based on ETA drift
#'     during optimization, and resetting the optimization.  A
#'     p-value of 0 indicates the THETAs never reset.  A p-value of 1
#'     indicates the THETAs always reset and is not allowed.  The
#'     theta reset is checked at the beginning and when nearing a
#'     local minima.  The percent change in objective function where
#'     a theta reset check is initiated is controlled in
#'     \code{resetThetaCheckPer}.
#'
#' @param resetThetaCheckPer represents objective function
#'     \% percentage below which resetThetaP is checked.
#'
#' @param resetThetaFinalP represents the p-value for reseting the
#'     population mu-referenced THETA parameters based on ETA drift
#'     during optimization, and resetting the optimization one final time.
#'
#' @param resetHessianAndEta is a boolean representing if the
#'     individual Hessian is reset when ETAs are reset using the
#'     option \code{resetEtaP}.
#'
#' @param diagOmegaBoundUpper This represents the upper bound of the
#'     diagonal omega matrix.  The upper bound is given by
#'     diag(omega)*diagOmegaBoundUpper.  If
#'     \code{diagOmegaBoundUpper} is 1, there is no upper bound on
#'     Omega.
#'
#' @param diagOmegaBoundLower This represents the lower bound of the
#'     diagonal omega matrix.  The lower bound is given by
#'     diag(omega)/diagOmegaBoundUpper.  If
#'     \code{diagOmegaBoundLower} is 1, there is no lower bound on
#'     Omega.
#'
#' @param rhobeg Beginning change in parameters for bobyqa algorithm
#'     (trust region).  By default this is 0.2 or 20% of the initial
#'     parameters when the parameters are scaled to 1. rhobeg and
#'     rhoend must be set to the initial and final values of a trust
#'     region radius, so both must be positive with 0 < rhoend <
#'     rhobeg. Typically rhobeg should be about one tenth of the
#'     greatest expected change to a variable.  Note also that
#'     smallest difference abs(upper-lower) should be greater than or
#'     equal to rhobeg*2. If this is not the case then rhobeg will be
#'     adjusted.
#'
#' @param rhoend The smallest value of the trust region radius that
#'     is allowed. If not defined, then 10^(-sigdig-1) will be used.
#'
#' @param npt The number of points used to approximate the objective
#'     function via a quadratic approximation for bobyqa. The value
#'     of npt must be in the interval [n+2,(n+1)(n+2)/2] where n is
#'     the number of parameters in par. Choices that exceed 2*n+1 are
#'     not recommended. If not defined, it will be set to 2*n + 1
#' @param eval.max Number of maximum evaluations of the objective function
#'
#' @param iter.max Maximum number of iterations allowed.
#'
#' @param rel.tol Relative tolerance before nlminb stops.
#'
#' @param x.tol X tolerance for nlmixr2 optimizers
#'
#' @param abstol Absolute tolerance for nlmixr2 optimizer
#'
#' @param reltol  tolerance for nlmixr2
#'
#' @param gillK The total number of possible steps to determine the
#'     optimal forward/central difference step size per parameter (by
#'     the Gill 1983 method).  If 0, no optimal step size is
#'     determined.  Otherwise this is the optimal step size
#'     determined.
#'
#' @param gillRtol The relative tolerance used for Gill 1983
#'     determination of optimal step size.
#'
#' @param scaleType The scaling scheme for nlmixr2.  The supported types are:
#'
#' \itemize{
#' \item \code{nlmixr2}  In this approach the scaling is performed by the following equation:
#'
#'    v_{scaled} = (v_{current} - v_{init})/scaleC[i] + scaleTo
#'
#' The \code{scaleTo} parameter is specified by the \code{normType},
#' and the scales are specified by \code{scaleC}.
#'
#' \item \code{norm} This approach uses the simple scaling provided
#'     by the \code{normType} argument.
#'
#' \item \code{mult} This approach does not use the data
#' normalization provided by \code{normType}, but rather uses
#' multiplicative scaling to a constant provided by the \code{scaleTo}
#' argument.
#'
#'   In this case:
#'
#'   v_{scaled} = v_{current}/v_{init}*scaleTo
#'
#' \item \code{multAdd} This approach changes the scaling based on
#' the parameter being specified.  If a parameter is defined in an
#' exponential block (ie exp(theta)), then it is scaled on a
#' linearly, that is:
#'
#'   v_{scaled} = (v_{current}-v_{init}) + scaleTo
#'
#' Otherwise the parameter is scaled multiplicatively.
#'
#'    v_{scaled} = v_{current}/v_{init}*scaleTo
#'
#' }
#'
#' @param scaleC The scaling constant used with
#'     \code{scaleType=nlmixr2}.  When not specified, it is based on
#'     the type of parameter that is estimated.  The idea is to keep
#'     the derivatives similar on a log scale to have similar
#'     gradient sizes.  Hence parameters like log(exp(theta)) would
#'     have a scaling factor of 1 and log(theta) would have a scaling
#'     factor of ini_value (to scale by 1/value; ie
#'     d/dt(log(ini_value)) = 1/ini_value or scaleC=ini_value)
#'
#'    \itemize{
#'
#'    \item For parameters in an exponential (ie exp(theta)) or
#'    parameters specifying powers, boxCox or yeoJohnson
#'    transformations , this is 1.
#'
#'    \item For additive, proportional, lognormal error structures,
#'    these are given by 0.5*abs(initial_estimate)
#'
#'    \item Factorials are scaled by abs(1/digamma(inital_estimate+1))
#'
#'    \item parameters in a log scale (ie log(theta)) are transformed
#'    by log(abs(initial_estimate))*abs(initial_estimate)
#'
#'    }
#'
#'    These parameter scaling coefficients are chose to try to keep
#'    similar slopes among parameters.  That is they all follow the
#'    slopes approximately on a log-scale.
#'
#'    While these are chosen in a logical manner, they may not always
#'    apply.  You can specify each parameters scaling factor by this
#'    parameter if you wish.
#'
#' @param scaleC0 Number to adjust the scaling factor by if the initial
#'     gradient is zero.
#'
#' @param scaleCmax Maximum value of the scaleC to prevent overflow.
#'
#' @param scaleCmin Minimum value of the scaleC to prevent underflow.
#'
#' @param normType This is the type of parameter
#'     normalization/scaling used to get the scaled initial values
#'     for nlmixr2.  These are used with \code{scaleType} of.
#'
#'     With the exception of \code{rescale2}, these come
#'     from
#'     \href{https://en.wikipedia.org/wiki/Feature_scaling}{Feature
#'     Scaling}. The \code{rescale2} The rescaling is the same type
#'     described in the
#'     \href{http://apmonitor.com/me575/uploads/Main/optimization_book.pdf}{OptdesX}
#'     software manual.
#'
#'     In general, all all scaling formula can be described by:
#'
#'     v_{scaled} = (v_{unscaled}-C_{1})/C_{2}
#'
#'     Where
#'
#'
#'     The other data normalization approaches follow the following formula
#'
#'     v_{scaled} = (v_{unscaled}-C_{1})/C_{2};
#'
#' \itemize{
#'
#' \item \code{rescale2} This scales all parameters from (-1 to 1).
#'     The relative differences between the parameters are preserved
#'     with this approach and the constants are:
#'
#'     C_{1} = (max(all unscaled values)+min(all unscaled values))/2
#'
#'     C_{2} = (max(all unscaled values) - min(all unscaled values))/2
#'
#'
#' \item \code{rescale} or min-max normalization. This rescales all
#'     parameters from (0 to 1).  As in the \code{rescale2} the
#'     relative differences are preserved.  In this approach:
#'
#'     C_{1} = min(all unscaled values)
#'
#'     C_{2} = max(all unscaled values) - min(all unscaled values)
#'
#'
#' \item \code{mean} or mean normalization.  This rescales to center
#'     the parameters around the mean but the parameters are from 0
#'     to 1.  In this approach:
#'
#'     C_{1} = mean(all unscaled values)
#'
#'     C_{2} = max(all unscaled values) - min(all unscaled values)
#'
#' \item \code{std} or standardization.  This standardizes by the mean
#'      and standard deviation.  In this approach:
#'
#'     C_{1} = mean(all unscaled values)
#'
#'     C_{2} = sd(all unscaled values)
#'
#' \item \code{len} or unit length scaling.  This scales the
#'    parameters to the unit length.  For this approach we use the Euclidean length, that
#'    is:
#'
#'     C_{1} = 0
#'
#'     C_{2} = sqrt(v_1^2 + v_2^2 + ... + v_n^2)
#'
#'
#' \item \code{constant} which does not perform data normalization. That is
#'
#'     C_{1} = 0
#'
#'     C_{2} = 1
#'
#' }
#'
#' @param gillStep When looking for the optimal forward difference
#'     step size, this is This is the step size to increase the
#'     initial estimate by.  So each iteration the new step size =
#'     (prior step size)*gillStep
#'
#' @param gillFtol The gillFtol is the gradient error tolerance that
#'     is acceptable before issuing a warning/error about the gradient estimates.
#'
#' @param gillKcov The total number of possible steps to determine
#'     the optimal forward/central difference step size per parameter
#'     (by the Gill 1983 method) during the covariance step.  If 0,
#'     no optimal step size is determined.  Otherwise this is the
#'     optimal step size determined.
#'
#' @param gillStepCov When looking for the optimal forward difference
#'     step size, this is This is the step size to increase the
#'     initial estimate by.  So each iteration during the covariance
#'     step is equal to the new step size = (prior step size)*gillStepCov
#'
#' @param gillFtolCov The gillFtol is the gradient error tolerance
#'     that is acceptable before issuing a warning/error about the
#'     gradient estimates during the covariance step.
#'
#' @param rmatNorm A parameter to normalize gradient step size by the
#'     parameter value during the calculation of the R matrix
#'
#' @param smatNorm A parameter to normalize gradient step size by the
#'     parameter value during the calculation of the S matrix
#'
#' @param covGillF Use the Gill calculated optimal Forward difference
#'     step size for the instead of the central difference step size
#'     during the central difference gradient calculation.
#'
#' @param optGillF Use the Gill calculated optimal Forward difference
#'     step size for the instead of the central difference step size
#'     during the central differences for optimization.
#'
#' @param covSmall The covSmall is the small number to compare
#'     covariance numbers before rejecting an estimate of the
#'     covariance as the final estimate (when comparing sandwich vs
#'     R/S matrix estimates of the covariance).  This number controls
#'     how small the variance is before the covariance matrix is
#'     rejected.
#'
#' @param adjLik In nlmixr2, the objective function matches NONMEM's
#'     objective function, which removes a 2*pi constant from the
#'     likelihood calculation. If this is TRUE, the likelihood
#'     function is adjusted by this 2*pi factor.  When adjusted this
#'     number more closely matches the likelihood approximations of
#'     nlme, and SAS approximations.  Regardless of if this is turned
#'     on or off the objective function matches NONMEM's objective
#'     function.
#'
#' @param gradTrim The parameter to adjust the gradient to if the
#'     |gradient| is very large.
#'
#' @param gradCalcCentralSmall A small number that represents the value
#'     where |grad| < gradCalcCentralSmall where forward differences
#'     switch to central differences.
#'
#' @param gradCalcCentralLarge A large number that represents the value
#'     where |grad| > gradCalcCentralLarge where forward differences
#'     switch to central differences.
#'
#' @param etaNudge By default initial ETA estimates start at zero;
#'   Sometimes this doesn't optimize appropriately.  If this value is
#'   non-zero, when the n1qn1 optimization didn't perform
#'   appropriately, reset the Hessian, and nudge the ETA up by this
#'   value; If the ETA still doesn't move, nudge the ETA down by this
#'   value. By default this value is qnorm(1-0.05/2)*1/sqrt(3), the
#'   first of the Gauss Quadrature numbers times by the 0.95\% normal
#'   region. If this is not successful try the second eta nudge
#'   number (below).  If +-etaNudge2 is not successful, then assign
#'   to zero and do not optimize any longer
#'
#' @param etaNudge2 This is the second eta nudge.  By default it is
#'   qnorm(1-0.05/2)*sqrt(3/5), which is the n=3 quadrature point
#'   (excluding zero) times by the 0.95\% normal region
#'
#' @param maxOdeRecalc Maximum number of times to reduce the ODE
#'     tolerances and try to resolve the system if there was a bad
#'     ODE solve.
#'
#' @param odeRecalcFactor The factor to increase the rtol/atol with
#'     bad ODE solving.
#'
#' @param repeatGillMax If the tolerances were reduced when
#'     calculating the initial Gill differences, the Gill difference
#'     is repeated up to a maximum number of times defined by this
#'     parameter.
#'
#' @param stickyRecalcN The number of bad ODE solves before reducing
#'     the atol/rtol for the rest of the problem.
#'
#' @param nRetries If FOCEi doesn't fit with the current parameter
#'     estimates, randomly sample new parameter estimates and restart
#'     the problem.  This is similar to 'PsN' resampling.
#'
#' @param eventFD Finite difference step for forward or central
#'     difference estimation of event-based gradients
#'
#' @param eventType Event gradient type for dosing events; Can be
#'   "gill", "central" or "forward"
#'
#' @param gradProgressOfvTime This is the time for a single objective
#'     function evaluation (in seconds) to start progress bars on gradient evaluations
#'
#' @param singleOde This option allows a single ode model to include
#'   the PK parameter information instead of splitting it into a
#'   function and a rxode2 model
#'
#' @param badSolveObjfAdj The objective function adjustment when the
#'   ODE system cannot be solved.  It is based on each individual bad
#'   solve.
#'
#' @inheritParams configsaem
#' @inheritParams rxode2::rxSolve
#' @inheritParams minqa::bobyqa
#'
#' @details
#'
#' Note this uses the R's L-BFGS-B in \code{\link{optim}} for the
#' outer problem and the BFGS \code{\link[n1qn1]{n1qn1}} with that
#' allows restoring the prior individual Hessian (for faster
#' optimization speed).
#'
#' However the inner problem is not scaled.  Since most eta estimates
#' start near zero, scaling for these parameters do not make sense.
#'
#' This process of scaling can fix some ill conditioning for the
#' unscaled problem.  The covariance step is performed on the
#' unscaled problem, so the condition number of that matrix may not
#' be reflective of the scaled problem's condition-number.
#'
#' @author Matthew L. Fidler
#'
#' @return The control object that changes the options for the FOCEi
#'   family of estimation methods
#'
#' @seealso \code{\link{optim}}
#' @seealso \code{\link[n1qn1]{n1qn1}}
#' @seealso \code{\link[rxode2]{rxSolve}}
#' @export
##' @param odeRecalcFactor
##' @param maxOdeRecalc
foceiControl <- function(sigdig = 3, ...,
                         epsilon = NULL, # 1e-4,
                         maxInnerIterations = 1000,
                         maxOuterIterations = 5000,
                         n1qn1nsim = NULL,
                         method = c("liblsoda", "lsoda", "dop853"),
                         transitAbs = NULL, atol = NULL, rtol = NULL,
                         atolSens = NULL, rtolSens = NULL,
                         ssAtol = NULL, ssRtol = NULL, ssAtolSens = NULL, ssRtolSens = NULL,
                         minSS = 10L, maxSS = 1000L,
                         maxstepsOde = 500000L, hmin = 0L, hmax = NA_real_, hini = 0, maxordn = 12L, maxords = 5L, cores,
                         covsInterpolation = c("locf", "linear", "nocb", "midpoint"),
                         print = 1L,
                         printNcol = floor((getOption("width") - 23) / 12),
                         scaleTo = 1.0,
                         scaleObjective = 0,
                         normType = c("rescale2", "mean", "rescale", "std", "len", "constant"),
                         scaleType = c("nlmixr2", "norm", "mult", "multAdd"),
                         scaleCmax = 1e5,
                         scaleCmin = 1e-5,
                         scaleC = NULL,
                         scaleC0 = 1e5,
                         derivEps = rep(20 * sqrt(.Machine$double.eps), 2),
                         derivMethod = c("switch", "forward", "central"),
                         derivSwitchTol = NULL,
                         covDerivMethod = c("central", "forward"),
                         covMethod = c("r,s", "r", "s", ""),
                         hessEps = (.Machine$double.eps)^(1 / 3),
                         eventFD = sqrt(.Machine$double.eps),
                         eventType = c("gill", "central", "forward"),
                         centralDerivEps = rep(20 * sqrt(.Machine$double.eps), 2),
                         lbfgsLmm = 7L,
                         lbfgsPgtol = 0,
                         lbfgsFactr = NULL,
                         eigen = TRUE,
                         addPosthoc = TRUE,
                         diagXform = c("sqrt", "log", "identity"),
                         sumProd = FALSE,
                         optExpression = TRUE,
                         ci = 0.95,
                         useColor = crayon::has_color(),
                         boundTol = NULL,
                         calcTables = TRUE,
                         noAbort = TRUE,
                         interaction = TRUE,
                         cholSEtol = (.Machine$double.eps)^(1 / 3),
                         cholAccept = 1e-3,
                         resetEtaP = 0.15,
                         resetThetaP = 0.05,
                         resetThetaFinalP = 0.15,
                         diagOmegaBoundUpper = 5, # diag(omega) = diag(omega)*diagOmegaBoundUpper; =1 no upper
                         diagOmegaBoundLower = 100, # diag(omega) = diag(omega)/diagOmegaBoundLower; = 1 no lower
                         cholSEOpt = FALSE,
                         cholSECov = FALSE,
                         fo = FALSE,
                         covTryHarder = FALSE,
                         ## Ranking based on run 025
                         ## L-BFGS-B: 20970.53 (2094.004    429.535)
                         ## bobyqa: 21082.34 (338.677    420.754)
                         ## lbfgsb3* (modified for tolerances):
                         ## nlminb: 20973.468 (755.821    458.343)
                         ## mma: 20974.20 (Time: Opt: 3000.501 Cov: 467.287)
                         ## slsqp: 21023.89 (Time: Opt: 460.099; Cov: 488.921)
                         ## lbfgsbLG: 20974.74 (Time: Opt: 946.463; Cov:397.537)
                         outerOpt = c("nlminb", "bobyqa", "lbfgsb3c", "L-BFGS-B", "mma", "lbfgsbLG", "slsqp", "Rvmmin"),
                         innerOpt = c("n1qn1", "BFGS"),
                         ##
                         rhobeg = .2,
                         rhoend = NULL,
                         npt = NULL,
                         ## nlminb
                         rel.tol = NULL,
                         x.tol = NULL,
                         eval.max = 4000,
                         iter.max = 2000,
                         abstol = NULL,
                         reltol = NULL,
                         resetHessianAndEta = FALSE,
                         stateTrim = Inf,
                         gillK = 10L,
                         gillStep = 4,
                         gillFtol = 0,
                         gillRtol = sqrt(.Machine$double.eps),
                         gillKcov = 10L,
                         gillStepCov = 2,
                         gillFtolCov = 0,
                         rmatNorm = TRUE,
                         smatNorm = TRUE,
                         covGillF = TRUE,
                         optGillF = TRUE,
                         covSmall = 1e-5,
                         adjLik = TRUE, ## Adjust likelihood by 2pi for FOCEi methods
                         gradTrim = Inf,
                         maxOdeRecalc = 5,
                         odeRecalcFactor = 10^(0.5),
                         gradCalcCentralSmall = 1e-4,
                         gradCalcCentralLarge = 1e4,
                         etaNudge = qnorm(1-0.05/2)/sqrt(3),
                         etaNudge2=qnorm(1-0.05/2) * sqrt(3/5),
                         stiff,
                         nRetries = 3,
                         seed = 42,
                         resetThetaCheckPer = 0.1,
                         etaMat = NULL,
                         repeatGillMax = 3,
                         stickyRecalcN = 5,
                         gradProgressOfvTime = 10,
                         addProp = c("combined2", "combined1"),
                         singleOde = TRUE,
                         badSolveObjfAdj=100) {
  if (is.null(boundTol)) {
    boundTol <- 5 * 10^(-sigdig + 1)
  }
  if (is.null(epsilon)) {
    epsilon <- 10^(-sigdig - 1)
  }
  if (is.null(abstol)) {
    abstol <- 10^(-sigdig - 1)
  }
  if (is.null(reltol)) {
    reltol <- 10^(-sigdig - 1)
  }
  if (is.null(rhoend)) {
    rhoend <- 10^(-sigdig - 1)
  }
  if (is.null(lbfgsFactr)) {
    lbfgsFactr <- 10^(-sigdig - 1) / .Machine$double.eps
  }
  if (is.null(atol)) {
    atol <- 0.5 * 10^(-sigdig - 2)
  }
  if (is.null(rtol)) {
    rtol <- 0.5 * 10^(-sigdig - 2)
  }
  if (is.null(atolSens)) {
    atolSens <- 0.5 * 10^(-sigdig - 1.5)
  }
  if (is.null(rtolSens)) {
    rtolSens <- 0.5 * 10^(-sigdig - 1.5)
  }
  if (is.null(ssAtol)) {
    ssAtol <- 0.5 * 10^(-sigdig - 2)
  }
  if (is.null(ssRtol)) {
    ssRtol <- 0.5 * 10^(-sigdig - 2)
  }
  if (is.null(ssAtolSens)) {
    ssAtolSens <- 0.5 * 10^(-sigdig - 1.5)
  }
  if (is.null(ssRtolSens)) {
    ssRtolSens <- 0.5 * 10^(-sigdig - 1.5)
  }
  if (is.null(rel.tol)) {
    rel.tol <- 10^(-sigdig - 1)
  }
  if (is.null(x.tol)) {
    x.tol <- 10^(-sigdig - 1)
  }
  if (is.null(derivSwitchTol)) {
    derivSwitchTol <- 2 * 10^(-sigdig - 1)
  }
  ## if (is.null(gillRtol)){
  ##     ## FIXME: there is a way to calculate this according to the
  ##     ## Gill paper but it is buried in their optimization book.
  ##     gillRtol <- 10 ^ (-sigdig - 1);
  ## }
  .xtra <- list(...)
  if (is.null(transitAbs) && !is.null(.xtra$transit_abs)) { # nolint
    transitAbs <- .xtra$transit_abs # nolint
  }
  if (missing(covsInterpolation) && !is.null(.xtra$covs_interpolation)) { # nolint
    covsInterpolation <- .xtra$covs_interpolation # nolint
  }
  if (missing(maxInnerIterations) && !is.null(.xtra$max_iterations)) { # nolint
    maxInnerIterations <- .xtra$max_iterations # nolint
  }
  if (!missing(stiff) && missing(method)) {
    if (rxode2::rxIs(stiff, "logical")) {
      if (stiff) {
        method <- "lsoda"
        warning("stiff=TRUE has been replaced with method = \"lsoda\".")
      } else {
        method <- "dop853"
        warning("stiff=FALSE has been replaced with method = \"dop853\".")
      }
    }
  } else {
    if (inherits(method, "numeric")) {
      method <- as.integer(method)
    }
    if (!rxode2::rxIs(method, "integer")) {
      if (inherits(method, "character")) {
        method <- match.arg(method)
      } else {
        method <- "liblsoda"
        warning("could not figure out method, using 'liblsoda'")
      }
    }
  }
  ## .methodIdx <- c("lsoda"=1L, "dop853"=0L, "liblsoda"=2L);
  ## method <- as.integer(.methodIdx[method]);
  if (rxode2::rxIs(scaleType, "character")) {
    .scaleTypeIdx <- c("norm" = 1L, "nlmixr2" = 2L, "mult" = 3L, "multAdd" = 4L)
    scaleType <- as.integer(.scaleTypeIdx[match.arg(scaleType)])
  }
  if (rxode2::rxIs(eventType, "character")) {
    .eventTypeIdx <- c("gill" = 1L, "central" = 2L, "forward" = 3L)
    eventType <- as.integer(.eventTypeIdx[match.arg(eventType)])
  }
  if (rxode2::rxIs(normType, "character")) {
    .normTypeIdx <- c("rescale2" = 1L, "rescale" = 2L, "mean" = 3L, "std" = 4L, "len" = 5L, "constant" = 6)
    normType <- as.integer(.normTypeIdx[match.arg(normType)])
  }
  derivMethod <- match.arg(derivMethod)
  .methodIdx <- c("forward" = 0L, "central" = 1L, "switch" = 3L)
  derivMethod <- as.integer(.methodIdx[derivMethod])
  covDerivMethod <- .methodIdx[match.arg(covDerivMethod)]
  if (length(covsInterpolation) > 1) covsInterpolation <- covsInterpolation[1]
  if (!rxode2::rxIs(covsInterpolation, "integer")) {
    covsInterpolation <- tolower(match.arg(
      covsInterpolation,
      c("linear", "locf", "LOCF", "constant", "nocb", "NOCB", "midpoint")
    ))
  }

  ## if (covsInterpolation == "constant") covsInterpolation <- "locf";
  ## covsInterpolation  <- as.integer(which(covsInterpolation == c("linear", "locf", "nocb", "midpoint")) - 1);
  if (missing(cores)) {
    cores <- rxode2::rxCores()
  }
  if (missing(n1qn1nsim)) {
    n1qn1nsim <- 10 * maxInnerIterations + 1
  }
  if (length(covMethod) == 1) {
    if (covMethod == "") {
      covMethod <- 0L
    }
  }
  if (rxode2::rxIs(covMethod, "character")) {
    covMethod <- match.arg(covMethod)
    .covMethodIdx <- c("r,s" = 1L, "r" = 2L, "s" = 3L)
    covMethod <- .covMethodIdx[match.arg(covMethod)]
  }
  .outerOptTxt <- "custom"
  if (rxode2::rxIs(outerOpt, "character")) {
    outerOpt <- match.arg(outerOpt)
    .outerOptTxt <- outerOpt
    if (outerOpt == "bobyqa") {
      rxode2::rxReq("minqa")
      outerOptFun <- .bobyqa
      outerOpt <- -1L
    } else if (outerOpt == "nlminb") {
      outerOptFun <- .nlminb
      outerOpt <- -1L
    } else if (outerOpt == "mma") {
      outerOptFun <- .nloptr
      outerOpt <- -1L
    } else if (outerOpt == "slsqp") {
      outerOptFun <- .slsqp
      outerOpt <- -1L
    } else if (outerOpt == "lbfgsbLG") {
      outerOptFun <- .lbfgsbLG
      outerOpt <- -1L
    } else if (outerOpt == "Rvmmin") {
      outerOptFun <- .Rvmmin
      outerOpt <- -1L
    } else {
      .outerOptIdx <- c("L-BFGS-B" = 0L, "lbfgsb3c" = 1L)
      outerOpt <- .outerOptIdx[outerOpt]
      if (outerOpt == 1L) {
        rxode2::rxReq("lbfgsb3c")
      }
      outerOptFun <- NULL
    }
  } else if (is(outerOpt, "function")) {
    outerOptFun <- outerOpt
    outerOpt <- -1L
  }
  if (rxode2::rxIs(innerOpt, "character")) {
    .innerOptFun <- c("n1qn1" = 1L, "BFGS" = 2L)
    innerOpt <- setNames(.innerOptFun[match.arg(innerOpt)], NULL)
  }

  if (resetEtaP > 0 & resetEtaP < 1) {
    .resetEtaSize <- qnorm(1 - (resetEtaP / 2))
  } else if (resetEtaP <= 0) {
    .resetEtaSize <- Inf
  } else {
    .resetEtaSize <- 0
  }

  if (resetThetaP > 0 & resetThetaP < 1) {
    .resetThetaSize <- qnorm(1 - (resetThetaP / 2))
  } else if (resetThetaP <= 0) {
    .resetThetaSize <- Inf
  } else {
    stop("Cannot always reset THETAs")
  }
  if (resetThetaFinalP > 0 & resetThetaFinalP < 1) {
    .resetThetaFinalSize <- qnorm(1 - (resetThetaFinalP / 2))
  } else if (resetThetaP <= 0) {
    .resetThetaFinalSize <- Inf
  } else {
    stop("Cannot always reset THETAs")
  }
  if (inherits(addProp, "numeric")) {
    if (addProp == 1) {
      addProp <- "combined1"
    } else if (addProp == 2) {
      addProp <- "combined2"
    } else {
      stop("addProp must be 1, 2, \"combined1\" or \"combined2\"", call.=FALSE)
    }
  } else {
    addProp <- match.arg(addProp)
  }
  .ret <- list(
    maxOuterIterations = as.integer(maxOuterIterations),
    maxInnerIterations = as.integer(maxInnerIterations),
    method = method,
    transitAbs = transitAbs,
    atol = atol,
    rtol = rtol,
    atolSens = atolSens,
    rtolSens = rtolSens,
    ssAtol = ssAtol,
    ssRtol = ssRtol,
    ssAtolSens = ssAtolSens,
    ssRtolSens = ssRtolSens,
    minSS = minSS, maxSS = maxSS,
    maxstepsOde = maxstepsOde,
    hmin = hmin,
    hmax = hmax,
    hini = hini,
    maxordn = maxordn,
    maxords = maxords,
    cores = cores,
    covsInterpolation = covsInterpolation,
    n1qn1nsim = as.integer(n1qn1nsim),
    print = as.integer(print),
    lbfgsLmm = as.integer(lbfgsLmm),
    lbfgsPgtol = as.double(lbfgsPgtol),
    lbfgsFactr = as.double(lbfgsFactr),
    scaleTo = scaleTo,
    epsilon = epsilon,
    derivEps = derivEps,
    derivMethod = derivMethod,
    covDerivMethod = covDerivMethod,
    covMethod = covMethod,
    centralDerivEps = centralDerivEps,
    eigen = as.integer(eigen),
    addPosthoc = as.integer(addPosthoc),
    diagXform = match.arg(diagXform),
    sumProd = sumProd,
    optExpression = optExpression,
    outerOpt = as.integer(outerOpt),
    ci = as.double(ci),
    sigdig = as.double(sigdig),
    scaleObjective = as.double(scaleObjective),
    useColor = as.integer(useColor),
    boundTol = as.double(boundTol),
    calcTables = calcTables,
    printNcol = as.integer(printNcol),
    noAbort = as.integer(noAbort),
    interaction = as.integer(interaction),
    cholSEtol = as.double(cholSEtol),
    hessEps = as.double(hessEps),
    cholAccept = as.double(cholAccept),
    resetEtaSize = as.double(.resetEtaSize),
    resetThetaSize = as.double(.resetThetaSize),
    resetThetaFinalSize = as.double(.resetThetaFinalSize),
    diagOmegaBoundUpper = diagOmegaBoundUpper,
    diagOmegaBoundLower = diagOmegaBoundLower,
    cholSEOpt = as.integer(cholSEOpt),
    cholSECov = as.integer(cholSECov),
    fo = as.integer(fo),
    covTryHarder = as.integer(covTryHarder),
    outerOptFun = outerOptFun,
    ## bobyqa
    rhobeg = as.double(rhobeg),
    rhoend = as.double(rhoend),
    npt = npt,
    ## nlminb
    rel.tol = as.double(rel.tol),
    x.tol = as.double(x.tol),
    eval.max = eval.max,
    iter.max = iter.max,
    innerOpt = innerOpt,
    ## BFGS
    abstol = abstol,
    reltol = reltol,
    derivSwitchTol = derivSwitchTol,
    resetHessianAndEta = as.integer(resetHessianAndEta),
    stateTrim = as.double(stateTrim),
    gillK = as.integer(gillK),
    gillKcov = as.integer(gillKcov),
    gillRtol = as.double(gillRtol),
    gillStep = as.double(gillStep),
    gillStepCov = as.double(gillStepCov),
    scaleType = scaleType,
    normType = normType,
    scaleC = scaleC,
    scaleCmin = as.double(scaleCmin),
    scaleCmax = as.double(scaleCmax),
    scaleC0 = as.double(scaleC0),
    outerOptTxt = .outerOptTxt,
    rmatNorm = as.integer(rmatNorm),
    smatNorm = as.integer(smatNorm),
    covGillF = as.integer(covGillF),
    optGillF = as.integer(optGillF),
    gillFtol = as.double(gillFtol),
    gillFtolCov = as.double(gillFtolCov),
    covSmall = as.double(covSmall),
    adjLik = adjLik,
    gradTrim = as.double(gradTrim),
    gradCalcCentralSmall = as.double(gradCalcCentralSmall),
    gradCalcCentralLarge = as.double(gradCalcCentralLarge),
    etaNudge = as.double(etaNudge),
    etaNudge2=as.double(etaNudge2),
    maxOdeRecalc = as.integer(maxOdeRecalc),
    odeRecalcFactor = as.double(odeRecalcFactor),
    nRetries = nRetries,
    seed = seed,
    resetThetaCheckPer = resetThetaCheckPer,
    etaMat = etaMat,
    repeatGillMax = as.integer(repeatGillMax),
    stickyRecalcN = as.integer(max(1, abs(stickyRecalcN))),
    eventFD = eventFD,
    eventType = eventType,
    gradProgressOfvTime = gradProgressOfvTime,
    addProp = addProp,
    singleOde = singleOde,
    badSolveObjfAdj=badSolveObjfAdj,
    ...
  )
  if (!missing(etaMat) && missing(maxInnerIterations)) {
    warning("by supplying etaMat, assume you wish to evaluate at ETAs, so setting maxInnerIterations=0")
    .ret$maxInnerIterations <- 0L
    .ret$etaMat <- etaMat
  } else if (!is.null(etaMat)) {
    .ret$etaMat <- etaMat
  }
  .tmp <- .ret
  .tmp$maxsteps <- maxstepsOde
  .tmp <- do.call(rxode2::rxControl, .tmp)
  .ret$rxControl <- .tmp
  class(.ret) <- "foceiControl"
  return(.ret)
}

.ucminf <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  rxode2::rxReq("ucminf")
  .ctl <- control
  .ctl$stepmax <- control$rhobeg
  .ctl$maxeval <- control$maxOuterIterations
  .ctl <- .ctl[names(.ctl) %in% c("stepmax", "maxeval")]
  .ret <- ucminf::ucminf(par, fn, gr = NULL, ..., control = list(), hessian = 2)
  .ret$x <- .ret$par
  return(.ret)
}

.bobyqa <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- control
  if (is.null(.ctl$npt)) .ctl$npt <- length(par) * 2 + 1
  .ctl$iprint <- 0L
  .ctl <- .ctl[names(.ctl) %in% c("npt", "rhobeg", "rhoend", "iprint", "maxfun")]
  .ret <- minqa::bobyqa(par, fn,
    control = .ctl,
    lower = lower,
    upper = upper
  )
  .ret$x <- .ret$par
  .ret$message <- .ret$msg
  .ret$convergence <- .ret$ierr
  .ret$value <- .ret$fval
  return(.ret)
}

.lbfgsb3c <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .w <- which(names(control) %in% c("trace", "factr", "pgtol", "abstol", "reltol", "lmm", "maxit", "iprint"))
  .control <- control[.w]
  .ret <- lbfgsb3c::lbfgsb3c(par = as.vector(par), fn = fn, gr = gr, lower = lower, upper = upper, control = .control)
  .ret$x <- .ret$par
  return(.ret)
}

.lbfgsbO <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .control <- control[names(control) %in% c("trace", "factr", "pgtol", "abstol", "reltol", "lmm", "maxit", "iprint")]
  .w <- which(sapply(.control, is.null))
  .control <- .control[-.w]
  .ret <- optim(
    par = par, fn = fn, gr = gr, method = "L-BFGS-B",
    lower = lower, upper = upper,
    control = .control, hessian = FALSE
  )
  .ret$x <- .ret$par
  return(.ret)
}

.mymin <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .control <- control[names(control) %in% c(
    "eval.max", "iter.max", "trace", "abs.tol",
    "rel.tol", "x.tol", "xf.tol", "step.min", "step.max", "sing.tol", "scale.init", "diff.g"
  )]

  if (all(lower != -Inf) | all(upper != Inf)) {
    warning("Optimization: Boundaries not used in Nelder-Mead")
  }
  fit <- mymin(par, fn, control = .control)
  fit$message <- c("NON-CONVERGENCE", "NELDER_FTOL_REACHED")[1 + fit$convergence]
  fit$x <- fit$par
  return(fit)
}

.nlminb <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- control
  .ctl <- .ctl[names(.ctl) %in% c(
    "eval.max", "iter.max", "trace", "abs.tol", "rel.tol", "x.tol", "xf.tol", "step.min", "step.max", "sing.tol",
    "scale.inti", "diff.g"
  )]
  .ctl$trace <- 0
  .ret <- stats::nlminb(
    start = par, objective = fn, gradient = gr, hessian = NULL, control = .ctl,
    lower = lower, upper = upper
  )
  .ret$x <- .ret$par
  ## .ret$message   already there.
  ## .ret$convergence already there.
  return(.ret)
}

.Rvmmin <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  ## Also gives unreasonable estimates
  rxode2::rxReq("Rvmmin")
  .masked <- rep_len(1, length(par))
  .ctl <- list(
    maxit = control$maxOuterIterations,
    ## maxfevals
    trace = 0, dowarn = FALSE, checkgrad = FALSE, checkbounds = FALSE,
    keepinputpar = FALSE, eps = control$abstol
  )
  .ret <- Rvmmin::Rvmmin(par = par, fn = fn, gr = gr, lower = lower, upper = upper, bdmsk = .masked, control = list(), ...)
  .ret$x <- .ret$par
  .ret$message <- .ret$message
  return(.ret)
}

.nloptr <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ..., nloptrAlgoritm = "NLOPT_LD_MMA") {
  rxode2::rxReq("nloptr")
  .ctl <- list(
    algorithm = nloptrAlgoritm,
    xtol_rel = control$reltol,
    xtol_abs = rep_len(control$abstol, length(par)),
    ftol_abs = control$abstol,
    ftol_rel = control$reltol,
    print_level = 0,
    check_derivatives = FALSE,
    check_derivatives_print = FALSE,
    maxeval = control$maxOuterIterations
  )
  .ret <- nloptr::nloptr(
    x0 = par, eval_f = fn, eval_grad_f = gr,
    lb = lower, ub = upper,
    opts = .ctl
  )
  .ret$par <- .ret$solution
  .ret$x <- .ret$solution
  .ret$convergence <- .ret$status
  .ret$value <- .ret$objective
  return(.ret)
}

.bobyqaNLopt <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctl <- list(
    algorithm = "NLOPT_LN_BOBYQA",
    xtol_rel = control$reltol,
    xtol_abs = rep_len(control$abstol, length(par)),
    ftol_abs = control$abstol,
    ftol_rel = control$reltol,
    print_level = 0,
    check_derivatives = FALSE,
    check_derivatives_print = FALSE,
    maxeval = control$maxOuterIterations
  )
  .ret <- nloptr::nloptr(
    x0 = par, eval_f = fn,
    lb = lower, ub = upper,
    opts = .ctl
  )
  .ret$par <- .ret$solution
  .ret$x <- .ret$solution
  .ret$convergence <- .ret$status
  .ret$value <- .ret$objective
  return(.ret)
}

.slsqp <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  return(.nloptr(par, fn, gr, lower, upper, control, ..., nloptrAlgoritm = "NLOPT_LD_SLSQP"))
}

.lbfgsbLG <- function(par, fn, gr, lower = -Inf, upper = Inf, control = list(), ...) {
  .ctlLocal <- list(
    algorithm = "NLOPT_LD_LBFGS",
    xtol_rel = control$reltol,
    xtol_abs = rep_len(control$abstol, length(par)),
    ftol_abs = control$abstol,
    ftol_rel = control$reltol,
    print_level = 0,
    check_derivatives = FALSE,
    check_derivatives_print = FALSE,
    maxeval = control$maxOuterIterations
  )
  .ctl <- opts <- list(
    "algorithm" = "NLOPT_LD_AUGLAG",
    xtol_rel = control$reltol,
    xtol_abs = rep_len(control$abstol, length(par)),
    ftol_abs = control$abstol,
    ftol_rel = control$reltol,
    maxeval = control$maxOuterIterations,
    "local_opts" = .ctlLocal,
    "print_level" = 0
  )
  .ret <- nloptr::nloptr(
    x0 = par, eval_f = fn, eval_grad_f = gr,
    lb = lower, ub = upper,
    opts = .ctl
  )
  .ret$par <- .ret$solution
  .ret$x <- .ret$solution
  .ret$convergence <- .ret$status
  .ret$value <- .ret$objective
  return(.ret)
}

#' Get the THETA/ETA lines from rxode2 UI
#'
#' @param rxui This is the rxode2 ui object
#' @return The theta/eta lines
#' @author Matthew L. Fidler
#' @noRd
.uiGetThetaEta <- function(rxui) {
  .iniDf <- rxui$iniDf
  .w <- which(!is.na(.iniDf$ntheta))
  .thetas <- lapply(.w, function(i) {
    eval(parse(text=paste0("quote(", .iniDf$name[i], " <- THETA[", .iniDf$ntheta[i],"])")))
  })
  .etas <- NULL
  .i2 <- .iniDf[-.w, ]
  if (length(.i2$name) > 0) {
    .i2 <- .i2[.i2$neta1 == .i2$neta2, ]
    .etas <- lapply(seq_along(.i2$name), function(i) {
      eval(parse(text=paste0("quote(", .i2$name[i], " <- ETA[", .i2$neta1[i], "])")))
    })
  }
  c(.thetas, .etas)
}

#' Get the THETA/ETA params from the rxode2 UI
#'
#' @param rxui This is the rxode2 ui object
#' @return The params eirxode2 UI
#' @author Matthew L. Fidler
#' @noRd
.uiGetThetaEtaParams <- function(rxui, str=FALSE) {
  .iniDf <- rxui$iniDf
  .w <- which(!is.na(.iniDf$ntheta))
  .thetas <- vapply(.w, function(i) {
    paste0("THETA[", .iniDf$ntheta[i],"]")
  }, character(1), USE.NAMES=FALSE)
  .etas <- NULL
  .i2 <- .iniDf[-.w, ]
  if (length(.i2$name) > 0) {
    .i2 <- .i2[.i2$neta1 == .i2$neta2, ]
    .etas <- vapply(seq_along(.i2$name), function(i) {
      paste0("ETA[", .i2$neta1[i],"]")
    }, character(1), USE.NAMES=FALSE)
  }
  .str <- paste(c(.thetas, .etas, rxui$covariates), collapse=", ")
  if (str) {
    paste0("params(", .str, ")")
  } else {
    eval(parse(text=paste0("quote(params(", .str, "))")))
  }
}

#' @export
rxUiGet.foceiParams <- function(x, ...) {
  .ui <- x[[1]]
  .uiGetThetaEtaParams(.ui, str=TRUE)
}

#' @export
rxUiGet.foceiCmtPreModel <- function(x, ...) {
  .ui <- x[[1]]
  .state <- rxode2::rxState(.ui$mv0)
  if (length(.state) == 0) return("")
  paste(paste0("cmt(", .state, ")"), collapse="\n")
}

# This handles the errors for focei
.createFoceiLineObject <- function(x, line) {
  .predDf <- get("predDf", x)
  if (line > nrow(.predDf)) {
    return(NULL)
  }
  .predLine <- .predDf[line, ]
  .ret <- list(x, .predLine)
  class(.ret) <- c(paste(.predLine$distribution), "rxGetDistributionFoceiLines")
  .ret
}

#' This is a S3 method for getting the distribution lines for a base rxode2 focei problem
#'
#' @param line Parsed rxode2 model environment
#' @return Lines for the focei. This is based
#'   on the idea that the focei parameters are defined
#' @author Matthew Fidler
#' @keywords internal
#' @export
rxGetDistributionFoceiLines <- function(line) {
  UseMethod("rxGetDistributionFoceiLines")
}

#' @export
rxGetDistributionFoceiLines.norm <- function(line) {
  env <- line[[1]]
  pred1 <- line[[2]]
  rxode2::.handleSingleErrTypeNormOrTFoceiBase(env, pred1)
}

#' @export
rxGetDistributionFoceiLines.t <- function(line) {
  stop("t isn't supported yet")
}

#' @export
rxGetDistributionFoceiLines.default  <- function(line) {
  stop("Distribution not supported")
}

#' @export
rxGetDistributionFoceiLines.rxUi <- function(line) {
  .predDf <- get("predDf", line)
  lapply(seq_along(.predDf$cond), function(c){
    .mod <- .createFoceiLineObject(line, c)
    rxGetDistributionFoceiLines(.mod)
  })
}

#' @export
rxUiGet.foceiModel0 <- function(x, ...) {
  .f <- x[[1]]
  rxode2::rxCombineErrorLines(.f, errLines=rxGetDistributionFoceiLines(.f),
                              prefixLines=.uiGetThetaEta(.f),
                              paramsLine=NA, #.uiGetThetaEtaParams(.f),
                              modelVars=TRUE,
                              cmtLines=FALSE,
                              dvidLine=FALSE)
}
#attr(rxUiGet.foceiModel0, "desc") <- "FOCEi model base"

.foceiPrune <- function(x, fullModel=TRUE) {
  .x <- x[[1]]
  .x <- .x$foceiModel0[[-1]]
  .env <- new.env(parent = emptyenv())
  .env$.if <- NULL
  .env$.def1 <- NULL
  if (fullModel) {
    .malert("pruning branches ({.code if}/{.code else}) of full model...")
  } else {
    .malert("pruning branches ({.code if}/{.code else})...")
  }
  .ret <- rxode2::.rxPrune(.x, envir = .env)
  .mv <- rxode2::rxModelVars(.ret)
  ## Need to convert to a function
  if (rxode2::.rxIsLinCmt() == 1L) {
    .vars <- c(.mv$params, .mv$lhs, .mv$slhs)
    .mv <- rxode2::.rxLinCmtGen(length(.mv$state), .vars)
  }
  .msuccess("done")
  rxode2::rxNorm(.mv)
}

.loadSymengine <- function(newmod, promoteLinSens = TRUE, fullModel = FALSE) {
  if (fullModel) {
    .malert("loading full model into {.pkg symengine} environment...")
  } else {
    .malert("loading into {.pkg symengine} environment...")
  }
  rxode2::rxS(newmod, TRUE, promoteLinSens = promoteLinSens)
}

#' @export
rxUiGet.loadPruneSens <- function(x, ...) {
  .loadSymengine(.foceiPrune(x), promoteLinSens = TRUE)
}
#attr(rxUiGet.loadPruneSens, "desc") <- "load sensitivity with linCmt() promoted"

#' @export
rxUiGet.loadPrune <- function(x, ...) {
  .loadSymengine(.foceiPrune(x), promoteLinSens = FALSE)
}
#attr(rxUiGet.loadPrune, "desc") <- "load sensitivity without linCmt() promoted"

.sensEtaOrTheta <- function(s, theta=FALSE) {
  .etaVars <- NULL
  if (theta && exists("..maxTheta", s)) {
    .etaVars <- paste0("THETA_", seq(1, s$..maxTheta), "_")
  } else if (exists("..maxEta", s)) {
    .etaVars <- paste0("ETA_", seq(1, s$..maxEta), "_")
  }
  if (length(.etaVars) == 0L) {
    stop("cannot identify parameters for sensitivity analysis\n   with nlmixr2 an 'eta' initial estimate must use '~'", call. = FALSE)
  }
  .stateVars <- rxode2::rxState(s)
  rxode2::.rxJacobian(s, c(.stateVars, .etaVars))
  rxode2::.rxSens(s, .etaVars)
  s
}

#' @export
rxUiGet.foceiEtaS <- function(x, ..., theta=FALSE) {
  .s <- rxUiGet.loadPruneSens(x, ...)
  .sensEtaOrTheta(.s)
}
#attr(rxUiGet.foceiEtaS, "desc") <- "Get symengine environment with eta sensitivities"


#' @export
rxUiGet.foceiThetaS <- function(x, ..., theta=FALSE) {
  .s <- rxUiGet.loadPruneSens(x, ...)
  .sensEtaOrTheta(.s, theta=TRUE)
}
#attr(rxUiGet.foceiEtaS, "desc") <- "Get symengine environment with eta sensitivities"

#' @export
rxUiGet.foceiHdEta <- function(x, ...) {
  .s <- rxUiGet.foceiEtaS(x)
  .stateVars <- rxode2::rxState(.s)
  # FIXME: take out pred.minus.dv
  .predMinusDv <- rxode2::rxGetControl(x[[1]], "predMinusDv", TRUE)
  .grd <- rxode2::rxExpandFEta_(
    .stateVars, .s$..maxEta,
    ifelse(.predMinusDv, 1L, 2L)
  )
  if (rxode2::.useUtf()) {
    .malert("calculate \u2202(f)/\u2202(\u03B7)")
  } else {
    .malert("calculate d(f)/d(eta)")
  }
  rxode2::rxProgress(dim(.grd)[1])
  on.exit({
    rxode2::rxProgressAbort()
  })
  .any.zero <- FALSE
  .all.zero <- TRUE
  .ret <- apply(.grd, 1, function(x) {
    .l <- x["calc"]
    .l <- eval(parse(text = .l))
    .ret <- paste0(x["dfe"], "=", rxode2::rxFromSE(.l))
    .zErr <- suppressWarnings(try(as.numeric(get(x["dfe"], .s)), silent = TRUE))
    if (identical(.zErr, 0)) {
      .any.zero <<- TRUE
    } else if (.all.zero) {
      .all.zero <<- FALSE
    }
    rxode2::rxTick()
    return(.ret)
  })
  if (.all.zero) {
    stop("none of the predictions depend on 'ETA'", call. = FALSE)
  }
  if (.any.zero) {
    warning("some of the predictions do not depend on 'ETA'", call. = FALSE)
  }
  .s$..HdEta <- .ret
  .s$..pred.minus.dv <- .predMinusDv
  rxode2::rxProgressStop()
  .s
}
attr(rxUiGet.foceiHdEta, "desc") <- "Generate the d(err)/d(eta) values for FO related methods"


#' Finalize inner rxode2 based on symengine saved info
#'
#' @param .s Symengine/rxode2 object
#' @return Nothing
#' @author Matthew L Fidler
#' @noRd
.rxFinalizeInner <- function(.s, sum.prod = FALSE,
                             optExpression = TRUE) {
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  .r <- get("rx_r_", envir = .s)
  .r <- paste0("rx_r_=", rxode2::rxFromSE(.r))
  .yj <- paste(get("rx_yj_", envir = .s))
  .yj <- paste0("rx_yj_~", rxode2::rxFromSE(.yj))
  .lambda <- paste(get("rx_lambda_", envir = .s))
  .lambda <- paste0("rx_lambda_~", rxode2::rxFromSE(.lambda))
  .hi <- paste(get("rx_hi_", envir = .s))
  .hi <- paste0("rx_hi_~", rxode2::rxFromSE(.hi))
  .low <- paste(get("rx_low_", envir = .s))
  .low <- paste0("rx_low_~", rxode2::rxFromSE(.low))
  .ddt <- .s$..ddt
  if (is.null(.ddt)) .ddt <- character(0)
  .sens <- .s$..sens
  if (is.null(.sens)) .sens <- character(0)
  .s$..inner <- paste(c(
    .ddt,
    .sens,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .s$..HdEta,
    .r,
    .s$..REta,
    .s$..stateInfo["statef"],
    .s$..stateInfo["dvid"],
    ""
  ), collapse = "\n")
  if (sum.prod) {
    .malert("stabilizing round off errors in inner problem...")
    .s$..inner <- rxode2::rxSumProdModel(.s$..inner)
    .msuccess("done")
  }
  if (optExpression) {
    .s$..inner <- rxode2::rxOptExpr(.s$..inner, "inner model")
  }
}

#' @export
rxUiGet.foceiEnv <- function(x, ...) {
  .s <- rxUiGet.foceiHdEta(x, ...)
  .stateVars <- rxode2::rxState(.s)
  .grd <- rxode2::rxExpandFEta_(.stateVars, .s$..maxEta, FALSE)
  if (rxode2::.useUtf()) {
    .malert("calculate \u2202(R\u00B2)/\u2202(\u03B7)")
  } else {
    .malert("calculate d(R^2)/d(eta)")
  }
  rxode2::rxProgress(dim(.grd)[1])
  on.exit({
    rxode2::rxProgressAbort()
  })
  .ret <- apply(.grd, 1, function(x) {
    .l <- x["calc"]
    .l <- eval(parse(text = .l))
    .ret <- paste0(x["dfe"], "=", rxode2::rxFromSE(.l))
    rxode2::rxTick()
    return(.ret)
  })

  .s$..REta <- .ret
  rxode2::rxProgressStop()
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  .rxFinalizeInner(.s, .sumProd, .optExpression)
  .rxFinalizePred(.s, .sumProd, .optExpression)
  .s$..outer <- NULL
  .s
}
#attr(rxUiGet.foceiEnv, "desc") <- "Get the focei environment"

#' @export
rxUiGet.foceEnv <- function(x, ...) {
  .s <- rxUiGet.foceiHdEta(x, ...)
  .s$..REta <- NULL
  ## Take etas from rx_r
  eval(parse(text = rxode2::rxRepR0_(.s$..maxEta)))
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  .rxFinalizeInner(.s, .sumProd, .optExpression)
  .rxFinalizePred(.s, .sumProd, .optExpression)
  .s$..outer <- NULL
  .s
}
#attr(rxUiGet.foceEnv, "desc") <- "Get the foce environment"


#' @export
rxUiGet.getEBEEnv <- function(x, ...) {
  .s <- rxUiGet.loadPrune(x, ...)
  .s$..inner <- NULL
  .s$..outer <- NULL
  .sumProd <- rxode2::rxGetControl(x[[1]], "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(x[[1]], "optExpression", TRUE)
  .rxFinalizePred(.s, .sumProd, .optExpression)
  .s
}
#attr(rxUiGet.getEBEEnv, "desc") <- "Get the EBE environment"

.toRxParam <- ""
.toRxDvidCmt <- ""

.toRx <- function(x, msg) {
  if (is.null(x)) {
    return(NULL)
  }
  .malert(msg)
  .ret <- rxode2::rxode2(paste(.toRxParam, x, .toRxDvidCmt))
  .msuccess("done")
  return(.ret)
}

.nullInt <- function(x) {
  if (rxode2::rxIs(x, "integer") || rxode2::rxIs(x, "numeric")) {
    return(as.integer(x))
  } else {
    return(integer(0))
  }
}


.rxFinalizePred <- function(.s, sum.prod = FALSE,
                            optExpression = TRUE) {
  .prd <- get("rx_pred_", envir = .s)
  .prd <- paste0("rx_pred_=", rxode2::rxFromSE(.prd))
  .r <- get("rx_r_", envir = .s)
  .r <- paste0("rx_r_=", rxode2::rxFromSE(.r))
  .yj <- paste(get("rx_yj_", envir = .s))
  .yj <- paste0("rx_yj_~", rxode2::rxFromSE(.yj))
  .lambda <- paste(get("rx_lambda_", envir = .s))
  .lambda <- paste0("rx_lambda_~", rxode2::rxFromSE(.lambda))
  .hi <- paste(get("rx_hi_", envir = .s))
  .hi <- paste0("rx_hi_~", rxode2::rxFromSE(.hi))
  .low <- paste(get("rx_low_", envir = .s))
  .low <- paste0("rx_low_~", rxode2::rxFromSE(.low))
  .lhs0 <- .s$..lhs0
  if (is.null(.lhs0)) .lhs0 <- ""
  .lhs <- .s$..lhs
  if (is.null(.lhs)) .lhs <- ""
  .ddt <- .s$..ddt
  if (is.null(.ddt)) .ddt <- ""
  .s$..pred <- paste(c(
    .s$..stateInfo["state"],
    .lhs0,
    .ddt,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .r,
    .lhs,
    .s$..stateInfo["statef"],
    .s$..stateInfo["dvid"],
    "tad=tad()",
    "dosenum=dosenum()",
    ""
  ), collapse = "\n")
  .s$..pred.nolhs <- paste(c(
    .s$..stateInfo["state"],
    .lhs0,
    .ddt,
    .yj,
    .lambda,
    .hi,
    .low,
    .prd,
    .r,
    .s$..stateInfo["statef"],
    .s$..stateInfo["dvid"],
    ""
  ), collapse = "\n")
  if (sum.prod) {
    .malert("stabilizing round off errors in predictions or EBE model...")
    .s$..pred <- rxode2::rxSumProdModel(.s$..pred)
    .msuccess("done")
  }
  if (optExpression) {
    .s$..pred <- rxode2::rxOptExpr(.s$..pred, "EBE model")
  }
}

.innerInternal <- function(ui, s) {
  if (exists("..maxTheta", s)) {
    .eventTheta <- rep(0L, s$..maxTheta)
  } else {
    .eventTheta <- integer()
  }
  if (exists("..maxEta", s)) {
    .eventEta <- rep(0L, s$..maxEta)
  } else {
    .eventEta <- integer()
  }
  for (.v in s$..eventVars) {
    .vars <- as.character(get(.v, envir = s))
    .vars <- rxGetModel(paste0("rx_lhs=", .vars))$params
    for (.v2 in .vars) {
      .reg <- rex::rex(start, "ETA_", capture(any_numbers), "_", end)
      if (regexpr(.reg, .v2) != -1) {
        .num <- as.numeric(sub(.reg, "\\1", .v2))
        .eventEta[.num] <- 1L
      }
      .reg <- rex::rex(start, "THETA_", capture(any_numbers), "_", end)
      if (regexpr(.reg, .v2) != -1) {
        .num <- as.numeric(sub(.reg, "\\1", .v2))
        .eventTheta[.num] <- 1L
      }
    }
  }
  pred.opt <- NULL
  inner <- .toRx(s$..inner, "compiling inner model...")
  .sumProd <- rxode2::rxGetControl(ui, "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(ui, "optExpression", TRUE)
  .predMinusDv <- rxode2::rxGetControl(ui, "predMinusDv", TRUE)
  if (!is.null(inner)) {
    if (.sumProd) {
      .malert("stabilizing round off errors in FD model...")
      s$..pred.nolhs <- rxSumProdModel(s$..pred.nolhs)
      .msuccess("done")
    }
    if (.optExpression) {
      s$..pred.nolhs <- rxode2::rxOptExpr(s$..pred.nolhs, "FD model")
    }
    s$..pred.nolhs <- paste(c(
      paste0("params(", paste(inner$params, collapse = ","), ")"),
      s$..pred.nolhs
    ), collapse = "\n")
    pred.opt <- s$..pred.nolhs
  }
  .ret <- list(
    inner = inner,
    predOnly = .toRx(s$..pred, "compiling EBE model..."),
    extra.pars = s$..extraPars,
    outer = .toRx(s$..outer),
    predNoLhs = .toRx(pred.opt, "compiling events FD model..."),
    theta = NULL,
    ## warn=.zeroSens,
    pred.minus.dv = .predMinusDv,
    log.thetas = .nullInt(s$..extraTheta[["exp"]]),
    log.etas = .nullInt(s$..extraEta[["exp"]]),
    extraProps = s$..extraTheta,
    eventTheta = .eventTheta,
    eventEta = .eventEta
    ## ,
    ## cache.file=cache.file
  )
  class(.ret) <- "foceiModelList"
  .ret
}

#' @export
rxUiGet.focei <- function(x, ...) {
  .s <- rxUiGet.foceiEnv(x, ...)
  .innerInternal(x[[1]], .s)
}
#attr(rxUiGet.focei, "desc") <- "Get the FOCEi foceiModelList object"

#' @export
rxUiGet.foce <- function(x, ...) {
  .s <- rxUiGet.foceEnv(x, ...)
  .innerInternal(x[[1]], .s)
}
#attr(rxUiGet.foce, "desc") <- "Get the FOCE foceiModelList object"


#' @export
rxUiGet.ebe <- function(x, ...) {
  .s <- rxUiGet.getEBEEnv(x, ...)
  .innerInternal(x[[1]], .s)
}
#attr(rxUiGet.ebe, "desc") <- "Get the EBE foceiModelList object"

#' @export
rxUiGet.foceiModelDigest <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- get("iniDf", .ui)
  .sumProd <- rxode2::rxGetControl(.ui, "sumProd", FALSE)
  .optExpression <- rxode2::rxGetControl(.ui, "optExpression", TRUE)
  .predMinusDv   <- rxode2::rxGetControl(.ui, "predMinusDv", TRUE)
  digest::digest(c(all(is.na(.iniDf$neta1)),
                   rxode2::rxGetControl(.ui, "interaction", 1L),
                   .sumProd, .optExpression, .predMinusDv,
                   rxode2::rxGetControl(.ui, "addProp", getOption("rxode2.addProp", "combined2")),
                   .ui$lstExpr))
}
#attr(rxUiGet.foceiModelDigest, "desc") <- "Get the md5 digest for the focei model"
#' @export
rxUiGet.foceiModelCache <- function(x, ...) {
  file.path(rxode2::rxTempDir(),
            paste0("focei-", rxUiGet.foceiModelDigest(x, ...), ".qs"))
}
#attr(rxUiGet.foceiModelCache, "desc") <- "Get the focei cache file for a model"

#' @export
rxUiGet.foceiModel <- function(x, ...) {
  .cacheFile <- rxUiGet.foceiModelCache(x, ...)
  if (file.exists(.cacheFile)) {
    return(qs::qread(.cacheFile))
  }
  .ui <- x[[1]]
  .iniDf <- get("iniDf", .ui)
  if (all(is.na(.iniDf$neta1))) {
    .ret <- rxUiGet.ebe(x, ...)
  } else {
    if (rxode2::rxGetControl(.ui, "interaction", 1L)) {
      .ret <- rxUiGet.focei(x, ...)
    } else {
      .ret <- rxUiGet.foce(x, ...)
    }
  }
  qs::qsave(.ret, .cacheFile)
  .ret
}
# attr(rxUiGet.foceiModel, "desc") <- "Get focei model object"

#' @export
rxUiGet.foceiFixed <- function(x, ...) {
  .x <- x[[1]]
  .df <- get("iniDf", .x)
  .dft <- .df[!is.na(.df$ntheta), ]
  .fix <- .dft$fix
  .dft <- .df[is.na(.df$ntheta), ]
  c(.fix, .dft$fix)
}
#attr(rxUiGet.foFixed, "desc") <- "focei theta fixed vector"

#' @export
rxUiGet.foceiEtaNames <- function(x, ...) {
  .x <- x[[1]]
  .df <- get("iniDf", .x)
  .dft <- .df[is.na(.df$ntheta), ]
  .dft[.dft$neta1 == .dft$neta2, "name"]
}
#attr(rxUiGet.foceiEtaNames, "desc") <- "focei eta names"


#' This assigns the tolerances based on a different tolerance for the sensitivity equations
#'
#' It will update and modify the control inside of the UI.
#'
#' It also updates the predNeq that is needed for numeric derivatives
#'
#' @param ui rxode2 UI object
#' @param env focei environment for solving
#' @return Called for side effects
#' @author Matthew L. Fidler
#' @noRd
.foceiOptEnvAssignTol <- function(ui, env) {
  .len <- length(env$model$pred.nolhs$state)
  rxode2::rxAssignControlValue(ui, "predNeq", .len)
  .atol <- rep(rxode2::rxGetControl(ui, "atol", 5e-06), .len)
  .rtol <- rep(rxode2::rxGetControl(ui, "rtol", 5e-06), .len)
  .ssAtol <- rep(rxode2::rxGetControl(ui, "ssAtol", 5e-06), .len)
  .ssRtol <- rep(rxode2::rxGetControl(ui, "ssRtol", 5e-06), .len)
  if (!is.null(env$model$inner)) {
    .len2 <- length(env$model$inner$state) - .len
    .defSens <- 5e-06 / 2
    .atol <- c(.atol, rep(rxode2::rxGetControl(ui, "atolSens", .defSens), .len2))
    .rtol <- c(.rtol, rep(rxode2::rxGetControl(ui, "rtolSens", .defSens), .len2))
    .ssAtol <- c(.ssAtol, rep(rxode2::rxGetControl(ui, "ssAtolSens", .defSens), .len2))
    .ssRtol <- c(.ssAtol, rep(rxode2::rxGetControl(ui, "ssRtolSens", .defSens), .len2))
  }
  rxode2::rxAssignControlValue(ui, "atol", .atol)
  rxode2::rxAssignControlValue(ui, "rtol", .rtol)
  rxode2::rxAssignControlValue(ui, "ssAtol", .ssAtol)
  rxode2::rxAssignControlValue(ui, "ssRtol", .ssRtol)
}

#'  This sets up the initial omega/eta estimates and the boundaries for the whole system
#'
#' @param ui rxode2 UI object
#' @param env focei solving environment
#' @return NoHing, called for side effecs
#' @author Matthew L. Fidler
#' @noRd
.foceiOptEnvSetupBounds <- function(ui, env) {
  .iniDf <- ui$iniDf
  .w <- which(!is.na(.iniDf$ntheta))
  .lower <- .iniDf$lower[.w]
  .upper <- .iniDf$upper[.w]
  env$thetaIni <- ui$theta
  env$thetaIni <- setNames(env$thetaIni, paste0("THETA[", seq_along(env$thetaIni), "]"))
  rxode2::rxAssignControlValue(ui, "nfixed", sum(ui$iniDf$fix))
  .mixed <- !is.null(env$etaNames)
  if (.mixed && length(env$etaNames) == 0L) .mixed <- FALSE
  if (!.mixed) {
    rxode2::rxAssignControlValue(ui, "nomega", 0)
    rxode2::rxAssignControlValue(ui, "neta", 0)
    env$xType <- -1
    rxode2::rxAssignControlValue(ui, "ntheta", length(ui$iniDf$lower))
  } else {
    .om0 <- ui$omega
    .diagXform <- rxode2::rxGetControl(ui, "diagXform", "sqrt")
    env$rxInv <- rxode2::rxSymInvCholCreate(mat = .om0, diag.xform = .diagXform)
    env$xType <- env$rxInv$xType
    .om0a <- .om0
    .om0a <- .om0a / rxode2::rxGetControl(ui, "diagOmegaBoundLower", 100)
    .om0b <- .om0
    .om0b <- .om0b * rxode2::rxGetControl(ui, "diagOmegaBoundUpper", 5)
    .om0a <- rxode2::rxSymInvCholCreate(mat = .om0a, diag.xform = .diagXform)
    .om0b <- rxode2::rxSymInvCholCreate(mat = .om0b, diag.xform = .diagXform)
    .omdf <- data.frame(a = .om0a$theta, m = env$rxInv$theta, b = .om0b$theta, diag = .om0a$theta.diag)
    .omdf$lower <- with(.omdf, ifelse(a > b, b, a))
    .omdf$lower <- with(.omdf, ifelse(lower == m, -Inf, lower))
    .omdf$lower <- with(.omdf, ifelse(!diag, -Inf, lower))
    .omdf$upper <- with(.omdf, ifelse(a < b, b, a))
    .omdf$upper <- with(.omdf, ifelse(upper == m, Inf, upper))
    .omdf$upper <- with(.omdf, ifelse(!diag, Inf, upper))
    rxode2::rxAssignControlValue(ui, "nomega", length(.omdf$lower))
    rxode2::rxAssignControlValue(ui, "neta", sum(.omdf$diag))
    rxode2::rxAssignControlValue(ui, "ntheta", length(.lower))
    .lower <- c(.lower, .omdf$lower)
    .upper <- c(.upper, .omdf$upper)
  }
  env$lower <- .lower
  env$upper <- .upper
  env$etaMat <- rxode2::rxGetControl(ui, "etaMat", NULL)
  env
}

#' Setup the scaleC
#'
#' @param ui rxode2 UI
#' @param env Focei setup environment
#' @return NoHing called for side effects
#' @author Matthew L. Fidler
#' @noRd
.foceiOptEnvSetupScaleC <- function(ui, env) {
  .controlScaleC <- rxode2::rxGetControl(ui, "scaleC", NULL)
  .len <- length(env$lower)
  if (is.null(.controlScaleC)) {
    .scaleC <- rep(NA_real_, .len)
  } else {
    .scaleC <- as.double(.controlScaleC)
  }
  .lenC <- length(.scaleC)
  if (.len > .lenC) {
    .scaleC <- c(.scaleC, rep(NA_real_, .len - .lenC))
  } else if (.len < .lenC) {
    .scaleC <- .scaleC[seq(1, .lenC)]
    warning("scaleC control option has more options than estimated population parameters, please check.")
  }

  .ini <- ui$iniDf
  .ini <- .ini[!is.na(.ini$err), c("est", "err", "ntheta")]
  for (.i in seq_along(.ini$err)) {
    if (is.na(.scaleC[.ini$ntheta[.i]])) {
      if (any(.ini$err[.i] == c("boxCox", "yeoJohnson", "pow2", "tbs", "tbsYj"))) {
        .scaleC[.ini$ntheta[.i]] <- 1
      } else if (any(.ini$err[.i] == c("prop", "add", "norm", "dnorm", "logn", "dlogn", "lnorm", "dlnorm"))) {
        .scaleC[.ini$ntheta[.i]] <- 0.5 * abs(.ini$est[.i])
      }
    }
  }
  .muRefCurEval <- ui$muRefCurEval
  .ini <- ui$iniDf
  for (.i in seq_along(.muRefCurEval$parameter)) {
    .curEval <- .muRefCurEval$curEval[.i]
    .par <- .muRefCurEval$parameter[.i]
    .w <- which(.ini$name == .par)
    if (length(.w) == 1) {
      if (!is.na(.ini$ntheta[.w])) {
        .j <- .ini$ntheta[.w]
        if (is.na(.scaleC[.j])) {
          # These have similar deriavtes on a log scale.
          if (.curEval == "exp") {
            # Hence D(S("log(exp(x))"}, "x")
            .scaleC[.j] <- 1 # log scaled
          } else if (.curEval == "factorial") {
            # Hence 1/D(S("log(factorial(x))"}, "x"):
            .scaleC[.j] <- abs(1 / digamma(.ini$est[.j] + 1))
          } else if (.curEval == "gamma") {
            #1/D(log(gamma(x)), x)
            .scaleC[.j] <- abs(1 / digamma(.ini$est[.j]))
          } else if (.curEval == "log") {
            #1/D(log(log(x)), x)
            .scaleC[.j] <- log(abs(.ini$est[.j])) * abs(.ini$est[.j])
          } else if (.curEval == "logit") {
            # 1/D(log(logit(x, a, b)))
            .a <- .muRefCurEval$low[.i]
            .b <- .muRefCurEval$hi[.i]
            .x <- .ini$est[.j]
            .scaleC[.j] <- -1.0*(-.a + .x)^2*(-1.0 + 1.0*(-.a + .b)/(-.a + .x))*log(abs(-1.0 + 1.0*(-.a + b)/(-.a + .x)))/(-.a + .b)
          } else if (.curEval == "expit") {
            # 1/D(log(expit(x, a, b)))
            .a <- .muRefCurEval$low[.i]
            .b <- .muRefCurEval$hi[.i]
            .x <- .ini$est[.j]
            .scaleC[.j] <- 1.0*exp(.x)*(1.0 + exp(-.x))^2*(.a + 1.0*(-.a + .b)/(1.0 + exp(-.x)))/(-.a + .b)
          } else if (.curEval == "probitInv") {
            .a <- .muRefCurEval$low[.i]
            .b <- .muRefCurEval$hi[.i]
            .x <- .ini$est[.j]
            .scaleC[.j] <- 1.4142135623731*exp(0.5*.x^2)*sqrt(pi)*(.a + 0.5*(-.a + .b)*(1.0 + erf(0.707106781186547*x)))/(-.a + .b)
          } else if (.curEval == "probit") {
            .a <- .muRefCurEval$low[.i]
            .b <- .muRefCurEval$hi[.i]
            .x <- .ini$est[.j]
            erfinvF  <- function(y) {
              if(abs(y) > 1) return(NA_real_)
              sqrt(qchisq(abs(y),1)/2) * sign(y)
            }
            .scaleC[.j] <- sqrt(2)*(-.a+.b)*erfinvF(-1+2*(-.a+.x)/(-.a+.b))/sqrt(pi)/2*exp(((erfinvF(-1+2*(-.a+.x)/(-.a+.b))) ^ 2))
          }
        }
      }
    }
  }
  env$scaleC <- .scaleC
}

#' This sets up the transformation bounds and indexes and bounds for inner.cpp
#'
#' Note that the C code assumes the index starts at 1
#'
#' @param ui rxode2 ui environment
#' @param env focei environment used for solving
#' @return Nothing called for side effects
#' @author Matthew L. Fidler
#' @noRd
.foceiOptEnvSetupTransformIndexs <- function(ui, env) {
  .muRefCurEval <- ui$muRefCurEval
  .ini <- ui$iniDf
  .ini <- .ini[!is.na(.ini$ntheta), c("ntheta", "name")]
  names(.ini)[2] <- "parameter"
  .transform <- merge(.ini, .muRefCurEval)
  .transform <- .transform[order(.transform$ntheta), ]

  env$logThetasF <- .transform[which(.transform$curEval == "exp"), "ntheta"]

  env$logitThetasF <- .transform[which(.transform$curEval == "expit"), "ntheta"]
  env$logitThetasLowF <- .transform[which(.transform$curEval == "expit"), "low"]
  env$logitThetasHiF <- .transform[which(.transform$curEval == "expit"), "hi"]

  env$probitThetasF <- .transform[which(.transform$curEval == "probitInv"), "ntheta"]
  env$probitThetasLowF <- .transform[which(.transform$curEval == "probitInv"), "low"]
  env$probitThetasHiF <- .transform[which(.transform$curEval == "probitInv"), "hi"]
}

# focei.mu.ref
# eta# and the corresponding theta number

#' @export
rxUiGet.foceiMuRefVector <- function(x, ...) {
  .ui <- x[[1]]
  .iniDf <- .ui$iniDf
  .muRefDataFrame <- .ui$muRefDataFrame
  .w <- which(!is.na(.iniDf$ntheta))
  .i2 <- .iniDf[-.w, ]
  if (length(.i2$name) > 0) {
    .i2 <- .i2[.i2$neta1 == .i2$neta2, ]
    .i2 <- .i2[order(.i2$neta1), ]
    vapply(seq_along(.i2$neta1), function(i){
      if (.i2$fix[i]) return(-1L)
      .name <- .i2$name[i]
      .w <- which(.muRefDataFrame$eta == .name)
      if (length(.w) != 1) return(-1L)
      .name <- .muRefDataFrame$theta[.w]
      .w <- which(.iniDf$name == .name)
      if (length(.w) != 1) return(-1L)
      if (.iniDf$fix[.w]) return(-1L)
      .iniDf$ntheta[.w] - 1L
    }, integer(1))
  } else {
    integer(0)
  }
}
#attr(rxUiGet.foceiMuRefVector, "desc") <- "focei mu ref vector"

#' @export
rxUiGet.foceiSkipCov <- function(x, ...) {
  .ui <- x[[1]]
  .maxTheta <- max(.ui$iniDf$ntheta, na.rm=TRUE)
  .theta <- .ui$iniDf[!is.na(.ui$iniDf$ntheta), ]
  .skipCov <- rep(FALSE, .maxTheta)
  .skipCov[which(!is.na(.theta$err))] <- TRUE
  .skipCov[.theta$fix] <- TRUE
  .skipCov
}
#attr(rxUiGet.foceiSkipCov, "desc") <- "what covariance elements to skip"

#'  Setup the skip covariate function
#'
#'
#' @param ui rxode2 parsed function
#' @param env environment
#' @return Nothing called for side effects.
#' @author Matthew L. Fidler
#' @noRd
.foceiSetupSkipCov <- function(ui, env) {
  env$skipCov <- rxode2::rxGetControl(ui, "skipCov", NULL)
  if (is.null(env$skipCov)) {
    env$skipCov <- ui$foceiSkipCov
  }
  .maxTheta <- max(ui$iniDf$ntheta, na.rm=TRUE)
  if (length(env$skipCov) > .maxTheta) {
    if (all(env$skipCov[-seq_len(.maxTheta)])) {
      assign("skipCov",env$skipCov[seq_len(.maxTheta)], env)
    }
  }
  if (length(env$skipCov) != .maxTheta) {
    stop("'skipCov' improperly specified", call.=FALSE)
  }
}

.foceiOptEnvLik <- function(ui, env) {
  #if (!exists("noLik", envir = env)){
  if (!exists("model", envir=env)) {
    env$model <- rxUiGet.foceiModel(list(ui))
  }
  #} else {
    #env$model <- rxUiGet.ebe(list(ui))
  #}
  .foceiOptEnvAssignTol(ui, env)
  .foceiOptEnvSetupBounds(ui, env)
  .foceiOptEnvSetupScaleC(ui, env)
  .foceiOptEnvSetupTransformIndexs(ui, env)
  .foceiSetupSkipCov(ui, env)
  env$control <- get("control", envir=ui)
  env$control$nF <- 0
  env$control$printTop <- TRUE
  env
}

#' @export
rxUiGet.foceiOptEnv <- function(x, ...) {
  .x <- x[[1]]
  if (exists("foceiEnv", envir=.x)) {
    .env <- get("foceiEnv", envir=.x)
    rm("foceiEnv", envir=.x)
  } else {
    .env <- new.env(parent=emptyenv())
  }
  .env$etaNames <- rxUiGet.foceiEtaNames(x, ...)
  .env$thetaFixed <- rxUiGet.foceiFixed(x, ...)
  rxode2::rxAssignControlValue(.x, "foceiMuRef", .x$foceiMuRefVector)
  .env$adjLik <- rxode2::rxGetControl(.x, "adjLik", TRUE)
  .env$diagXformInv <- c("sqrt" = ".square", "log" = "exp", "identity" = "identity")[rxode2::rxGetControl(.x, "diagXform", "sqrt")]
  .env$thetaNames <- .x$iniDf[!is.na(.x$iniDf$ntheta), "name"]
  # FIXME is ODEmodel needed?
  .env$ODEmodel <- TRUE
  .foceiOptEnvLik(.x, .env)
  .env
}
attr(rxUiGet.foceiOptEnv, "desc") <- "Get focei optimization environment"
#' This function process the data for use in focei
#'
#' The $origData is the data that is fed into the focei before modification
#' The $dataSav is the data saved for focei
#'
#' @param data Input dataset
#' @param env focei environment where focei family is run
#' @param ui rxode2 ui
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.foceiPreProcessData <- function(data, env, ui) {
  env$origData <- data
  .covNames <- ui$covariates
  colnames(data) <- vapply(names(data), function(x) {
      if (any(x == .covNames)) {
        return(x)
      } else {
        return(toupper(x))
      }
  }, character(1))
  if (is.null(data$ID)) stop('"ID" not found in data')
  if (is.null(data$DV)) stop('"DV" not found in data')
  if (is.null(data$EVID)) data$EVID <- 0
  if (is.null(data$AMT)) data$AMT <- 0
  ## Make sure they are all double amounts.
  for (.v in c("TIME", "AMT", "DV", .covNames)) {
    if (!any(names(data) == .v)) {
      stop("missing '", .v, "' in data", call.=FALSE)
    }
    data[[.v]] <- as.double(data[[.v]])
  }
  env$dataSav <- as.data.frame(rxode2::etTrans(inData=data, obj=ui$mv0,
                                               addCmt=TRUE, dropUnits=TRUE,
                                               allTimeVar=TRUE, keepDosingOnly=FALSE))
}

.thetaReset <- new.env(parent = emptyenv())
#' Internal focei fit function in R
#'
#' @param .ret Internal focei environment
#' @return Modified focei environment with fit information (from C++)
#' @author Matthew L. Fidler
#' @noRd
.foceiFitInternal <- function(.ret) {
  this.env <- new.env(parent=emptyenv())
  assign("err", "theta reset", this.env)
  .thetaReset$thetaNames <- .ret$thetaNames
  while (this.env$err == "theta reset") {
    assign("err", "", this.env)
    .ret0 <- tryCatch(
    {
      foceiFitCpp_(.ret)
    },
    error = function(e) {
      if (regexpr("theta reset", e$message) != -1) {
        assign("zeroOuter", FALSE, this.env)
        assign("zeroGrad", FALSE, this.env)
        if (regexpr("theta reset0", e$message) != -1) {
          assign("zeroGrad", TRUE, this.env)
        }  else if (regexpr("theta resetZ", e$message) != -1) {
          assign("zeroOuter", TRUE, this.env)
        }
        assign("err", "theta reset", this.env)
      } else {
        assign("err", e$message, this.env)
      }
    })
    if (this.env$err == "theta reset") {
      .nm <- names(.ret$thetaIni)
      .ret$thetaIni <- setNames(.thetaReset$thetaIni + 0.0, .nm)
      .ret$rxInv$theta <- .thetaReset$omegaTheta
      .ret$control$printTop <- FALSE
      .ret$etaMat <- .thetaReset$etaMat
      .ret$control$etaMat <- .thetaReset$etaMat
      .ret$control$maxInnerIterations <- .thetaReset$maxInnerIterations
      .ret$control$nF <- .thetaReset$nF
      .ret$control$gillRetC <- .thetaReset$gillRetC
      .ret$control$gillRet <- .thetaReset$gillRet
      .ret$control$gillRet <- .thetaReset$gillRet
      .ret$control$gillDf <- .thetaReset$gillDf
      .ret$control$gillDf2 <- .thetaReset$gillDf2
      .ret$control$gillErr <- .thetaReset$gillErr
      .ret$control$rEps <- .thetaReset$rEps
      .ret$control$aEps <- .thetaReset$aEps
      .ret$control$rEpsC <- .thetaReset$rEpsC
      .ret$control$aEpsC <- .thetaReset$aEpsC
      .ret$control$c1 <- .thetaReset$c1
      .ret$control$c2 <- .thetaReset$c2
      if (this.env$zeroOuter) {
        message("Posthoc reset")
        .ret$control$maxOuterIterations <- 0L
      } else if (this.env$zeroGrad) {
        message("Theta reset (zero gradient values); Switch to bobyqa")
        rxode2::rxReq("minqa")
        .ret$control$outerOptFun <- .bobyqa
        .ret$control$outerOpt <- -1L
      } else {
        message("Theta reset (ETA drift)")
      }
    }
    if (this.env$err != "") {
      stop(this.env$err)
    } else {
      return(.ret0)
    }
  }
}
#'  Restart the estimation if it wasn't successful by moving the parameters (randomly)
#'
#' @param .ret0 Fit
#' @param .ret Input focei environment
#' @param control Control represents the foceiControl to restart the fit
#' @return final focei fit, may still not work
#' @author Matthew L. Fidler
#' @noRd
.nlmixrFoceiRestartIfNeeded <- function(.ret0, .ret, control) {
  .n <- 1
  .est0 <- .ret$thetaIni
  lower <- .ret$lower
  upper <- .ret$upper

  while (inherits(.ret0, "try-error") && control$maxOuterIterations != 0 && .n <= control$nRetries) {
    ## Maybe change scale?
    message(sprintf("Restart %s", .n))
    .ret$control$nF <- 0
    .estNew <- .est0 + 0.2 * .n * abs(.est0) * stats::runif(length(.est0)) - 0.1 * .n
    .estNew <- sapply(
      seq_along(.est0),
      function(.i) {
        if (.ret$thetaFixed[.i]) {
          return(.est0[.i])
        } else if (.estNew[.i] < lower[.i]) {
          return(lower + (.Machine$double.eps)^(1 / 7))
        } else if (.estNew[.i] > upper[.i]) {
          return(upper - (.Machine$double.eps)^(1 / 7))
        } else {
          return(.estNew[.i])
        }
      }
    )
    .ret$thetaIni <- .estNew
    .ret0 <- try(.foceiFitInternal(.ret))
    .n <- .n + 1
  }
  .ret0
}

.foceiFamilyControl <- function(env, ...) {
  .ui <- env$ui
  .control <- env$control
  if (is.null(.control)) {
    .control <- foceiControl()
  }
  if (!inherits(.control, "foceiControl")){
    .control <- do.call(nlmixr2::foceiControl, .control)
  }
  assign("control", .control, envir=.ui)
}

#' Get the cmt() and dvid() lines
#'
#' @param ui rxode UI
#' @return cmt() and dvid() string
#' @author Matthew L. Fidler
#' @noRd
.foceiToCmtLinesAndDvid <- function(ui) {
  .cmtLines <- ui$cmtLines
  paste(c("", vapply(seq_along(.cmtLines),
         function(i){deparse1(.cmtLines[[i]])},
         character(1), USE.NAMES=FALSE),
         deparse1(ui$dvidLine)),
        collapse="\n")
}

#' Setup the par history information
#'
#' @param .ret Return data
#' @return Nothing called for side effects
#' @author Matthew L. Fidler
#' @noRd
.foceiSetupParHistData <- function(.ret) {
  if (exists("parHistData", envir=.ret)) {
    .tmp <- .ret$parHistData
    .tmp <- .tmp[.tmp$type == "Unscaled", names(.tmp) != "type"]
    .iter <- .tmp$iter
    .tmp <- .tmp[, names(.tmp) != "iter"]
    ## .ret$parHistStacked <- data.frame(stack(.tmp), iter = .iter)
    ## names(.ret$parHistStacked) <- c("val", "par", "iter")
    .ret$parHist <- data.frame(iter = .iter, .tmp)
  }
}


.foceiFamilyReturn0 <- function(env, ui, ..., method=NULL, est="none") {
  assignInMyNamespace(".toRxParam", paste0(.uiGetThetaEtaParams(ui, TRUE), "\n",
                                           ui$foceiCmtPreModel, "\n"))
  assignInMyNamespace(".toRxDvidCmt", .foceiToCmtLinesAndDvid(ui))
  .control <- ui$control
  .env <- ui$foceiOptEnv
  .data <- env$data
  .foceiPreProcessData(.data, .env, ui)
  .ret0 <- try(.foceiFitInternal(.env))
  .ret0 <- .nlmixrFoceiRestartIfNeeded(.ret0, .env, .control)
  if (inherits(.ret0, "try-error")) {
    stop("Could not fit data\n  ", attr(.ret0, "condition")$message, call.=FALSE)
  }
  .ret <- .ret0
  if (!is.null(method))
    .ret$method <- method
  .ret$ui <- ui
  .foceiSetupParHistData(.ret)
  if (!all(is.na(ui$iniDf$neta1))) {
    .etas <- .ret$ranef
    .thetas <- .ret$fixef
    .pars <- .Call(`_nlmixr2_nlmixr2Parameters`, .thetas, .etas)
    .ret$shrink <- .Call(`_nlmixr2_calcShrinkOnly`, .ret$omega, .pars$eta.lst, length(.etas$ID))
  }
  .updateParFixed(.ret)
  if (!exists("table", .ret)) {
    .ret$table <- tableControl()
  }
  .nlmixr2FitUpdateParams(.ret)
  .ret$IDlabel <- rxode2::.getLastIdLvl()
  if (exists("skipTable", envir=.ret)) {
    if (.ret$skipTable) {
      .control$calcTables <- FALSE
    }
  }
  assign("est", est, envir=.ret)
  assign("skipCov", .env$skipCov, envir=.ret)
  nmObjHandleModelObject(.ret$model, .ret)
  nmObjHandleControlObject(get("control", envir=.ret), .ret)
  if (.control$calcTables) {
    .ret <- addTable(.ret, updateObject="no", keep=.ret$table$keep, drop=.ret$table$drop,
                     table=.ret$table)
  }
  if (exists("saem", .env)) {
    .saem <- get("saem", envir=.env)
    .saemCfg <- attr(.saem, "saem.cfg")
    # Delete unneeded variables
    .saemCfg2 <- list()
    for (.v in c("i1", "nphi1", "nphi0", "N", "ntotal", "ix_endpnt", "y", "nmc", "niter", "opt", "inits", "Mcovariables")) {
      .saemCfg2[[.v]] <- .saemCfg[[.v]]
    }
    attr(.saem, "saem.cfg") <- .saemCfg2
    rm(list="saem", envir=.env)
    .env$saem0 <- .saem
  }
  for (.item in c("origData", "phiM", "parHist", "saem0")) {
    if (exists(.item, .env)) {
      .obj <- get(.item, envir=.env)
      .size <- object.size(.obj)
      .objC <- qs::qserialize(.obj)
      .size2 <- object.size(.objC)
      if (.size2 < .size) {
        .size0 <- (.size - .size2)
        .malert("compress {  .item } in nlmixr2 object, save { .size0 }" )
        assign(.item, .objC, envir=.env)
      }
    }
  }
  for (.item in c("adj", "adjLik", "diagXformInv", "etaMat", "etaNames",
                  "fullTheta", "scaleC", "gillRet", "gillRetC",
                  "logitThetasF", "logitThetasHiF", "logitThetasLowF", "logThetasF",
                  "lower", "noLik", "objf", "OBJF", "probitThetasF", "probitThetasHiF", "probitThetasLowF",
                  "rxInv", "scaleC", "se", "skipCov", "thetaFixed", "thetaIni", "thetaNames", "upper",
                  "xType", "IDlabel", "ODEmodel",
                  # times
                  "optimTime", "setupTime", "covTime",
                  "parHistData", "dataSav", "theta")) {
    if (exists(.item, .env)) {
      rm(list=.item, envir=.env)
    }
  }
  .ret
}

.foceiFamilyReturn <- function(env, ui, ..., method=NULL, est="none") {
  .envReset <- new.env(parent=emptyenv())
  .envReset$reset <- TRUE
  .envReset$env <- new.env(parent=emptyenv())
  .envReset$ui <- ui
  .envReset$method <- method
  .envReset$est <- est
  lapply(ls(envir = env, all.names = TRUE), function(item) {
    assign(item, get(item, envir = env), envir = .envReset$env)
  })
  .envReset$cacheReset <- FALSE
  while (.envReset$reset) {
    .envReset$reset <- FALSE
    tryCatch({
      .envReset$ret <- .foceiFamilyReturn0(env=env, ui=.envReset$ui, method=.envReset$method, est=.envReset$est)
    },
    error=function(e) {
      if (regexpr("not provided by package", e$message) != -1) {
        if (.envReset$cacheReset) {
          .malert("unsuccessful cache reset; try manual reset with 'rxode2::rxClean()'")
          stop(e)
        } else {
          # reset
          rm(list=ls(envir = env, all.names = TRUE), envir=env)
          lapply(ls(envir = .envReset$env, all.names = TRUE), function(item) {
            assign(item, get(item, envir = .envReset$env), envir = env)
          })
          gc()
          .minfo("try resetting cache")
          rxode2::rxClean()
          .envReset$cacheReset <- TRUE
          .envReset$reset <- TRUE
          .msuccess("done")
        }
      } else {
        stop(e)
      }
    })
  }
  .envReset$ret
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.focei <- function(env, ...) {
  .ui <- env$ui
  .foceiFamilyControl(env, ...)
  on.exit({rm("control", envir=.ui)})
  .foceiFamilyReturn(env, .ui, ..., est="focei")
}


#'@rdname nlmixr2Est
#'@export
nlmixr2Est.foce <- function(env, ...) {
  .ui <- env$ui
  .foceiFamilyControl(env, ...)
  rxode2::rxAssignControlValue(.ui, "interaction", 0L)
  on.exit({rm("control", envir=.ui)})
  env$est <- "foce"
  .foceiFamilyReturn(env, .ui, ..., est="focei")
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.posthoc <- function(env, ...) {
  .ui <- env$ui
  .foceiFamilyControl(env, ...)
  rxode2::rxAssignControlValue(.ui, "interaction", 0L)
  rxode2::rxAssignControlValue(.ui, "covMethod", 0L)
  rxode2::rxAssignControlValue(.ui, "maxOuterIterations", 0L)
  on.exit({rm("control", envir=.ui)})
  env$est <- "posthoc"
  .foceiFamilyReturn(env, .ui, ..., est="posthoc")
}

#' Add objective function line to the return object
#'
#' @param ret Return object
#' @param objDf Objective function data frame to add
#' @return Nothing, called for side effects
#' @author Matthew L. Fidler
#' @noRd
.addObjDfToReturn <- function(ret, objDf) {
  if (inherits(ret, "nlmixr2FitData")) {
    ret <- attr(class(ret), ".foceiEnv")
  }
  .objDf1 <- get("objDf", ret)
  if (any(names(.objDf1) == "Condition Number")) {
    if (!any(names(objDf) == "Condition Number")) {
      objDf[["Condition Number"]] <- NA_real_
    }
  } else if (any(names(objDf) == "Condition Number")) {
    if (!any(names(.objDf1) == "Condition Number")) {
      .objDf1[["Condition Number"]] <- NA_real_
    }
  }
  assign("objDf", rbind(.objDf1, objDf), envir=ret)
}

#'@rdname nlmixr2Est
#'@export
nlmixr2Est.foi <- function(env, ...) {
  .ui <- env$ui
  .foceiFamilyControl(env, ...)
  .control <- .ui$control
  rxode2::rxAssignControlValue(.ui, "interaction", 0L)
  rxode2::rxAssignControlValue(.ui, "covMethod", 0L)
  rxode2::rxAssignControlValue(.ui, "fo", TRUE)
  rxode2::rxAssignControlValue(.ui, "boundTol", 0)
  on.exit({rm("control", envir=.ui)})
  env$skipTable <- TRUE
  .ret <- .foceiFamilyReturn(env, .ui, ...)
  .objDf <- .ret$objDf
  .ui <- .ret$ui
  .foceiFamilyControl(env, ...)
  rxode2::rxAssignControlValue(.ui, "interaction", 1L)
  rxode2::rxAssignControlValue(.ui, "maxOuterIterations", 0L)
  .ret <- .foceiFamilyReturn(env, .ui, ..., method="FO", est="foi")
  .addObjDfToReturn(.ret, .objDf)
  .ret
}


#'@rdname nlmixr2Est
#'@export
nlmixr2Est.fo <- function(env, ...) {
  .ui <- env$ui
  .foceiFamilyControl(env, ...)
  .control <- .ui$control
  rxode2::rxAssignControlValue(.ui, "interaction", 0L)
  rxode2::rxAssignControlValue(.ui, "covMethod", 0L)
  rxode2::rxAssignControlValue(.ui, "fo", TRUE)
  rxode2::rxAssignControlValue(.ui, "boundTol", 0)
  on.exit({rm("control", envir=.ui)})
  env$skipTable <- TRUE
  .ret <- .foceiFamilyReturn(env, .ui, ...)
  .objDf <- .ret$objDf
  .ui <- .ret$ui
  .foceiFamilyControl(env, ...)
  rxode2::rxAssignControlValue(.ui, "interaction", 0L)
  rxode2::rxAssignControlValue(.ui, "maxOuterIterations", 0L)
  .ret <- .foceiFamilyReturn(env, .ui, ..., method="FO", est="fo")
  .addObjDfToReturn(.ret, .objDf)
  .ret
}
#'@rdname nlmixr2Est
#'@export
nlmixr2Est.output <- function(env, ...) {
  .ui <- env$ui
  .foceiFamilyControl(env, ...)
  rxode2::rxAssignControlValue(.ui, "interaction", 0L)
  rxode2::rxAssignControlValue(.ui, "maxOuterIterations", 0L)
  rxode2::rxAssignControlValue(.ui, "maxInnerIterations", 0L)
  on.exit({rm("control", envir=.ui)})
  if (!exists("est", envir=env)) env$est <- "posthoc"
  .foceiFamilyReturn(env, .ui, ..., est=env$est)
}


nlmixr2CreateOutputFromUi <- function(ui, data=NULL, control=NULL, table=NULL, env=NULL, est="none") {
  if (inherits(ui, "function")) {
    ui <- rxode2::rxode2(ui)
  }
  if (!inherits(ui, "rxUi")) {
    stop("the first argument needs to be from rxode2 ui", call.=FALSE)
  }
  if (inherits(env, "environment")) {
    assign("foceiEnv", env, envir=ui)
  }
  if (!inherits(data, "data.frame")) {
    stop("the 'data' argument must be a data.frame", call.=FALSE)
  }

  .env <- new.env(parent=emptyenv())
  .env$ui <- ui
  .env$data <- data
  .env$control <- control
  .env$table <- table
  .env$est <- est
  class(.env) <- c("output", "nlmixr2Est")
  nlmixr2Est(.env)
}

