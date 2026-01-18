.foceiControlInternal <- c("genRxControl", "resetEtaSize",
                           "resetThetaSize", "resetThetaFinalSize",
                           "outerOptFun", "outerOptTxt", "skipCov",
                           "foceiMuRef", "predNeq", "nfixed", "nomega",
                           "neta", "ntheta", "nF", "printTop", "needOptimHess")

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
#' }
#'
#' @param sigdigTable Significant digits in the final output table.
#'   If not specified, then it matches the significant digits in the
#'   `sigdig` optimization algorithm.  If `sigdig` is NULL, use 3.
#'
#' @param epsilon Precision of estimate for n1qn1 optimization.
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
#'     value.  By default this is 0 (meaning do not scale)
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
#' @param hessEps is a double value representing the epsilon for the
#'   Hessian calculation. This is used for the R matrix calculation.
#'
#' @param hessEpsLlik is a double value representing the epsilon for
#'   the Hessian calculation when doing focei generalized
#'   log-likelihood estimation.  This is used for the R matrix
#'   calculation.
#'
#' @param optimHessType The hessian type for when calculating the
#'   individual hessian by numeric differences (in generalized
#'   log-likelihood estimation).  The options are "central", and
#'   "forward".  The central differences is what R's `optimHess()`
#'   uses and is the default for this method. (Though the "forward" is
#'   faster and still reasonable for most cases).  The Shi21 cannot be
#'   changed for the Gill83 algorithm with the optimHess in a
#'   generalized likelihood problem.
#'
#' @param optimHessCovType The hessian type for when calculating the
#'   individual hessian by numeric differences (in generalized
#'   log-likelihood estimation).  The options are "central", and
#'   "forward".  The central differences is what R's `optimHess()`
#'   uses.  While this takes longer in optimization, it is more
#'   accurate, so for calculating the covariance and final likelihood,
#'   the central differences are used. This also uses the modified
#'   Shi21 method
#'
#' @param shi21maxOuter The maximum number of steps for the
#'   optimization of the forward-difference step size.  When not zero,
#'   use this instead of Gill differences.
#'
#' @param shi21maxInner The maximum number of steps for the
#'   optimization of the individual Hessian matrices in the
#'   generalized likelihood problem. When 0, un-optimized finite differences
#'   are used.
#'
#' @param shi21maxInnerCov The maximum number of steps for the
#'   optimization of the individual Hessian matrices in the
#'   generalized likelihood problem for the covariance step. When 0,
#'   un-optimized finite differences are used.
#'
#' @param shi21maxFD The maximum number of steps for the optimization
#'   of the forward difference step size when using dosing events (lag
#'   time, modeled duration/rate and bioavailability)
#'
#' @param centralDerivEps Central difference tolerances.  This is a
#'   numeric vector of relative difference and absolute difference.
#'   The central/forward difference step size h is calculated as:
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
#'
#' @param iovXform This is the transformation used on the diagonal
#'     of the `iov`. The possibilities are:
#'
#' \itemize{
#'
#'  \item \code{sd} Estimate the IOV as the standard deviation for IOV
#'
#'  \item \code{var} Estimate the IOV as the variance for IOV.
#'
#'  \item \code{logsd} Estimate the IOV as the log(sd) instead of sd.
#'
#'  \item \code{logvar} Estimate the IOV as the log(var) instead of variance.
#'
#' }
#'
#' @param sumProd Is a boolean indicating if the model should change
#'     multiplication to high precision multiplication and sums to
#'     high precision sums using the PreciseSums package.  By default
#'     this is \code{FALSE}.
#'
#'
#' @param optExpression Optimize the rxode2 expression to speed up
#'     calculation. By default this is turned on.
#'
#' @param literalFix boolean, substitute fixed population values as
#'   literals and re-adjust ui and parameter estimates after
#'   optimization; Default is `TRUE`.
#'
#' @param literalFixRes boolean, substitute fixed population values as
#'   literals and re-adjust ui and parameter estimates after
#'   optimization; Default is `TRUE`.
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
#'     adjusted. (bobyqa)
#'
#' @param rhoend The smallest value of the trust region radius that
#'     is allowed. If not defined, then 10^(-sigdig-1) will be used. (bobyqa)
#'
#' @param npt The number of points used to approximate the objective
#'     function via a quadratic approximation for bobyqa. The value
#'     of npt must be in the interval [n+2,(n+1)(n+2)/2] where n is
#'     the number of parameters in par. Choices that exceed 2*n+1 are
#'     not recommended. If not defined, it will be set to 2*n + 1. (bobyqa)
#'
#' @param eval.max Number of maximum evaluations of the objective function (nlmimb)
#'
#' @param iter.max Maximum number of iterations allowed (nlmimb)
#'
#' @param rel.tol Relative tolerance before nlminb stops (nlmimb).
#'
#' @param x.tol X tolerance for nlmixr2 optimizer
#'
#' @param abstol Absolute tolerance for nlmixr2 optimizer (BFGS)
#'
#' @param reltol  tolerance for nlmixr2 (BFGS)
#'
#' @param gillK The total number of possible steps to determine the
#'     optimal forward/central difference step size per parameter (by
#'     the Gill 1983 method).  If 0, no optimal step size is
#'     determined.  Otherwise this is the optimal step size
#'     determined.
#'
#' @param gillKcovLlik The total number of possible steps to determine
#'   the optimal forward/central difference step per parameter when
#'   using the generalized focei log-likelihood method (by the Gill
#'   1986 method).  If 0, no optimal step size is
#'   determined. Otherwise this is the optimal step size is determined
#'
#' @param gillRtol The relative tolerance used for Gill 1983
#'     determination of optimal step size.
#'
#' @param scaleType The scaling scheme for nlmixr2.  The supported types are:
#'
#' \itemize{
#' \item \code{nlmixr2}  In this approach the scaling is performed by the following equation:
#'
#'    \deqn{v_{scaled}}{Vscaled} = (\deqn{v_{current} - v_{init}}{Vcurrent - Vinit})*scaleC[i] + scaleTo
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
#'   \deqn{v_{scaled}}{Vscaled} = \deqn{v_{current}}{Vcurrent}/\deqn{v_{init}}{Vinit}*scaleTo
#'
#' \item \code{multAdd} This approach changes the scaling based on
#' the parameter being specified.  If a parameter is defined in an
#' exponential block (ie exp(theta)), then it is scaled on a
#' linearly, that is:
#'
#'   \deqn{v_{scaled}}{Vscaled} = (\deqn{v_{current}-v_{init}}{Vcurrent-Vinit}) + scaleTo
#'
#' Otherwise the parameter is scaled multiplicatively.
#'
#'    \deqn{v_{scaled}}{Vscaled} = \deqn{v_{current}}{Vcurrent}/\deqn{v_{init}}{Vinit}*scaleTo
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
#'    \item Factorials are scaled by abs(1/digamma(initial_estimate+1))
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
#'     \deqn{v_{scaled}}{Vscaled} = (\deqn{v_{unscaled}-C_{1}}{Vunscaled-C1})/\deqn{C_{2}}{C2}
#'
#'
#'     Where
#'
#'
#'     The other data normalization approaches follow the following formula
#'
#'     \deqn{v_{scaled}}{Vscaled} = (\deqn{v_{unscaled}-C_{1}}{Vunscaled-C1})/\deqn{C_{2}}{C2}
#'
#' \itemize{
#'
#' \item \code{rescale2} This scales all parameters from (-1 to 1).
#'     The relative differences between the parameters are preserved
#'     with this approach and the constants are:
#'
#'     \deqn{C_{1}}{C1} = (max(all unscaled values)+min(all unscaled values))/2
#'
#'     \deqn{C_{2}}{C2} = (max(all unscaled values) - min(all unscaled values))/2
#'
#'
#' \item \code{rescale} or min-max normalization. This rescales all
#'     parameters from (0 to 1).  As in the \code{rescale2} the
#'     relative differences are preserved.  In this approach:
#'
#'     \deqn{C_{1}}{C1} = min(all unscaled values)
#'
#'     \deqn{C_{2}}{C2} = max(all unscaled values) - min(all unscaled values)
#'
#'
#' \item \code{mean} or mean normalization.  This rescales to center
#'     the parameters around the mean but the parameters are from 0
#'     to 1.  In this approach:
#'
#'     \deqn{C_{1}}{C1} = mean(all unscaled values)
#'
#'     \deqn{C_{2}}{C2} = max(all unscaled values) - min(all unscaled values)
#'
#' \item \code{std} or standardization.  This standardizes by the mean
#'      and standard deviation.  In this approach:
#'
#'     \deqn{C_{1}}{C1} = mean(all unscaled values)
#'
#'     \deqn{C_{2}}{C2} = sd(all unscaled values)
#'
#' \item \code{len} or unit length scaling.  This scales the
#'    parameters to the unit length.  For this approach we use the Euclidean length, that
#'    is:
#'
#'     \deqn{C_{1}}{C1} = 0
#'
#'     \deqn{C_{2}}{C2} = \deqn{\sqrt(v_1^2 + v_2^2 + \cdots + v_n^2)}{sqrt(v_1^2 + v_2^2 + ... + v_n^2)}
#'
#'
#' \item \code{constant} which does not perform data normalization. That is
#'
#'     \deqn{C_{1}}{C1} = 0
#'
#'     \deqn{C_{2}}{C2} = 1
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
#' @param gillStepCovLlik Same as above but during generalized focei
#'   log-likelihood
#'
#' @param gillFtolCov The gillFtol is the gradient error tolerance
#'     that is acceptable before issuing a warning/error about the
#'     gradient estimates during the covariance step.
#'
#' @param gillFtolCovLlik Same as above but applied during generalized
#'   log-likelihood estimation.
#'
#' @param rmatNorm A parameter to normalize gradient step size by the
#'     parameter value during the calculation of the R matrix
#'
#' @param rmatNormLlik A parameter to normalize gradient step size by
#'   the parameter value during the calculation of the R matrix if you
#'   are using generalized log-likelihood Hessian matrix.
#'
#' @param smatNorm A parameter to normalize gradient step size by the
#'     parameter value during the calculation of the S matrix
#'
#' @param smatNormLlik A parameter to normalize gradient step size by
#'   the parameter value during the calculation of the S matrix if you
#'   are using the generalized log-likelihood.
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
#' @param eventType Event gradient type for dosing events; Can be
#'   "central" or "forward"
#'
#' @param gradProgressOfvTime This is the time for a single objective
#'     function evaluation (in seconds) to start progress bars on gradient evaluations
#'
#' @param badSolveObjfAdj The objective function adjustment when the
#'   ODE system cannot be solved.  It is based on each individual bad
#'   solve.
#'
#' @param compress Should the object have compressed items
#'
#' @param etaMat Eta matrix for initial estimates or final estimates
#'   of the ETAs.
#'
#'   This can also be a fit to take use the final estimation estimates
#'   and use them as the initial eta value of the next fit.
#'
#'   By default, it will be the initial values of the etas from the
#'   last fit (if supplied) or missing, meaning all ETAs start at
#'   zero (`NULL`)
#'
#'   When this value is `NA`, the initial ETA estimates are not taken
#'   from the last fit.
#'
#' @param addProp specifies the type of additive plus proportional
#'   errors, the one where standard deviations add (combined1) or the
#'   type where the variances add (combined2).
#'
#' The combined1 error type can be described by the following equation:
#'
#'   \deqn{y = f + (a + b\times f^c) \times \varepsilon}{y = f + (a + b*f^c)*err}
#'
#' The combined2 error model can be described by the following equation:
#'
#'  \deqn{y = f + \sqrt{a^2 + b^2\times f^{2\times c}} \times \varepsilon}{y = f + sqrt(a^2 + b^2*(f^c)^2)*err}
#'
#'  Where:
#'
#'  - y represents the observed value
#'
#'  - f represents the predicted value
#'
#'  - a  is the additive standard deviation
#'
#'  - b is the proportional/power standard deviation
#'
#'  - c is the power exponent (in the proportional case c=1)
#'
#' @param odeRecalcFactor The ODE recalculation factor when ODE
#'   solving goes bad, this is the factor the rtol/atol is reduced
#'
#' @param rxControl `rxode2` ODE solving options during fitting, created with `rxControl()`
#'
#' @param fallbackFD Fallback to the finite differences if the
#'   sensitivity equations do not solve.
#'
#' @param smatPer A percentage representing the number of failed
#'   parameter gradients for each individual (which are replaced with
#'   the overall gradient for the parameter) out of the total number
#'   of gradients parameters (ie `ntheta*nsub`) before the S matrix is
#'   considered to be a bad matrix.
#'
#' @param sdLowerFact A factor for multiplying the estimate by when
#'   the lower estimate is zero and the error is known to represent a
#'   standard deviation of a parameter (like add.sd, prop.sd, pow.sd,
#'   lnorm.sd, etc).  When zero, no factor is applied.  If your
#'   initial estimate is 0.15 and your lower bound is zero, then the
#'   lower bound would be assumed to be 0.00015.
#'
#' @param zeroGradFirstReset boolean, when `TRUE` if the first
#'   gradient is zero, reset the zero gradient to
#'   `sqrt(.Machine$double.eps)` to get past the bad initial estimate,
#'   otherwise error (and possibly reset), when `FALSE` error when the
#'   first gradient is zero.  When `NA` on the last reset, have the
#'   zero gradient ignored, otherwise error and look for another
#'   value.  Default is `TRUE`
#'
#' @param zeroGradRunReset boolean, when `TRUE` if a gradient is zero,
#'   reset the zero gradient to `sqrt(.Machine$double.eps)` to get
#'   past the bad estimate while running.  Otherwise error (and
#'   possibly reset). Default is `TRUE`
#'
#' @param zeroGradBobyqa boolean, when `TRUE` if a gradient is zero,
#'   the reset will change the method to the gradient free bobyqa
#'   method. When `NA`, the zero gradient will change to bobyqa only
#'   when the first gradient is zero.  Default is `TRUE`
#'
#' @param mceta Integer indicating the type of Monte Carlo sampling to
#'   perform for the best initial ETA estimate (based on
#'   `omega`). When:
#'
#'   - `-1` the last eta is used for the optimization (default)
#'
#'   - `0` eta=0 is used for each inner optimization
#'
#'  For the rest of the `mceta`, each parameter's inner objective
#'  function is calculated and the eta set with the best objective
#'  function is used.  With these further options:
#'
#'   - `1` the last eta and eta=0 are used
#'
#'   - `2` the last eta and eta=0 are used, as well as 1 sampled eta
#'   from the omega matrix
#'
#'   - `n` the last eta and eta=0 are used, as well as n-1 sampled
#'   etas from the omega matrix
#'
#' @param nAGQ Number of Gauss-Hermite Adaptive Quadrature points to
#'   take.  When `nAGQ=0`, the AGQ is not used.  With `nAGQ=1`, this
#'   is equivalent to the Laplace method. The adaptive quadrature
#'   expands every node for each of the ETAs, so it can be quite
#'   expensive with a large amount of ETAs.  Once the EBE is obtained
#'   for a subject, you will have nAGQ^neta additional function
#'   evaluations for even nAGQ numbers and (nAGQ^neta)-1 additional
#'   function evaluations for odd nAGQ numbers.
#'
#' @param agqLow The lower bound for adaptive quadrature
#'   log-likelihood. By default this is -Inf; in the original nlmixr's
#'   gnlmm it was -700.
#'
#' @param agqHi The upper bound for adaptive quadrature
#'   log-likelihood.  By default this is Inf; in the original nlmixr's
#'   gnlmm was 400.
#'
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
#' @references
#'
#' Gill, P.E., Murray, W., Saunders, M.A., & Wright,
#' M.H. (1983). Computing Forward-Difference Intervals for Numerical
#' Optimization. Siam Journal on Scientific and Statistical Computing,
#' 4, 310-321.
#'
#' Shi, H.M., Xie, Y., Xuan, M.Q., & Nocedal, J. (2021). Adaptive
#' Finite-Difference Interval Estimation for Noisy Derivative-Free
#' Optimization.
#'
#' @family Estimation control
#' @export
foceiControl <- function(sigdig = 3, #
                         ...,
                         epsilon = NULL, # 1e-4,
                         maxInnerIterations = 1000, #
                         maxOuterIterations = 5000, #
                         n1qn1nsim = NULL, #
                         print = 1L, #
                         printNcol = floor((getOption("width") - 23) / 12), #
                         scaleTo = 1.0, #
                         scaleObjective = 0, #
                         normType = c("rescale2", "mean", "rescale", "std", "len", "constant"), #
                         scaleType = c("nlmixr2", "norm", "mult", "multAdd"), #
                         scaleCmax = 1e5, #
                         scaleCmin = 1e-5, #
                         scaleC = NULL, #
                         scaleC0 = 1e5, #
                         derivEps = rep(20 * sqrt(.Machine$double.eps), 2), #
                         derivMethod = c("switch", "forward", "central"), #
                         derivSwitchTol = NULL, #
                         covDerivMethod = c("central", "forward"), #
                         covMethod = c("r,s", "r", "s", ""), #
                         # norm of weights = 1/0.225
                         #hessEps = (1/0.225*.Machine$double.eps)^(1 / 4), #
                         hessEps =(.Machine$double.eps)^(1/3),
                         #hessEpsLlik =(1/0.225*.Machine$double.eps)^(1/4),
                         hessEpsLlik =(.Machine$double.eps)^(1/3),
                         optimHessType = c("central", "forward"),
                         optimHessCovType=c("central", "forward"),
                         eventType = c("central", "forward"), #
                         centralDerivEps = rep(20 * sqrt(.Machine$double.eps), 2), #
                         lbfgsLmm = 7L, #
                         lbfgsPgtol = 0, #
                         lbfgsFactr = NULL, #
                         eigen = TRUE, #
                         diagXform = c("sqrt", "log", "identity"), #
                         iovXform = c("sd", "var", "logsd", "logvar"), #
                         sumProd = FALSE, #
                         optExpression = TRUE,#
                         literalFix=TRUE,
                         literalFixRes=TRUE,
                         ci = 0.95, #
                         useColor = crayon::has_color(), #
                         boundTol = NULL, #
                         calcTables = TRUE,#
                         noAbort = TRUE, #
                         interaction = TRUE, #
                         cholSEtol = (.Machine$double.eps)^(1 / 3), #
                         cholAccept = 1e-3, #
                         resetEtaP = 0.15, #
                         resetThetaP = 0.05, #
                         resetThetaFinalP = 0.15, #
                         diagOmegaBoundUpper = 5, # diag(omega) = diag(omega)*diagOmegaBoundUpper; =1 no upper
                         diagOmegaBoundLower = 100, # diag(omega) = diag(omega)/diagOmegaBoundLower; = 1 no lower
                         cholSEOpt = FALSE, #
                         cholSECov = FALSE, #
                         fo = FALSE, #
                         covTryHarder = FALSE, #
                         ## Ranking based on run 025
                         ## L-BFGS-B: 20970.53 (2094.004    429.535)
                         ## bobyqa: 21082.34 (338.677    420.754)
                         ## lbfgsb3* (modified for tolerances):
                         ## nlminb: 20973.468 (755.821    458.343)
                         ## mma: 20974.20 (Time: Opt: 3000.501 Cov: 467.287)
                         ## slsqp: 21023.89 (Time: Opt: 460.099; Cov: 488.921)
                         ## lbfgsbLG: 20974.74 (Time: Opt: 946.463; Cov:397.537)
                         outerOpt = c("nlminb",
                                      "bobyqa",
                                      "lbfgsb3c",
                                      "L-BFGS-B",
                                      "mma",
                                      "lbfgsbLG",
                                      "slsqp",
                                      "Rvmmin"), #
                         innerOpt = c("n1qn1", "BFGS"), #
                         ##
                         rhobeg = .2, #
                         rhoend = NULL, #
                         npt = NULL, #
                         ## nlminb
                         rel.tol = NULL, #
                         x.tol = NULL, #
                         eval.max = 4000, #
                         iter.max = 2000, #
                         abstol = NULL, #
                         reltol = NULL, #
                         resetHessianAndEta = FALSE, #
                         stateTrim = Inf, #
                         shi21maxOuter = 0L,
                         shi21maxInner = 20L,
                         shi21maxInnerCov =20L,
                         shi21maxFD=20L,
                         gillK = 10L, #
                         gillStep = 4, #
                         gillFtol = 0, #
                         gillRtol = sqrt(.Machine$double.eps), #
                         gillKcov = 10L, #
                         #gillKcovLlik = 20L,
                         gillKcovLlik = 10L,
                         gillStepCovLlik = 4.5,
                         #gillStepCovLlik = 2,
                         gillStepCov = 2, #
                         gillFtolCov = 0, #
                         gillFtolCovLlik = 0, #
                         rmatNorm = TRUE, #
                         #rmatNormLlik= FALSE, #
                         rmatNormLlik= TRUE, #
                         smatNorm = TRUE, #
                         ## smatNormLlik = FALSE,
                         smatNormLlik = TRUE,
                         covGillF = TRUE, #
                         optGillF = TRUE, #
                         covSmall = 1e-5, #
                         adjLik = TRUE, ## Adjust likelihood by 2pi for FOCEi methods
                         gradTrim = Inf, #
                         maxOdeRecalc = 5, #
                         odeRecalcFactor = 10^(0.5), #
                         gradCalcCentralSmall = 1e-4, #
                         gradCalcCentralLarge = 1e4, #
                         etaNudge = qnorm(1-0.05/2)/sqrt(3), #
                         etaNudge2=qnorm(1-0.05/2) * sqrt(3/5), #
                         nRetries = 3, #
                         seed = 42, #
                         resetThetaCheckPer = 0.1, #
                         etaMat = NULL, #
                         repeatGillMax = 1,#
                         stickyRecalcN = 4, #
                         gradProgressOfvTime = 10, #
                         addProp = c("combined2", "combined1"),
                         badSolveObjfAdj=100, #
                         compress=FALSE, #
                         rxControl=NULL,
                         sigdigTable=NULL,
                         fallbackFD=FALSE,
                         smatPer=0.6,
                         sdLowerFact=0.001,
                         zeroGradFirstReset=TRUE,
                         zeroGradRunReset=TRUE,
                         zeroGradBobyqa=TRUE,
                         mceta=-1L,
                         nAGQ=0,
                         agqLow=-Inf,
                         agqHi=Inf) { #
  if (!is.null(sigdig)) {
    checkmate::assertNumeric(sigdig, lower=1, finite=TRUE, any.missing=TRUE, len=1)
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
    if (is.null(rel.tol)) {
      rel.tol <- 10^(-sigdig - 1)
    }
    if (is.null(x.tol)) {
      x.tol <- 10^(-sigdig - 1)
    }
    if (is.null(derivSwitchTol)) {
      derivSwitchTol <- 2 * 10^(-sigdig - 1)
    }
  }
  if (is.null(sigdigTable)) {
    if (is.null(sigdig)) {
      sigdigTable <- 3L
    } else {
      sigdigTable <- sigdig
    }
  } else {
    checkmate::assertNumeric(sigdigTable, lower=1, finite=TRUE, any.missing=TRUE, len=1)
  }

  checkmate::assertNumeric(epsilon, lower=0, finite=TRUE, any.missing=FALSE, len=1)
  checkmate::assertIntegerish(maxInnerIterations, lower=0, any.missing=FALSE, len=1)
  checkmate::assertIntegerish(maxInnerIterations, lower=0, any.missing=FALSE, len=1)
  checkmate::assertIntegerish(maxOuterIterations, lower=0, any.missing=FALSE, len=1)

  checkmate::assertNumeric(sdLowerFact, lower=0, finite=TRUE, upper=0.1, any.missing=FALSE, len=1)

  if (is.null(n1qn1nsim)) {
    n1qn1nsim <- 10 * maxInnerIterations + 1
  }
  checkmate::assertIntegerish(n1qn1nsim, len=1, lower=1, any.missing=FALSE)
  checkmate::assertIntegerish(print, len=1, lower=0, any.missing=FALSE)
  checkmate::assertIntegerish(printNcol, len=1, lower=1, any.missing=FALSE)
  checkmate::assertNumeric(scaleTo, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(scaleObjective, len=1, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(scaleCmax, lower=0, any.missing=FALSE, len=1)
  checkmate::assertNumeric(scaleCmin, lower=0, any.missing=FALSE, len=1)
  if (!is.null(scaleC)) {
    checkmate::assertNumeric(scaleC, lower=0, any.missing=FALSE)
  }
  checkmate::assertNumeric(scaleC0, lower=0, any.missing=FALSE, len=1)
  checkmate::assertNumeric(derivEps, lower=0, len=2, any.missing=FALSE)
  checkmate::assertNumeric(derivSwitchTol, lower=0, len=1, any.missing=FALSE)
  if (checkmate::testIntegerish(covTryHarder, lower=0, upper=1, any.missing=FALSE, len=1)) {
    covTryHarder <- as.integer(covTryHarder)
  } else {
    checkmate::assertLogical(covTryHarder, any.missing=FALSE, len=1)
    covTryHarder <- as.integer(covTryHarder)
  }

  checkmate::assertNumeric(rhobeg, lower=0, len=1, finite=TRUE, any.missing=FALSE)
  checkmate::assertNumeric(rhoend, lower=0, len=1, finite=TRUE, any.missing=FALSE)
  if (rhoend >= rhobeg) {
    stop("the trust region method needs '0 < rhoend < rhobeg'",
         call.=FALSE)
  }
  if (!is.null(npt)) {
    checkmate::assertIntegerish(npt, lower=1, len=1, any.missing=FALSE)
  }
  checkmate::assertIntegerish(eval.max, lower=1, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(iter.max, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(rel.tol, lower=0, len=1, any.missing=FALSE, finite=TRUE)
  checkmate::assertNumeric(x.tol, lower=0, len=1, any.missing=FALSE, finite=TRUE)
  checkmate::assertNumeric(abstol, lower=0, len=1, any.missing=FALSE, finite=TRUE)
  checkmate::assertNumeric(reltol, lower=0, len=1, any.missing=FALSE, finite=TRUE)

  checkmate::assertIntegerish(gillK, lower=0, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(gillKcov, lower=0, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(gillKcovLlik, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(gillStep, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(gillStepCov, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(gillStepCovLlik, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(gillFtol, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(gillFtolCov, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(gillFtolCovLlik, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(gillRtol, lower=0, len=1, any.missing=FALSE, finite=TRUE)
  # gillRtolCov is calculated in the `inner.cpp`
  if (!checkmate::testIntegerish(rmatNorm, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(rmatNorm, any.missing=FALSE, len=1)
  }
  rmatNorm <- as.integer(rmatNorm)
  if (!checkmate::testIntegerish(rmatNormLlik, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(rmatNormLlik, any.missing=FALSE, len=1)
  }
  rmatNormLlik <- as.integer(rmatNormLlik)
  if (!checkmate::testIntegerish(smatNorm, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(smatNorm, any.missing=FALSE, len=1)
  }
  smatNorm <- as.integer(smatNorm)
  if (!checkmate::testIntegerish(smatNormLlik, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(smatNormLlik, any.missing=FALSE, len=1)
  }
  smatNormLlik <- as.integer(smatNormLlik)
  if (!checkmate::testIntegerish(covGillF, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(covGillF, any.missing=FALSE, len=1)
  }
  covGillF <- as.integer(covGillF)
  if (!checkmate::testIntegerish(optGillF, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(optGillF, any.missing=FALSE, len=1)
  }
  optGillF <- as.integer(optGillF)

  checkmate::assertNumeric(hessEps, lower=0, any.missing=FALSE, len=1)
  checkmate::assertNumeric(hessEpsLlik, lower=0, any.missing=FALSE, len=1)
  checkmate::assertNumeric(centralDerivEps, lower=0, any.missing=FALSE, len=2)

  checkmate::assertIntegerish(lbfgsLmm, lower=1L, any.missing=FALSE, len=1)
  lbfgsLmm <- as.integer(lbfgsLmm)
  checkmate::assertNumeric(lbfgsPgtol, lower=0, any.missing=FALSE, len=1)
  checkmate::assertNumeric(lbfgsFactr, lower=0, any.missing=FALSE, len=1)
  if (!checkmate::testIntegerish(eigen, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(eigen, any.missing=FALSE, len=1)
  }
  eigen <- as.integer(eigen)

  checkmate::assertLogical(sumProd, any.missing=FALSE, len=1)
  checkmate::assertLogical(optExpression, any.missing=FALSE, len=1)
  checkmate::assertLogical(literalFix, any.missing=FALSE, len=1)
  checkmate::assertLogical(literalFixRes, any.missing=FALSE, len=1)

  checkmate::assertNumeric(ci, any.missing=FALSE, len=1, lower=0, upper=1)
  checkmate::assertLogical(useColor, any.missing=FALSE, len=1)
  checkmate::assertNumeric(boundTol, lower=0, any.missing=FALSE, len=1)

  checkmate::assertLogical(calcTables, len=1, any.missing=FALSE)
  if(!checkmate::testIntegerish(noAbort, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(noAbort, len=1, any.missing=FALSE)
  }
  noAbort <- as.integer(noAbort)

  if (!checkmate::testIntegerish(interaction, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(interaction, len=1, any.missing=FALSE)
  }
  interaction <- as.integer(interaction)

  checkmate::assertNumeric(cholSEtol, lower=0, any.missing=FALSE, len=1)
  checkmate::assertNumeric(cholAccept, lower=0, any.missing=FALSE, len=1)

  ## .methodIdx <- c("lsoda"=1L, "dop853"=0L, "liblsoda"=2L);
  ## method <- as.integer(.methodIdx[method]);
  if (checkmate::testIntegerish(scaleType, len=1, lower=1, upper=4, any.missing=FALSE)) {
    scaleType <- as.integer(scaleType)
  } else {
    .scaleTypeIdx <- c("norm" = 1L, "nlmixr2" = 2L, "mult" = 3L, "multAdd" = 4L)
    scaleType <- setNames(.scaleTypeIdx[match.arg(scaleType)], NULL)
  }

  if (checkmate::testIntegerish(optimHessType, len=1, lower=1, upper=3, any.missing=FALSE)) {
    optimHessType <- as.integer(optimHessType)
  } else {
    .optimHessTypeIdx <- c("central" = 1L, "forward" = 3L)
    optimHessType <- setNames(.optimHessTypeIdx[match.arg(optimHessType)], NULL)
  }

  if (checkmate::testIntegerish(optimHessCovType, len=1, lower=1, upper=3, any.missing=FALSE)) {
    optimHessCovType <- as.integer(optimHessCovType)
  } else {
    .optimHessCovTypeIdx <- c("central" = 1L, "forward" = 3L)
    optimHessCovType <- setNames(.optimHessCovTypeIdx[match.arg(optimHessCovType)], NULL)
  }
  if (checkmate::testIntegerish(eventType, len=1, lower=1, upper=3, any.missing=FALSE)) {
    eventType <- as.integer(eventType)
  } else {
    .eventTypeIdx <- c("central" = 2L, "forward" = 3L)
    eventType <- setNames(.eventTypeIdx[match.arg(eventType)], NULL)
  }

  .normTypeIdx <- c("rescale2" = 1L, "rescale" = 2L, "mean" = 3L, "std" = 4L, "len" = 5L, "constant" = 6L)
  if (checkmate::testIntegerish(normType, len=1, lower=1, upper=6, any.missing=FALSE)) {
    normType <- as.integer(normType)
  } else {
    normType <- setNames(.normTypeIdx[match.arg(normType)], NULL)
  }
  .methodIdx <- c("forward" = 0L, "central" = 1L, "switch" = 3L)
  if (checkmate::testIntegerish(derivMethod, len=1, lower=0L, upper=3L, any.missing=FALSE)) {
    derivMethod <- as.integer(derivMethod)
  } else {
    derivMethod <- match.arg(derivMethod)
    derivMethod <- setNames(.methodIdx[derivMethod], NULL)
  }
  if (checkmate::testIntegerish(covDerivMethod, len=1, lower=0L, upper=3L, any.missing=FALSE)) {
    covDerivMethod <- as.integer(covDerivMethod)
  } else {
    covDerivMethod <- match.arg(covDerivMethod)
    covDerivMethod <- setNames(.methodIdx[covDerivMethod], NULL)
  }
  if (checkmate::testIntegerish(covMethod, len=1, lower=0L, upper=3L, any.missing=FALSE)) {
    covMethod <- as.integer(covMethod)
  } else if (rxode2::rxIs(covMethod, "character")) {
    if (all(covMethod == "")) {
      covMethod <- 0L
    } else {
      covMethod <- match.arg(covMethod)
      .covMethodIdx <- c("r,s" = 1L, "r" = 2L, "s" = 3L)
      covMethod <- setNames(.covMethodIdx[match.arg(covMethod)], NULL)
    }
  }
  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% .foceiControlInternal)]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste
    (paste0("'", .bad, "'", sep=""), collapse=", "),
    call.=FALSE)
  }
  .skipCov <- NULL
  if (!is.null(.xtra$skipCov)) {
    .skipCov <- .xtra$skipCov
  }
  .outerOptTxt <- "custom"
  if (!is.null(.xtra$outerOptTxt)) {
    .outerOptTxt <- .xtra$outerOptTxt
  }
  outerOptFun <- NULL
  if (!is.null(.xtra$outerOptFun)) {
    outerOptFun <- .xtra$outerOptFun
  } else if (rxode2::rxIs(outerOpt, "character")) {
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
    } else {
      if (checkmate::testIntegerish(outerOpt, lower=0, upper=1, len=1)) {
        outerOpt <- as.integer(outerOpt)
      } else {
        .outerOptIdx <- c("L-BFGS-B" = 0L, "lbfgsb3c" = 1L)
        outerOpt <- .outerOptIdx[outerOpt]
        if (outerOpt == 1L) {
          rxode2::rxReq("lbfgsb3c")
        }
      }
      outerOptFun <- NULL
    }
  } else if (is(outerOpt, "function")) {
    outerOptFun <- outerOpt
    outerOpt <- -1L
  }
  if (checkmate::testIntegerish(innerOpt, lower=1, upper=2, len=1)) {
    innerOpt <- as.integer(innerOpt)
  } else {
    .innerOptFun <- c("n1qn1" = 1L, "BFGS" = 2L)
    innerOpt <- setNames(.innerOptFun[match.arg(innerOpt)], NULL)
  }
  if (!is.null(.xtra$resetEtaSize)) {
    .resetEtaSize <- .xtra$resetEtaSize
  } else {
    checkmate::assertNumeric(resetEtaP, lower=0, upper=1, len=1)
    if (resetEtaP > 0 & resetEtaP < 1) {
      .resetEtaSize <- qnorm(1 - (resetEtaP / 2))
    } else if (resetEtaP <= 0) {
      .resetEtaSize <- Inf
    } else {
      .resetEtaSize <- 0
    }
  }
  if (!is.null(.xtra$resetThetaSize)) {
    .resetThetaSize <- .xtra$resetThetaSize
  } else {
    checkmate::assertNumeric(resetThetaP, lower=0, upper=1, len=1)
    if (resetThetaP > 0 & resetThetaP < 1) {
      .resetThetaSize <- qnorm(1 - (resetThetaP / 2))
    } else if (resetThetaP <= 0) {
      .resetThetaSize <- Inf
    } else {
      stop("cannot always reset THETAs", call.=FALSE)
    }
  }
  if (!is.null(.xtra$resetThetaFinalSize)) {
    .resetThetaFinalSize <- .xtra$resetThetaFinalSize
  } else {
    checkmate::assertNumeric(resetThetaFinalP, lower=0, upper=1, len=1)
    if (resetThetaFinalP > 0 & resetThetaFinalP < 1) {
      .resetThetaFinalSize <- qnorm(1 - (resetThetaFinalP / 2))
    } else if (resetThetaP <= 0) {
      .resetThetaFinalSize <- Inf
    } else {
      stop("cannot always reset THETAs", call.=FALSE)
    }
  }
  if (checkmate::testIntegerish(addProp, lower=1, upper=1, len=1)) {
    addProp <- c("combined1", "combined2")[addProp]
  } else {
    addProp <- match.arg(addProp)
  }
  checkmate::assertLogical(compress, any.missing=FALSE, len=1)
  if (!is.null(.xtra$genRxControl)) {
    genRxControl <- .xtra$genRxControl
  } else {
    genRxControl <- FALSE
    if (is.null(rxControl)) {
      rxControl <- rxode2::rxControl(sigdig=sigdig,
                                     maxsteps=500000L)
      genRxControl <- TRUE
    } else if (is.list(rxControl)) {
      rxControl <- do.call(rxode2::rxControl, rxControl)
    }
    if (!inherits(rxControl, "rxControl")) {
      stop("rxControl needs to be ode solving options from rxode2::rxControl()",
           call.=FALSE)
    }
  }
  checkmate::assertNumeric(diagOmegaBoundUpper, lower=1, len=1, any.missing=FALSE, finite=TRUE)
  checkmate::assertNumeric(diagOmegaBoundLower, lower=1, len=1, any.missing=FALSE, finite=TRUE)

  if (!checkmate::testIntegerish(cholSEOpt, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(cholSEOpt, any.missing=FALSE, len=1)
  }
  cholSEOpt <- as.integer(cholSEOpt)

  if (!checkmate::testIntegerish(cholSECov, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(cholSECov, any.missing=FALSE, len=1)
  }
  cholSECov <- as.integer(cholSECov)

  if (!checkmate::testIntegerish(fo, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(fo, any.missing=FALSE, len=1)
  }
  fo <- as.integer(fo)

  if (!checkmate::testIntegerish(resetHessianAndEta, lower=0, upper=1, any.missing=FALSE, len=1)) {
    checkmate::assertLogical(resetHessianAndEta, any.missing=FALSE, len=1)
  }
  resetHessianAndEta <- as.integer(resetHessianAndEta)
  checkmate::assertNumeric(stateTrim, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(covSmall, lower=0, any.missing=FALSE, finite=TRUE)
  checkmate::assertLogical(adjLik, any.missing=FALSE, len=1)
  checkmate::assertNumeric(gradTrim, any.missing=FALSE, len=1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing=FALSE, len=1)
  checkmate::assertNumeric(odeRecalcFactor, len=1, lower=1, any.missing=FALSE)
  checkmate::assertNumeric(gradCalcCentralSmall, len=1, lower=0, any.missing=FALSE, finite=TRUE)
  checkmate::assertNumeric(gradCalcCentralLarge, len=1, lower=0, any.missing=FALSE, finite=TRUE)
  checkmate::assertNumeric(etaNudge, len=1, lower=0, any.missing=FALSE, finite=TRUE)
  checkmate::assertNumeric(etaNudge2, len=1, lower=0, any.missing=FALSE, finite=TRUE)
  checkmate::assertIntegerish(nRetries, lower=0, any.missing=FALSE)
  if (!is.null(seed)) {
    checkmate::assertIntegerish(seed, any.missing=FALSE, min.len=1)
  }
  checkmate::assertNumeric(resetThetaCheckPer, lower=0, upper=1, any.missing=FALSE, finite=TRUE)
  checkmate::assertIntegerish(repeatGillMax, any.missing=FALSE, lower=0, len=1)
  checkmate::assertIntegerish(stickyRecalcN, any.missing=FALSE, lower=0, len=1)
  checkmate::assertNumeric(gradProgressOfvTime, any.missing=FALSE, lower=0, len=1)
  checkmate::assertNumeric(badSolveObjfAdj, any.missing=FALSE, len=1)
  checkmate::assertLogical(fallbackFD, any.missing=FALSE, len=1)
  checkmate::assertLogical(zeroGradFirstReset, any.missing=TRUE, len=1)
  checkmate::assertLogical(zeroGradRunReset, any.missing=FALSE, len=1)
  checkmate::assertLogical(zeroGradBobyqa, any.missing=TRUE, len=1)

  checkmate::assertIntegerish(shi21maxOuter, lower=0, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(shi21maxInner, lower=0, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(shi21maxInnerCov, lower=0, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(shi21maxFD, lower=0, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(mceta, lower=-1, len=1,any.missing=FALSE)

  checkmate::assertNumeric(smatPer, any.missing=FALSE, lower=0, upper=1, len=1)
  checkmate::assertIntegerish(nAGQ, lower=0, len=1, any.missing=FALSE)
  checkmate::assertNumeric(agqHi, len=1, any.missing=FALSE)
  checkmate::assertNumeric(agqLow, len=1, any.missing=FALSE)
  .ret <- list(
    maxOuterIterations = as.integer(maxOuterIterations),
    maxInnerIterations = as.integer(maxInnerIterations),
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
    eigen = eigen,
    diagXform = match.arg(diagXform),
    iovXform = match.arg(iovXform),
    sumProd = sumProd,
    optExpression = optExpression,
    literalFix=literalFix,
    literalFixRes=literalFixRes,
    outerOpt = as.integer(outerOpt),
    ci = as.double(ci),
    sigdig = as.double(sigdig),
    sigdigTable=sigdigTable,
    scaleObjective = as.double(scaleObjective),
    useColor = useColor,
    boundTol = as.double(boundTol),
    calcTables = calcTables,
    printNcol = as.integer(printNcol),
    noAbort = noAbort,
    interaction = interaction,
    cholSEtol = as.double(cholSEtol),
    hessEps = as.double(hessEps),
    hessEpsLlik = as.double(hessEpsLlik),
    optimHessType=optimHessType,
    optimHessCovType=optimHessCovType,
    cholAccept = as.double(cholAccept),
    resetEtaSize = as.double(.resetEtaSize),
    resetThetaSize = as.double(.resetThetaSize),
    resetThetaFinalSize = as.double(.resetThetaFinalSize),
    diagOmegaBoundUpper = diagOmegaBoundUpper,
    diagOmegaBoundLower = diagOmegaBoundLower,
    cholSEOpt = cholSEOpt,
    cholSECov = cholSECov,
    fo = fo,
    covTryHarder = covTryHarder,
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
    resetHessianAndEta = resetHessianAndEta,
    stateTrim = as.double(stateTrim),
    gillK = as.integer(gillK),
    gillKcov = as.integer(gillKcov),
    gillKcovLlik = as.integer(gillKcovLlik),
    gillRtol = as.double(gillRtol),
    gillStep = as.double(gillStep),
    gillStepCov = as.double(gillStepCov),
    gillStepCovLlik = as.double(gillStepCovLlik),
    scaleType = scaleType,
    normType = normType,
    scaleC = scaleC,
    scaleCmin = as.double(scaleCmin),
    scaleCmax = as.double(scaleCmax),
    scaleC0 = as.double(scaleC0),
    outerOptTxt = .outerOptTxt,
    rmatNorm = rmatNorm,
    rmatNormLlik = rmatNormLlik,
    smatNorm = smatNorm,
    smatNormLlik = smatNormLlik,
    covGillF = covGillF,
    optGillF = optGillF,
    gillFtol = as.double(gillFtol),
    gillFtolCov = as.double(gillFtolCov),
    gillFtolCovLlik = as.double(gillFtolCovLlik),
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
    eventType = eventType,
    gradProgressOfvTime = gradProgressOfvTime,
    addProp = addProp,
    badSolveObjfAdj=badSolveObjfAdj,
    compress=compress,
    rxControl=rxControl,
    genRxControl=genRxControl,
    skipCov=.skipCov,
    fallbackFD=fallbackFD,
    shi21maxOuter=shi21maxOuter,
    shi21maxInner=shi21maxInner,
    shi21maxInnerCov=shi21maxInnerCov,
    shi21maxFD=shi21maxFD,
    smatPer=smatPer,
    sdLowerFact=sdLowerFact,
    zeroGradFirstReset=zeroGradFirstReset,
    zeroGradRunReset=zeroGradRunReset,
    zeroGradBobyqa=zeroGradBobyqa,
    mceta=as.integer(mceta),
    nAGQ=as.integer(nAGQ),
    agqHi=as.double(agqHi),
    agqLow=as.double(agqLow)
  )
  if (length(etaMat) == 1L && is.na(etaMat)) {
    .ret$etaMat <- NA
  } else if (!is.null(etaMat)) {
    .doWarn <- TRUE
    if (inherits(etaMat, "nlmixr2FitCore")) {
      etaMat <- etaMat$etaMat
      .doWarn <- FALSE
    }
    if (.doWarn && missing(maxInnerIterations)) {
      warning(sprintf("using 'etaMat' assuming 'maxInnerIterations=%d', set 'maxInnerIterations' explicitly to avoid this warning", maxInnerIterations))
    }
    checkmate::assertMatrix(etaMat, mode="double", any.missing=FALSE, min.rows=1, min.cols=1)
    .ret$etaMat <- etaMat
  }
  class(.ret) <- "foceiControl"
  .ret
}

#' @export
rxUiDeparse.foceiControl <- function(object, var) {
  .ret <- foceiControl()
  .outerOpt <- character(0)
  if (object$outerOpt == -1L && object$outerOptTxt == "custom") {
    warning("functions for `outerOpt` cannot be deparsed, reset to default",
            call.=FALSE)
  } else if (!(object$outerOptTxt %in% c("nlminb", "stats::optimize"))) {
    .outerOpt <- paste0("outerOpt=", deparse1(object$outerOptTxt))
  }
  .w <- .deparseDifferent(.ret, object, .foceiControlInternal)
  if (length(.w) == 0 && length(.outerOpt) == 0) {
    return(str2lang(paste0(var, " <- foceiControl()")))
  }
  .n <- names(.ret)[.w]
  .n <- .n[.n != "outerOpt"]
  .retD <- c(vapply(.n, function(x) {
    .val <- .deparseShared(x, object[[x]])
    if (!is.na(.val)) {
      return(.val)
    }
    if (x == "innerOpt") {
      .innerOptFun <- c("n1qn1" = 1L, "BFGS" = 2L)
      paste0("innerOpt =", deparse1(names(.innerOptFun[which(object[[x]] == .innerOptFun)])))
    } else if (x %in% c("optimHessType", "optimHessCovType")) {
      .methodIdx <- c("central" = 1L, "forward" = 3L)
      paste0(x, " =", deparse1(names(.methodIdx[which(object[[x]] == .methodIdx)])))
    } else if (x == "eventType") {
      .methodIdx <- c("central" = 2L, "forward" = 3L)
      paste0(x, " =", deparse1(names(.methodIdx[which(object[[x]] == .methodIdx)])))
    } else if (x %in% c("derivMethod", "covDerivMethod")) {
      .methodIdx <- c("forward" = 0L, "central" = 2L, "switch" = 3L)
      paste0(x, " =", deparse1(names(.methodIdx[which(object[[x]] == .methodIdx)])))
    } else if (x == "covMethod") {
      if (object[[x]] == 0L) {
        paste0(x, " = \"\"")
      } else {
        .covMethodIdx <- c("r,s" = 1L, "r" = 2L, "s" = 3L)
        paste0(x, " =", deparse1(names(.covMethodIdx[which(object[[x]] == .covMethodIdx)])))
      }
    } else {
      paste0(x, "=", deparse1(object[[x]]))
    }
  }, character(1)), .outerOpt)
  str2lang(paste(var, " <- foceiControl(", paste(.retD, collapse=","),")"))
}
