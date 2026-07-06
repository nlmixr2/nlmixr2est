.foceiControlInternal <- c("genRxControl", "resetEtaSize", "foceType",
                           "resetThetaSize", "resetThetaFinalSize",
                           "outerOptFun", "outerOptTxt", "skipCov",
                           "foceiMuRef", "foceiMuCovEta", "predNeq", "nfixed", "nomega",
                           "neta", "ntheta", "nF", "printTop", "needOptimHess",
                           "iterPrintControl", "est", "foceiMuModel", "foceiMuGroupTheta",
                           "foceiMuGroupEta", "foceiMuGroupCovStart", "foceiMuGroupCovCount",
                           "foceiMuGroupCovTheta", "foceiMuGroupCovUserFixed",
                           "foceiMuGroupCovBounded",
                           "foceiMuGroupCovData", "foceiMuGroupTol", "foceiMuGroupMaxCycles")

#' Control Options for FOCEi
#'
#' @param sigdig Optimization significant digits; controls the inner/outer
#'   optimization tolerance (\code{10^-sigdig}), ODE solver tolerance
#'   (\code{0.5*10^(-sigdig-2)}, or \code{0.5*10^(-sigdig-1.5)} for
#'   sensitivity/steady-state with liblsoda), and boundary check tolerance
#'   (\code{5*10^(-sigdig+1)}).
#'
#' @param sigdigTable Significant digits in the final output table.
#'   If not specified, then it matches the significant digits in the
#'   `sigdig` optimization algorithm.  If `sigdig` is NULL, use 3.
#'
#' @param epsilon Precision of estimate for n1qn1 optimization.
#'
#' @inheritParams iterPrintParams
#'
#' @param scaleTo Scale the initial parameter estimate to this value.
#'     By default this is 1.  When zero or below, no scaling is performed.
#'
#' @param scaleObjective Scale the initial objective function to this
#'     value.  By default this is 0 (meaning do not scale)
#'
#' @param derivEps Forward difference tolerances (relative, absolute); step
#'     size \code{h = abs(x)*derivEps[1] + derivEps[2]}.
#'
#' @param derivMethod Derivative method for the outer problem: "switch",
#'     "central", or "forward". "switch" starts forward and toggles to
#'     central when \code{abs(delta(OFV)) <= derivSwitchTol}.
#'
#' @param derivSwitchTol The tolerance to switch forward to central
#'     differences.
#'
#' @param covDerivMethod indicates the method for calculating the
#'     derivatives while calculating the covariance components
#'     (Hessian and S).
#'
#' @param covMethod Method for calculating covariance, where R is the
#'     Hessian and S the sum of individual gradient cross-products (at the
#'     empirical Bayes estimates): \code{"r,s"} sandwich
#'     (\code{solve(R)\%*\%S\%*\%solve(R)}), \code{"r"} Hessian-based
#'     (\code{solve(R)}), \code{"s"} cross-product-based (\code{solve(S)}), or
#'     \code{""} to skip the covariance step.
#'
#' @param covType covariance R-matrix (Hessian) source, \code{"fd"} (default) or
#'     \code{"analytic"}.  \code{"analytic"} uses the exact analytic
#'     observed-information R-matrix and can additionally return the residual and
#'     Omega standard errors.  It applies to FOCEI and FOCE fits with additive,
#'     proportional, or combined additive-plus-proportional error, and covers
#'     mu-referenced, covariate, and other non-mu-referenced structural parameters
#'     (and non-mu-referenced etas) as well as SD-scale inter-occasion variability.
#'     It defaults to \code{covMethod = "r"} (the observed-information \eqn{R^{-1}},
#'     computed even when \code{covMethod = ""}); an explicit \code{covMethod = "r,s"}
#'     or \code{"s"} is honored, with the analytic R feeding the native
#'     finite-difference sandwich / S-matrix.  It emits a message and falls back to
#'     the finite-difference Hessian for anything out of scope (FO, \code{nAGQ > 1},
#'     censoring, DV-transformed error, bounded-parameter transforms, a structural
#'     theta shared by two etas, or non-SD \code{iovXform}) -- and, for a
#'     pure-proportional variance that vanishes at a near-zero model prediction,
#'     where the observed information is ill-conditioned.
#'
#' @param covSolveTol absolute/relative ODE tolerance for the augmented-sensitivity
#'     solves behind \code{covType="analytic"}.  \code{NULL} (default) derives a
#'     tight tolerance from \code{sigdig}; supply a number to override it.
#'
#' @param covFull controls the shape of \code{fit$cov}.  \code{FALSE} (default)
#'     installs only the structural-theta block (the NONMEM-matched theta
#'     covariance, matching the historical finite-difference \code{fit$cov} shape
#'     for backwards compatibility); \code{TRUE} installs the full theta + residual
#'     sigma + Omega covariance -- assembled analytically for
#'     \code{covType="analytic"}, or by central finite differences of the objective
#'     over the same parameter set for \code{covType="fd"} (perturbing Omega on the
#'     variance-covariance scale, with the per-parameter Gill (1983) step and the
#'     5-point/4-point stencils that \code{foceiCalcR} uses).  The theta standard
#'     errors are identical either way.
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
#' @param optimHessType Hessian type for numeric-difference individual
#'   Hessians in generalized log-likelihood estimation: "central" (matches
#'   R's `optimHess()`, default) or "forward" (faster).
#'
#' @param optimHessCovType Hessian type for numeric-difference individual
#'   Hessians used for the covariance step/final likelihood: "central"
#'   (more accurate, used here) or "forward".
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
#' @param centralDerivEps Central difference tolerances (relative,
#'   absolute); step size \code{h = abs(x)*derivEps[1] + derivEps[2]}.
#'
#' @param lbfgsLmm An integer giving the number of BFGS updates
#'     retained in the "L-BFGS-B" method, It defaults to 7.
#'
#' @param lbfgsPgtol Projected-gradient convergence tolerance for
#'     "L-BFGS-B": iteration stops when
#'     \code{max(| proj g_i |) <= lbfgsPgtol}. Defaults to `0` (check
#'     suppressed).
#'
#' @param lbfgsFactr Convergence factor for "L-BFGS-B": converges when the
#'     objective reduction is within \code{lbfgsFactr * .Machine$double.eps}.
#'     Default `1e10` (~4 sigdigs, \code{2e-6}).
#'
#' @param diagXform Transformation used on the diagonal of
#'     \code{chol(solve(omega))} (the FOCEi-estimated parameters): one of
#'     \code{"sqrt"} (default), \code{"log"}, or \code{"identity"}.
#'
#' @param iovXform Transformation used on the diagonal of the IOV: one of
#'     \code{"sd"}, \code{"var"}, \code{"logsd"}, or \code{"logvar"}.
#'
#' @param sumProd Is a boolean indicating if the model should change
#'     multiplication to high precision multiplication and sums to
#'     high precision sums using the PreciseSums package.  By default
#'     this is \code{FALSE}.
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
#' @param noAbort Boolean to indicate if you should abort the FOCEi
#'     evaluation if it runs into troubles.  (default TRUE)
#'
#' @param interaction Boolean indicate FOCEi should be used (TRUE)
#'     instead of FOCE (FALSE)
#'
#' @param foce Controls how FOCE (\code{interaction = FALSE}) evaluates the
#'     residual variance R in the inner objective; ignored for FOCEi.  Either
#'     \code{"nonmem"} (default) or \code{"foce+"}:
#'
#'     \itemize{
#'
#'     \item \code{"nonmem"} freezes R at the \code{eta = 0} population
#'     prediction and holds it constant across the inner optimization, matching
#'     NONMEM's FOCE.  Advantage: reproduces NONMEM FOCE objective and standard
#'     errors, and an ODE model agrees with its closed-form (\code{linCmt})
#'     equivalent.  Disadvantage: R ignores the individual (conditional)
#'     heteroscedasticity, so it can be slightly less accurate than
#'     \code{"foce+"} for proportional/combined error.
#'
#'     \item \code{"foce+"} evaluates R at the current conditional
#'     \code{eta} (the live variance), keeping the truncated FOCE
#'     inner gradient.  Advantage: uses the conditional variance and
#'     is a bit more accurate than NONMEM's FOCE in some cases.
#'     Disadvantage: does not match NONMEM FOCE and is unsupported by
#'     \code{covType = "analytic"} (falls back to the
#'     finite-difference covariance).  This was the FOCE behavior in
#'     \pkg{nlmixr2est} 6.0.1 and earlier. This does not use the
#'     gradient of \code{eta} like the full \code{focei} method, so it
#'     is not as accurate as \code{focei}.
#'
#'     }
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
#' @param resetEtaP P-value for resetting an individual ETA to 0 during
#'     optimization, based on a z-test of \code{chol(omega^-1) \%*\% eta}
#'     or \code{eta/sd(allEtas)}. `0` = never reset, `1` = always reset.
#'
#' @param resetThetaP P-value for resetting mu-referenced THETAs based on
#'     ETA drift, checked at the start and near a local minimum (see
#'     \code{resetThetaCheckPer}). `0` = never reset; `1` is not allowed.
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
#' @param muModel Selects the mu-referenced-FOCEI-family regression variant
#'     for theta/eta in a mu-ref covariate relationship (see
#'     \code{muRefCovAlg}): \code{"none"} (default, ordinary FOCEI);
#'     \code{"lin"} (\code{mufocei}/\code{mufoce}/\code{muagq}/
#'     \code{mulaplace}: population theta and covariate coefficient(s) per
#'     mu-ref-covariate group are excluded from the outer optimizer and
#'     re-derived in C++ by closed-form OLS regression of each subject's
#'     back-calculated value on the covariate(s), residual becomes that
#'     subject's eta; repeats until convergence, see \code{muModelTol}/
#'     \code{muModelMaxCycles}); or \code{"irls"}
#'     (\code{irlsfocei}/\code{irlsfoce}/\code{irlsagq}/\code{irlslaplace}:
#'     same mechanism, reweighted by inner-optimization curvature).
#'
#'     A mu-ref-covariate theta with a finite bound falls back to ordinary
#'     bounded outer-optimizer handling with a warning (a bound on the
#'     group's population theta excludes the whole group; a bound on one
#'     covariate coefficient excludes only that covariate).
#'
#' @param muRefCovAlg When `TRUE` (default), algebraic expressions that can
#'     be mu-referenced are internally rewritten as mu-referenced
#'     covariates and restored after optimization. Mirrors
#'     \code{saemControl(muRefCovAlg=)}/\code{nlmeControl(muRefCovAlg=)};
#'     for \code{foceiControl()} only takes effect when
#'     \code{muModel != "none"}.
#'
#' @param muModelTol Convergence tolerance for the mu-referenced-FOCEI-family
#'     "re-optimize etas, then regress" cycle (\code{muModel != "none"}):
#'     repeats until the max mu-group theta change drops below this value
#'     or \code{muModelMaxCycles} is reached.
#'
#' @param muModelMaxCycles Maximum number of "re-optimize etas, regress"
#'     cycles per outer iteration (see \code{muModel}, \code{muModelTol}).
#'
#' @param diagOmegaBoundUpper Upper bound of the diagonal omega matrix, as
#'     \code{diag(omega)*diagOmegaBoundUpper}. `1` = no upper bound.
#'
#' @param diagOmegaBoundLower Lower bound of the diagonal omega matrix, as
#'     \code{diag(omega)/diagOmegaBoundLower}. `1` = no lower bound.
#'
#' @param rhobeg Initial trust region radius for the bobyqa outer optimizer
#'     (with `rhoend`, must satisfy `0 < rhoend < rhobeg`). Default `0.2`
#'     (20% of scaled parameters); adjusted upward if smaller than
#'     `abs(upper-lower)/2`. (bobyqa)
#'
#' @param rhoend Final trust region radius. If not defined,
#'     `10^(-sigdig-1)` is used. (bobyqa)
#'
#' @param npt Number of points for bobyqa's quadratic approximation to the
#'     objective; must be in `[n+2, (n+1)(n+2)/2]`. Defaults to `2*n + 1`.
#'     (bobyqa)
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
#' @param gillK Max steps to determine the optimal forward/central
#'     difference step size per parameter (Gill 1983). `0` = no optimal
#'     step size determined.
#'
#' @param gillKcovLlik Same as \code{gillK} but for the generalized focei
#'   log-likelihood method (Gill 1986).
#'
#' @param gillRtol The relative tolerance used for Gill 1983
#'     determination of optimal step size.
#'
#' @param scaleType The scaling scheme for nlmixr2: \code{"nlmixr2"}
#'     (default) scales as \code{(current-init)*scaleC[i] + scaleTo}, with
#'     \code{scaleTo} from \code{normType} and scales from \code{scaleC};
#'     \code{"norm"} uses the simple scaling from \code{normType};
#'     \code{"mult"} scales multiplicatively as \code{current/init*scaleTo};
#'     \code{"multAdd"} scales linearly (\code{(current-init)+scaleTo}) for
#'     parameters in an exponential block (e.g. \code{exp(theta)}) and
#'     multiplicatively otherwise.
#'
#' @param scaleC Scaling constant used with \code{scaleType="nlmixr2"};
#'     when not specified, chosen by parameter type to keep gradient sizes
#'     similar on a log scale: `1` for exp()-transformed/power/boxCox/
#'     yeoJohnson parameters, `0.5*abs(est)` for additive/proportional/
#'     lognormal error parameters, `abs(1/digamma(est+1))` for factorials,
#'     and `log(abs(est))*abs(est)` for log-scale parameters. May be set
#'     explicitly per parameter if these defaults don't apply well.
#'
#' @param scaleC0 Number to adjust the scaling factor by if the initial
#'     gradient is zero.
#'
#' @param scaleCmax Maximum value of the scaleC to prevent overflow.
#'
#' @param scaleCmin Minimum value of the scaleC to prevent underflow.
#'
#' @param normType Parameter normalization/scaling used to get scaled
#'     initial values for \code{scaleType}, of the form
#'     \code{Vscaled = (Vunscaled-C1)/C2} (see
#'     \href{https://en.wikipedia.org/wiki/Feature_scaling}{Feature Scaling};
#'     \code{rescale2} follows the
#'     \href{http://apmonitor.com/me575/uploads/Main/optimization_book.pdf}{OptdesX}
#'     manual): \code{"rescale2"} scales all parameters to (-1, 1);
#'     \code{"rescale"} (min-max) scales to (0, 1); \code{"mean"} centers on
#'     the mean with range (0, 1); \code{"std"} standardizes by mean/sd;
#'     \code{"len"} scales to unit (Euclidean) length; \code{"constant"}
#'     performs no normalization (\code{C1=0}, \code{C2=1}).
#'
#' @param gillStep When looking for the optimal forward difference
#'     step size, this is This is the step size to increase the
#'     initial estimate by.  So each iteration the new step size =
#'     (prior step size)*gillStep
#'
#' @param gillFtol The gillFtol is the gradient error tolerance that
#'     is acceptable before issuing a warning/error about the gradient estimates.
#'
#' @param gillKcov Max steps to determine the optimal forward/central
#'     difference step size per parameter (Gill 1983) during the
#'     covariance step. `0` = no optimal step size determined.
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
#' @param covSmall Small number used to compare covariance estimates
#'     (sandwich vs R/S matrix) before rejecting one as too small to be
#'     the final covariance estimate.
#'
#' @param adjLik When `TRUE`, adjusts the likelihood by the 2*pi constant
#'     nlmixr2's objective function otherwise omits (to match NONMEM),
#'     more closely matching nlme/SAS likelihood approximations. The
#'     objective function itself always matches NONMEM regardless.
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
#' @param etaNudge When n1qn1 optimization of an ETA (starting at zero)
#'   misbehaves, reset the Hessian and nudge the ETA up by this value, then
#'   down if it still doesn't move. Defaults to
#'   `qnorm(1-0.05/2)*1/sqrt(3)`. Falls back to \code{etaNudge2}, then to
#'   zero (stop optimizing) if unsuccessful.
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
#' @param indTolRelax When `TRUE` (default), only subjects whose ODE
#'     solve produced NaN/Inf have their tolerances relaxed, and the
#'     relaxed tolerance persists across optimizer calls (sticky).
#'     When `FALSE`, all subjects have their tolerances relaxed on
#'     each retry and tolerances are reset afterward.
#'
#' @param nRetries If FOCEi doesn't fit with the current parameter
#'     estimates, randomly sample new parameter estimates and restart
#'     the problem.  This is similar to 'PsN' resampling.
#'
#' @param eventType Event gradient type for dosing events; Can be
#'   "central" or "forward"
#'
#' @param eventSens How sensitivities of dosing/event parameters
#'   (absorption lag time, bioavailability, infusion rate and duration,
#'   etc.) are computed.  `"fd"` uses the legacy finite
#'   differences.  `"jump"` (the default) uses the analytic event ("jump")
#'   sensitivities provided by `rxode2`, which add accuracy and can speed
#'   up the gradient/Hessian by avoiding the extra finite-difference
#'   solves for these parameters.
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
#' @param etaMat Initial (or final) ETA estimates; can also be a prior fit,
#'   whose final ETAs are then used as initial values. By default, uses the
#'   last fit's ETAs if supplied, else all ETAs start at zero (`NULL`).
#'   `NA` disables reuse from a prior fit.
#'
#' @param addProp Type of additive-plus-proportional error: `"combined1"`,
#'   where standard deviations add:
#'   \deqn{y = f + (a + b\times f^c) \times \varepsilon}{y = f + (a + b*f^c)*err};
#'   or `"combined2"`, where variances add:
#'   \deqn{y = f + \sqrt{a^2 + b^2\times f^{2\times c}} \times \varepsilon}{y = f + sqrt(a^2 + b^2*(f^c)^2)*err}.
#'   Here y = observed, f = predicted, a = additive sd, b = proportional/power
#'   sd, c = power exponent (1 in the proportional case).
#'
#' @param odeRecalcFactor The ODE recalculation factor when ODE
#'   solving goes bad, this is the factor the rtol/atol is reduced
#'
#' @param rxControl `rxode2` ODE solving options during fitting, created with `rxControl()`
#'
#' @param fallbackFD Fallback to the finite differences if the
#'   sensitivity equations do not solve.
#'
#' @param smatPer Percentage of failed per-individual parameter gradients
#'   (replaced with the overall parameter gradient) out of the total
#'   (`ntheta*nsub`) above which the S matrix is considered bad.
#'
#' @param sdLowerFact Factor multiplying the estimate when the lower bound
#'   is zero for a standard-deviation error parameter (add.sd, prop.sd,
#'   etc); e.g. estimate 0.15 with lower bound 0 assumes a lower bound of
#'   0.00015. `0` disables this.
#'
#' @param zeroGradFirstReset When `TRUE` (default), reset a zero first
#'   gradient to `sqrt(.Machine$double.eps)` instead of erroring; `FALSE`
#'   errors; `NA` ignores it only on the last reset attempt.
#'
#' @param zeroGradRunReset When `TRUE` (default), reset a zero gradient
#'   encountered mid-run to `sqrt(.Machine$double.eps)` instead of erroring.
#'
#' @param zeroGradBobyqa When `TRUE` (default), a zero-gradient reset
#'   switches to the gradient-free bobyqa method; `NA` only does so for the
#'   first zero gradient.
#'
#' @param mceta Monte Carlo sampling for the best initial ETA estimate
#'   (based on `omega`): `-1` (default) uses the last eta; `0` uses eta=0
#'   for each inner optimization; for `n>0`, the last eta, eta=0, and n-1
#'   etas sampled from omega are each evaluated and the best (by inner
#'   objective) is used.
#'
#' @param nAGQ Number of Gauss-Hermite adaptive quadrature points. `0`
#'   disables AGQ; `1` is equivalent to Laplace. Cost grows quickly with
#'   ETAs: once the EBE is found, expect `nAGQ^neta` (even `nAGQ`) or
#'   `(nAGQ^neta)-1` (odd `nAGQ`) additional evaluations per subject.
#'
#' @param agqLow The lower bound for adaptive quadrature
#'   log-likelihood. By default this is -Inf; in the original nlmixr's
#'   gnlmm it was -700.
#'
#' @param agqHi The upper bound for adaptive quadrature
#'   log-likelihood.  By default this is Inf; in the original nlmixr's
#'   gnlmm was 400.
#'
#' @param boundedTransform When `TRUE` (default), bounded parameters are
#'   transformed for unbounded optimization methods and back-transformed
#'   for final estimates. `FALSE` optimizes on the original scale with
#'   bounds passed to the optimizer. `NA` transforms for optimization but
#'   skips the final back-transform.
#'
#' @param eventSens Controls how dosing/event-parameter (`alag`, `F`,
#'   `rate`, `dur`) sensitivities are computed for THETA/ETA gradients:
#'   `"jump"` (default) uses rxode2's analytic event sensitivities; `"fd"`
#'   uses the legacy finite-difference behavior.
#'
#' @param sensMethod Method used to compute the ODE parameter sensitivities:
#'   `"default"` (the default) defers to the global option
#'   `getOption("nlmixr2est.adjoint")` (itself `"forward"` by default);
#'   `"forward"` uses the classic variational (forward) sensitivity ODEs;
#'   `"adjoint"` uses the in-engine discrete adjoint with the matching adjoint
#'   (`s`) method.
#'
#' @inheritParams rxode2::rxSolve
#' @inheritParams minqa::bobyqa
#'
#' @details
#'
#' Uses R's L-BFGS-B (\code{\link{optim}}) for the outer problem and BFGS
#' \code{\link[n1qn1]{n1qn1}} (restoring the prior individual Hessian) for
#' the inner problem, which is left unscaled since eta estimates start near
#' zero. The covariance step is performed on the unscaled problem, so its
#' condition number may differ from the scaled problem's.
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
foceiControl <- function(sigdig = 4, #
                         ...,
                         epsilon = NULL, # 1e-4,
                         maxInnerIterations = 1000, #
                         maxOuterIterations = 5000, #
                         n1qn1nsim = NULL, #
                         print = 1L, #
                         printNcol = NULL, #
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
                         covType = c("analytic", "fd"), #
                         covSolveTol = NULL, #
                         covFull = TRUE, #
                         # norm of weights = 1/0.225
                         #hessEps = (1/0.225*.Machine$double.eps)^(1 / 4), #
                         hessEps =(.Machine$double.eps)^(1/3),
                         #hessEpsLlik =(1/0.225*.Machine$double.eps)^(1/4),
                         hessEpsLlik =(.Machine$double.eps)^(1/3),
                         optimHessType = c("central", "forward"),
                         optimHessCovType=c("central", "forward"),
                         eventType = c("central", "forward"), #
                         eventSens = c("jump", "fd"), #
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
                         useColor = NULL, #
                         boundTol = NULL, #
                         calcTables = TRUE,#
                         noAbort = TRUE, #
                         interaction = TRUE, #
                         foce = c("nonmem", "foce+"), #
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
                         outerOpt = c("lbfgsb3c",
                                      "nlminb",
                                      "bobyqa",
                                      "L-BFGS-B",
                                      "mma",
                                      "lbfgsbLG",
                                      "slsqp",
                                      "uobyqa",
                                      "newuoa"), #
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
                         muModel = c("none", "irls", "lin"), #
                         muRefCovAlg = TRUE, #
                         muModelTol = 1e-3, #
                         muModelMaxCycles = 10L, #
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
                         indTolRelax = TRUE, #
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
                         agqHi=Inf,
                         eventSens = c("jump", "fd"),
                         sensMethod = c("default", "forward", "adjoint"),
                         boundedTransform=TRUE) { #
  eventSens <- match.arg(eventSens)
  ## sensMethod: "forward" variational ODE parameter sensitivities; "adjoint"
  ## solves them with the in-engine discrete adjoint (matching s-method);
  ## "default" defers to getOption("nlmixr2est.adjoint").
  sensMethod <- match.arg(sensMethod)
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
  # Print args are absorbed/validated by iterPrintControl(); `iterPrintControl`
  # picked up from `...` handles the round-trip via do.call(foceiControl, .ctl).
  .iterPrintControl <- .absorbIterPrintControl(print = print,
                                               printNcol = printNcol,
                                               useColor = useColor,
                                               iterPrintControl = list(...)$iterPrintControl)
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

  foce <- match.arg(foce)
  ## FOCE (interaction=FALSE) residual-variance choice: 0="nonmem" (eta=0 frozen R),
  ## 1="foce+" (live conditional R).  Ignored when interaction=TRUE (FOCEi).
  foceType <- as.integer(foce == "foce+")

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
  ## How dosing/event-parameter (alag, F, rate, dur, ...) sensitivities are
  ## computed: "fd" (legacy finite differences) or "jump" (analytic jump/event
  ## sensitivities from rxode2).  "fd" is the backward-compatible default.
  eventSens <- match.arg(eventSens)

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
  .covMethodMissing <- missing(covMethod)
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
  covType <- match.arg(covType)
  if (!is.null(covSolveTol)) checkmate::assertNumeric(covSolveTol, len = 1, lower = 0,
                                                      finite = TRUE, any.missing = FALSE)
  checkmate::assertFlag(covFull)
  # covType="analytic" seams the analytic R-matrix into the R step.  Default to "r"
  # (the observed-information R^-1); an explicit "r,s"/"s" is honored -- the analytic R
  # then feeds the native finite-difference sandwich / S-matrix.  A missing method or ""
  # uses "r" so a covariance is always produced.
  if (identical(covType, "analytic")) {
    if (.covMethodMissing) {
      covMethod <- 2L
    } else if (identical(covMethod, 0L)) {
      warning("covType=\"analytic\" computes a covariance; using \"r\" instead of \"\"", call. = FALSE)
      covMethod <- 2L
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
    } else if (outerOpt == "uobyqa") {
      outerOptFun <- .uobyqa
      outerOpt <- -1L
    } else if (outerOpt == "newuoa") {
      outerOptFun <- .newuoa
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

  muModel <- match.arg(muModel)
  checkmate::assertLogical(muRefCovAlg, any.missing=FALSE, len=1)
  checkmate::assertNumeric(muModelTol, lower=0, len=1, any.missing=FALSE)
  checkmate::assertIntegerish(muModelMaxCycles, lower=1, len=1, any.missing=FALSE)
  muModelMaxCycles <- as.integer(muModelMaxCycles)

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
  checkmate::assertLogical(indTolRelax, any.missing=FALSE, len=1)
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
  checkmate::assertLogical(boundedTransform, len=1, any.missing=FALSE)
  .ret <- list(
    maxOuterIterations = as.integer(maxOuterIterations),
    maxInnerIterations = as.integer(maxInnerIterations),
    n1qn1nsim = as.integer(n1qn1nsim),
    iterPrintControl = .iterPrintControl,
    lbfgsLmm = as.integer(lbfgsLmm),
    lbfgsPgtol = as.double(lbfgsPgtol),
    lbfgsFactr = as.double(lbfgsFactr),
    scaleTo = scaleTo,
    epsilon = epsilon,
    derivEps = derivEps,
    derivMethod = derivMethod,
    covDerivMethod = covDerivMethod,
    covMethod = covMethod,
    covType = covType,
    covSolveTol = covSolveTol,
    covFull = covFull,
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
    boundTol = as.double(boundTol),
    calcTables = calcTables,
    noAbort = noAbort,
    interaction = interaction,
    foce = foce,
    foceType = foceType,
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
    muModel = muModel,
    muRefCovAlg = muRefCovAlg,
    muModelTol = as.double(muModelTol),
    muModelMaxCycles = muModelMaxCycles,
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
    indTolRelax = as.logical(indTolRelax),
    eventType = eventType,
    eventSens = eventSens,
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
    agqLow=as.double(agqLow),
    eventSens=eventSens,
    sensMethod=sensMethod,
    boundedTransform=boundedTransform
  )
  if (!is.null(.xtra$est)) {
    .ret$est <- .xtra$est
  }
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

.rxUiDeparseFoceiControl <- function(object, var, type="foceiControl") {
  .ret <- eval(str2lang(paste0(type, "()")))
  .outerOpt <- character(0)
  if (object$outerOpt == -1L && object$outerOptTxt == "custom") {
    warning("functions for `outerOpt` cannot be deparsed, reset to default",
            call.=FALSE)
  } else if (!(object$outerOptTxt %in% c(.ret$outerOptTxt, "stats::optimize"))) {
    .outerOpt <- paste0("outerOpt = ", deparse1(object$outerOptTxt))
  }
  .w <- .deparseDifferent(.ret, object, .foceiControlInternal)
  if (length(.w) == 0 && length(.outerOpt) == 0) {
    return(str2lang(paste0(var, " <- ", type, "()")))
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
      paste0("innerOpt = ", deparse1(names(.innerOptFun[which(object[[x]] == .innerOptFun)])))
    } else if (x %in% c("optimHessType", "optimHessCovType")) {
      .methodIdx <- c("central" = 1L, "forward" = 3L)
      paste0(x, " = ", deparse1(names(.methodIdx[which(object[[x]] == .methodIdx)])))
    } else if (x == "eventType") {
      .methodIdx <- c("central" = 2L, "forward" = 3L)
      paste0(x, " = ", deparse1(names(.methodIdx[which(object[[x]] == .methodIdx)])))
    } else if (x %in% c("derivMethod", "covDerivMethod")) {
      .methodIdx <- c("forward" = 0L, "central" = 1L, "switch" = 3L)
      paste0(x, " = ", deparse1(names(.methodIdx[which(object[[x]] == .methodIdx)])))
    } else if (x == "covMethod") {
      if (object[[x]] == 0L) {
        paste0(x, " = \"\"")
      } else {
        .covMethodIdx <- c("r,s" = 1L, "r" = 2L, "s" = 3L)
        paste0(x, " = ", deparse1(names(.covMethodIdx[which(object[[x]] == .covMethodIdx)])))
      }
    } else {
      paste0(x, " = ", deparse1(object[[x]]))
    }
  }, character(1)), .outerOpt)
  str2lang(paste(var, " <- ", type, "(", paste(.retD, collapse=", "), ")"))
}

#' @export
rxUiDeparse.foceiControl <- function(object, var) {
  .rxUiDeparseFoceiControl(object, var, type="foceiControl")
}
