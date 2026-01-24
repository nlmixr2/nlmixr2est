## saem_fit.R: population PK/PD modeling library
##
## Copyright (C) 2014 - 2016  Wenping Wang
##
## This file is part of nlmixr2.
##
## nlmixr2 is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## nlmixr2 is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with nlmixr2.  If not, see <http://www.gnu.org/licenses/>.

#' Configure an SAEM model
#'
#' Configure an SAEM model by generating an input list to the SAEM model function
#'
#' @param model a compiled saem model
#' @param data input data
#' @param inits initial values
#' @param mcmc a list of various mcmc options
#' @param ODEopt optional ODE solving options
#' @param seed seed for random number generator
#' @param distribution one of c("normal","poisson","binomial")
#' @param fixed a character vector of fixed effect only parameters (no random effects attached) to be fixed
#' @param DEBUG Integer determining if debugging is enabled
#' @param type indicates the type of optimization for the residuals; Can be one of c("nelder-mead", "newuoa")
#' @param lambdaRange This indicates the range that Box-Cox and Yeo-Johnson parameters are constrained to be;  The default is 3 indicating the range (-3,3)
#' @param powRange This indicates the range that powers can take for residual errors;  By default this is 10 indicating the range is c(1/10, 10) or c(0.1,10)
#' @inheritParams saemControl
#'
#' @return Returns a list neede for the saem fit procedure
#'
#' @details
#'    Fit a generalized nonlinear mixed-effect model by he Stochastic
#'    Approximation Expectation-Maximization (SAEM) algorithm
#'
#' @author Wenping Wang & Matthew Fidler
#' @examples
#' \donttest{
#'
#' # In this ODE system we simply specify the ODEs
#'
#' ode <- "d/dt(depot) =-KA*depot;
#'         d/dt(centr) = KA*depot - KE*centr;"
#' m1 <- rxode2(ode)
#'
#'
#' # In this ode System, we also specify the concentration as C2 = centr/V
#'
#' ode <- "C2 = centr/V;
#'       d/dt(depot) =-KA*depot;
#'       d/dt(centr) = KA*depot - KE*centr;"
#' m2 = rxode2(ode)
#'
#' PKpars <- function() {
#'   CL <- exp(lCL)
#'   V <- exp(lV)
#'   KA <- exp(lKA)
#'   KE <- CL / V
#' }
#'
#' PRED <- function() centr / V
#' PRED2 <- function() C2
#'
#' # You can also use the nlmixr2 UI to run this model and call the lower level functions
#'
#' one.compartment <- function() {
#' ini({
#'   tka <- 0.45 # Log Ka
#'   tcl <- 1 # Log Cl
#'   tv <- 3.45    # Log V
#'   eta.ka ~ 0.6
#'   eta.cl ~ 0.3
#'   eta.v ~ 0.1
#'   add.sd <- 0.7
#'   wt.est <- 0.0
#' })
#' model({
#'   ka <- exp(tka + eta.ka)
#'   cl <- exp(tcl + eta.cl)
#'   v <- exp(tv + eta.v + wt.est * WT)
#'   d/dt(depot) = -ka * depot
#'   d/dt(center) = ka * depot - cl / v * center
#'   cp = center / v
#'   cp ~ add(add.sd)
#' })
#' }
#' fit  <- nlmixr2(one.compartment, theo_sd, "saem")
#' fit
#'
#' }
#' @noRd
.configsaem <- function(model, data, inits,
                       mcmc = list(niter = c(200, 300), nmc = 3, nu = c(2, 2, 2)),
                       rxControl = list(atol = 1e-6, rtol = 1e-4, method = "lsoda", maxeval = 100000),
                       distribution = c("normal", "poisson", "binomial"),
                       seed = 99, fixedOmega = NULL, fixedOmegaValues=NULL,
                       parHistThetaKeep=NULL,
                       parHistOmegaKeep=NULL,
                       DEBUG = 0,
                       tol = 1e-4, itmax = 100L, type = c("nelder-mead", "newuoa"),
                       lambdaRange = 3, powRange = 10,
                       odeRecalcFactor=10^(0.5),
                       maxOdeRecalc=5L, nres,
                       perSa=0.75,
                       perNoCor=0.75,
                       perFixOmega=0.5,
                       perFixResid=0.75,
                       resFixed,
                       ue) {
  if (is.null(fixedOmega)) stop("requires fixedOmega", call.=FALSE)
  if (is.null(fixedOmegaValues)) stop("requires fixedOmegaValues", call.=FALSE)
  if (is.null(parHistThetaKeep)) stop("requires parHistThetaKeep", call.=FALSE)
  if (is.null(parHistOmegaKeep)) stop("requires parHistOmegaKeep", call.=FALSE)
  type.idx <- c("nelder-mead" = 1L, "newuoa" = 2L)
  type <- match.arg(type)
  type <- type.idx[type]
  force(rxControl)
  .env <- nlmixr2global$nlmixrEvalEnv$envir
  if (!is.environment(.env)) {
    .env <- parent.frame(1)
  }
  rxControl <- do.call(rxode2::rxControl, rxControl)
  rxControl$envir <- .env
  set.seed(seed)
  distribution.idx <- c("normal" = 1, "poisson" = 2, "binomial" = 3)
  distribution <- match.arg(distribution)
  distribution <- distribution.idx[distribution]
  .data <- data
  ## rxode2::rxTrans(data, model)
  data <- list(nmdat = data)

  neq <- attr(model$saem_mod, "neq")
  nlhs <- attr(model$saem_mod, "nlhs")
  inPars <- attr(model$saem_mod, "inPars")
  ninputpars <- length(inPars)
  opt <- optM <- c(list(neq = neq, nlhs = nlhs, inits = numeric(neq)),
    ninputpars = ninputpars, inPars = inPars
  )

  model$N.eta <- attr(model$saem_mod, "nrhs")
  model$nendpnt <- attr(model$saem_mod, "nendpnt")
  if (is.null(model$nendpnt)) {
    model$nendpnt <- 1
  }
  if (is.null(model$log.eta)) {
    model$log.eta <- rep(TRUE, model$N.eta)
  }
  if (is.null(model$omega)) {
    model$omega <- diag(model$N.eta)
  }
  if (is.null(model$res.mod)) {
    model$res.mod <- rep(1, model$nendpnt)
  }
  if (is.null(inits$omega)) {
    inits$omega <- rep(1, model$N.eta) * 4
  }
  if (is.null(inits$ares)) {
    inits$ares <- 10
  }
  if (is.null(inits$bres)) {
    inits$bres <- 1
  }
  if (is.null(inits$cres)) {
    inits$cres <- 1
  }
  if (is.null(inits$lres)) {
    inits$lres <- 1
  }
  if (is.null(mcmc$print)) {
    mcmc$print <- 1
  }
  if (model$N.eta - length(inits$theta) > 0) {
    inits$theta <- c(inits$theta, rep(NA_real_, model$N.eta - length(inits$theta)))
  }
  if (is.null(names(inits$theta))) {
    names(inits$theta) <- rep("", length(inits$theta))
  }
  .nt <- names(inits$theta)
  .nt <- .nt[!is.na(.nt)]
  inits.save <- inits
  inits$theta.fix <- matrix(names(inits$theta),
    byrow = TRUE,
    ncol = model$N.eta
    )
  inits$theta <- matrix(inits$theta, byrow = TRUE, ncol = model$N.eta)
  model$cov.mod <- 1 - is.na(inits$theta)
  data$N.covar <- nrow(inits$theta) - 1
  inits$theta[is.na(inits$theta)] <- 0

  ###  FIXME
  mcmc$stepsize <- 0:1
  mcmc$burn.in <- 300

  ###  FIXME: chk covars as char vec
  wh <- setdiff(c(model$covars, inPars), names(data$nmdat))
  if (length(wh)) {
    msg <- paste0("covariate(s) not found: ", paste(wh, collapse = ", "))
    stop(msg)
  }
  s <- subset(data$nmdat, EVID == 0)
  data$data <- as.matrix(s[, c("ID", "TIME", "DV", c(model$covars, inPars))])

  ###  chk for no obs records
  wh <- setdiff(unique(data$nmdat$ID), unique(data$data[, "ID"]))
  if (length(wh)) {
    msg <- paste0("No data with ID: ", paste(wh, collapse = ", "))
    stop(msg)
  }

  nphi <- model$N.eta
  mcov <- model$cov.mod
  covstruct <- model$omega

  check <- sum((covstruct - t(covstruct)) != 0)
  if (check) stop("illegal covstruct")
  check <- nphi - dim(covstruct)[1]
  if (check) stop("nphi and covstruct dim mismatch")

  check <- prod(mcov[1, ])
  if (check == 0) {
    print(mcov)
    stop("structural parameter(s) absent", call.=FALSE)
  }
  check <- nphi - dim(mcov)[2]
  if (check) stop("nphi and ncol(mcov) mismatch")
  check <- sum(dim(inits$theta) - dim(mcov) != 0)
  if (check) stop("initial theta's and mcov dim mismatch")
  check <- data$N.covar + 1 - dim(mcov)[1]
  if (check) stop("dim mcov and N.covar mismatch")

  check <- length(model$log.eta) - nphi
  if (check) stop("jlog length and nphi mismatch")

  check <- length(inits$omega) - nphi
  if (check) stop("length of omega inits and nphi mismatch")

  # check = mcmc$burn.in>sum(mcmc$niter)
  # if (check) stop("#burn-in exceeds niter")

  check <- prod(is.element(covstruct, c(0, 1)))
  if (check == 0) warning("non-zero value(s) in covstruct set to 1")
  covstruct[covstruct != 0] <- 1

  check <- prod(is.element(mcov, c(0, 1)))
  if (check == 0) warning("non-zero value(s) in mcov set to 1")
  mcov[mcov != 0] <- 1

  check <- sum(inits$theta[1, model$log.eta] <= 0)
  if (check) stop("illegal initial theta's")
  check <- sum(inits$omega <= 0)
  if (check) stop("illegal initial omega")
  # check = inits$sigma2<=0
  # if (check) stop("illegal initial sigma2")
  check <- sum(diag(covstruct) == 1)
  if (!check) stop("0 ETA's")
  y <- data$data[, "DV"]
  id <- data$data[, "ID"]
  check <- any(diff(unique(id)) != 1)
  if (check) stop("saem classic UI needs sequential ID. check your data")
  ntotal <- length(id)
  N <- length(unique(id))
  if (is.null(model$covars)) {
    covariables <- NULL
  } else {
    covariables <- unlist(stats::aggregate(.as.data.frame(data$data[, model$covars, drop = FALSE]),
                                           list(id),
                                           unique)[, -1, drop = FALSE])
  }
  if (!is.null(covariables)) {
    if (length(covariables) == N * data$N.covar) {
      dim(covariables) <- c(N, data$N.covar)
    } else {
      message("covariables")
      print(covariables)
      message("covars")
      print(model$covars)
      print(data$N.covar)
      stop("internal covariate mismatch for 'saem'",
           call.=FALSE)
    }
  }
  nb_measures <- table(id)
  ncov <- data$N.covar + 1
  nmc <- mcmc$nmc
  nM <- mcmc$nmc * N
  yM <- rep(y, nmc)
  mlen <- max(nb_measures)
  io <- t(sapply(nb_measures, function(x) rep(1:0, c(x, mlen - x))))
  ix <- rep(1:dim(io)[1], nmc)
  ioM <- io[ix, ]
  indioM <- grep(1, t(ioM)) - 1
  ## mPars <- if (ninputpars == 0) NULL else unlist(stats::aggregate(.as.data.frame(data$data[, inPars]), list(id), unique)[, -1])
  ## if (!is.null(mPars)) {
  ##   dim(mPars) <- c(N, ninputpars)
  ##   opt$mPars <- mPars
  ##   ix <- rep(1:dim(mPars)[1], nmc)
  ##   optM$mPars <- mPars[ix, ]
  ##   dim(optM$mPars) <- c(nmc * N, ninputpars)
  ## }

  if (is.null(data$nmdat$CMT)) data$nmdat$CMT <- 1 ## CHECKME
  if (any(is.na(data$nmdat$CMT))) {
    stop("'CMT' has NA(s)")
  }
  ## CHECKME
  form <- attr(model$saem_mod, "form")
  .nobs <- 0
  dat <- rxode2::etTrans(data$nmdat, attr(model$saem_mod, "rx"), addCmt=TRUE, dropUnits=TRUE, allTimeVar=TRUE,
                         addlKeepsCov = rxControl$addlKeepsCov, addlDropSs = rxControl$addlDropSs,
                         ssAtDoseTime = rxControl$ssAtDoseTime)
  .nobs <- attr(class(dat), ".rxode2.lst")$nobs
  dat <- as.data.frame(dat) # convert back evid=3 oddness...
  ## if(length(dat) !=7) stop("SAEM doesn't support time varying covariates yet.");
  .rx <- attr(model$saem_mod, "rx")
  .pars <- .rx$params
  .pars <- setNames(rep(1.1, length(.pars)), .pars)
  .pars <- .pars[is.na(match(names(.pars), inPars))]
  opt$.rx <- .rx
  opt$.pars <- .pars
  ## opt$.dat <- dat;
  dat <- .as.data.frame(dat[, -6])
  # Optimized: Create lookup table before vapply call for O(1) access
  .inParsLookup <- setNames(rep(TRUE, length(inPars)), inPars)
  names(dat) <- vapply(names(dat), function(n) {
    if (isTRUE(.inParsLookup[[n]])) return(n)
    return(toupper(n))
  }, character(1), USE.NAMES = FALSE)

  dat$ID <- as.integer(dat$ID)

  evt <- dat
  evt$ID <- evt$ID - 1
  ## r
  evtM <- evt[rep(1:dim(evt)[1], nmc), ]
  evtM$ID <- cumsum(c(FALSE, diff(evtM$ID) != 0))

  # i1:
  i1 <- grep(1, diag(covstruct))
  i0 <- grep(0, diag(covstruct))
  nphi1 <- sum(diag(covstruct))
  nphi0 <- nphi - nphi1
  na <- length(mcmc$stepsize)
  nlambda1 <- sum(mcov[, i1])
  nlambda0 <- sum(mcov[, i0])
  nlambda <- nlambda1 + nlambda0
  nd1 <- nphi1 + nlambda1 + 1
  nd2 <- nphi1 + nlambda1 + nlambda0
  nb_param <- nd2 + 1
  Mcovariables <- cbind(rep(1, N), covariables)[, 1:nrow(mcov)]
  dim(Mcovariables) <- c(length(Mcovariables) / nrow(mcov), nrow(mcov)) # FIXME

  # get fixed ix
  fixed <- inits$theta.fix != ""
  wh <- fixed[, i1][mcov[, i1] == 1]
  len <- length(wh)
  fixed.i1 <- (1:len)[wh] - 1
  wh <- fixed[, i0][mcov[, i0] == 1]
  len <- length(wh)
  fixed.i0 <- (1:len)[wh] - 1

  jlog1 <- grep(TRUE, model$log.eta)
  jcov <- grep(TRUE, apply(mcov, 1, sum) > 0)
  covstruct1 <- covstruct[i1, i1]
  dim(covstruct1) <- c(nphi1, nphi1)
  ind_cov <- grep(1, mcov[mcov > 0])

  mcov1 <- matrix(mcov[, i1], ncol = length(i1))
  mcov0 <- matrix(mcov[, i0], nrow = nrow(mcov), ncol = length(i0))
  ind_cov1 <- grep(1, mcov1[mcov1 > 0]) - 1
  ind_cov0 <- grep(1, mcov0[mcov0 > 0]) - 1

  pc <- apply(mcov, 2, sum)
  ipc <- cumsum(c(0, pc[1:(nphi - 1)])) + 1
  ipcl1 <- ipc[jlog1]
  for (x in jlog1) inits$theta[1, x] <- log(inits$theta[1, x])

  idx <- as.vector(mcov1 > 0)
  COV1 <- Mcovariables[, row(mcov1)[idx]]
  dim(COV1) <- c(N, sum(idx))
  COV21 <- crossprod(COV1)

  x <- mcov1 * col(mcov1)
  x <- sapply(x[idx], function(x) {
    ret <- rep(0, nphi1)
    ret[x] <- 1
    ret
  })
  LCOV1 <- t(x)
  dim(LCOV1) <- c(nlambda1, nphi1)
  pc1 <- apply(LCOV1, 2, sum)

  x1 <- diag(sum(idx))
  diag(x1) <- inits$theta[, i1][idx]
  MCOV1 <- x1 %*% LCOV1
  jcov1 <- grep(1, LCOV1) - 1

  idx <- as.vector(mcov0 > 0)
  COV0 <- Mcovariables[, row(mcov0)[idx]]
  dim(COV0) <- c(N, sum(idx))
  COV20 <- crossprod(COV0)

  x <- mcov0 * col(mcov0)
  x <- sapply(x[idx], function(x) {
    ret <- rep(0, nphi0)
    ret[x] <- 1
    ret
  })
  LCOV0 <- t(x)
  dim(LCOV0) <- c(nlambda0, nphi0)

  x1 <- diag(sum(idx))
  diag(x1) <- inits$theta[, i0][idx]
  if (dim(x1)[1] > 0) {
    MCOV0 <- x1 %*% LCOV0
  } else {
    MCOV0 <- matrix(x1, nrow = 0, ncol = dim(LCOV0)[2])
  }
  jcov0 <- grep(1, LCOV0) - 1

  mprior_phi1 <- Mcovariables %*% inits$theta[, i1, drop = FALSE]
  mprior_phi0 <- Mcovariables %*% inits$theta[, i0, drop = FALSE]

  Gamma2_phi1 <- diag(nphi1)
  diag(Gamma2_phi1) <- inits$omega[i1]
  Gamma2_phi0 <- diag(nphi0)
  diag(Gamma2_phi0) <- inits$omega[i0]

  Gamma2_phi1fixedIx <- fixedOmega[i1, i1, drop = FALSE]
  Gamma2_phi1fixedValues <- fixedOmegaValues[i1, i1, drop = FALSE]
  Gamma2_phi1fixed <- as.integer(any(Gamma2_phi1fixedIx == 1L))

  phiM <- matrix(0, N, nphi)
  phiM[, i1] <- mprior_phi1
  phiM[, i0] <- mprior_phi0
  phiM <- phiM[rep(1:N, nmc), , drop = FALSE]
  .tmp <- diag(sqrt(inits$omega))
  if (model$N.eta == 1) .tmp <- matrix(sqrt(inits$omega))
  .dim <- dimnames(ue)[[2]]
  .ue <- do.call("cbind",
                 lapply(names(model$log.eta),
                        function(n) {
                          if (n %in% .dim) return(ue[, n])
                          rep(1L, length(ue[, 1]))
                        }))
  dimnames(.ue) <- list(NULL, names(model$log.eta))
  .mat2 <- matrix(rnorm(phiM), dim(phiM))
  .ue <- .ue[rep(1:N, nmc),, drop = FALSE] * 1.0
  .mat2 <- .mat2 * .ue
  phiM <- phiM + .mat2 %*% .tmp
  # now replace with what is needed inside saem sampling

  # Since the .mat2 is adjusted for uninformative etas, the phiM stats
  # do not need to be adjusted
  mc.idx <- rep(1:N, nmc)
  statphi <- sapply(1:nphi, function(x) {
    tapply(phiM[, x], mc.idx, mean)
  })
  statphi11 <- statphi[, i1]
  dim(statphi11) <- c(N, length(i1))
  statphi01 <- statphi[, i0]
  dim(statphi01) <- c(N, length(i0))
  statphi12 <- crossprod(phiM[, i1])
  statphi02 <- crossprod(phiM[, i0])

  # x = mcov *cumsum(mcov)
  # x1 = cbind(x[,i1], x[,i0])
  # indiphi = order(x1[x1>0])	#FINDME

  niter <- sum(mcmc$niter)
  niter_phi0 <- round(niter * .5)
  nb_sa <- round(mcmc$niter[1] * perSa)
  nb_correl <- round(mcmc$niter[1] * perNoCor)
  nb_fixOmega <- round(mcmc$niter[1] * perFixOmega)
  nb_fixResid <- round(mcmc$niter[1] * perFixResid)
  va <- mcmc$stepsize
  vna <- mcmc$niter
  na <- length(va)
  pas <- 1 / (1:vna[1])^va[1]
  for (ia in 2:na) {
    end <- length(pas)
    k1 <- pas[end]^(-1 / va[ia])
    pas <- c(pas, 1 / ((k1 + 1):(k1 + vna[ia]))^va[ia])
  }
  pash <- c(rep(1, mcmc$burn.in), 1 / (1:niter))
  minv <- rep(1e-20, nphi)

  # preserve par order when printing iter history
  mcov[mcov == 1] <- 1:nlambda
  ilambda1 <- mcov[, i1]
  ilambda1 <- ilambda1[ilambda1 > 0] - 1
  ilambda0 <- mcov[, i0]
  ilambda0 <- ilambda0[ilambda0 > 0] - 1

  i1 <- i1 - 1
  i0 <- i0 - 1
  opt$distribution <- distribution
  opt$paramUpdate <- attr(model$saem_mod, "paramUpdate")
  optM$paramUpdate <- attr(model$saem_mod, "paramUpdate")
  opt$rxControl <- rxControl
  optM$rxControl <- rxControl
  cfg <- list(
    rxControl = rxControl,
    ue=.ue,
    inits = inits.save,
    nu = mcmc$nu,
    niter = niter,
    nb_sa = nb_sa,
    nb_correl = nb_correl,
    nb_fixOmega=nb_fixOmega,
    nb_fixResid=nb_fixResid,
    niter_phi0 = niter_phi0,
    nmc = nmc,
    coef_phi0 = .9638, # FIXME
    rmcmc = .5,
    coef_sa = .95,
    pas = pas,
    pash = pash,
    minv = minv,
    N = N,
    ntotal = ntotal,
    y = y,
    yM = yM,
    phiM = phiM,
    evt = as.matrix(evt),
    evtM = as.matrix(evtM),
    mlen = mlen,
    indioM = indioM,

    pc1 = pc1,
    covstruct1 = covstruct1,
    Mcovariables = Mcovariables,

    i1 = i1,
    i0 = i0,
    nphi1 = nphi1,
    nphi0 = nphi0,
    nlambda1 = nlambda1,
    nlambda0 = nlambda0,
    COV0 = COV0,
    COV1 = COV1,
    COV20 = COV20,
    COV21 = COV21,
    LCOV0 = LCOV0,
    LCOV1 = LCOV1,
    MCOV0 = MCOV0,
    MCOV1 = MCOV1,
    Gamma2_phi0 = Gamma2_phi0,
    Gamma2_phi1 = Gamma2_phi1,
    Gamma2_phi1fixed=Gamma2_phi1fixed,
    Gamma2_phi1fixedIx=Gamma2_phi1fixedIx,
    Gamma2_phi1fixedValues=Gamma2_phi1fixedValues,
    mprior_phi0 = mprior_phi0,
    mprior_phi1 = mprior_phi1,
    jcov0 = jcov0,
    jcov1 = jcov1,
    ind_cov0 = ind_cov0,
    ind_cov1 = ind_cov1,
    statphi11 = statphi11,
    statphi01 = statphi01,
    statphi02 = statphi02,
    statphi12 = statphi12,
    res.mod = model$res.mod,
    ares = inits$ares,
    bres = inits$bres,
    cres = inits$cres,
    lres = inits$lres,
    opt = opt,
    optM = optM,
    print = mcmc$print,
    distribution = distribution,
    parHistThetaKeep=parHistThetaKeep,
    parHistOmegaKeep=parHistOmegaKeep,
    seed = seed,
    fixed.i1 = fixed.i1,
    fixed.i0 = fixed.i0,
    ilambda1 = as.integer(ilambda1),
    ilambda0 = as.integer(ilambda0),
    nobs = .nobs,
    resFixed=resFixed)

  ## CHECKME
  s <- cfg$evt[cfg$evt[, "EVID"] == 0, "CMT"]
  cfg$opt$cmt_endpnt <- cfg$optM$cmt_endpnt <- sort(unique(s))
  cfg$nendpnt <- length(unique(s))
  if (model$nendpnt != cfg$nendpnt) {
    msg <- sprintf("mis-match in nbr endpoints in model & in data")
    stop(msg)
  }
  t <- unlist(split(1L:length(s), s))
  cfg$ysM <- rep(cfg$y[t], cfg$nmc)
  cfg$ix_sorting <- t - 1 # c-index for sorting by endpnt
  cfg$y_offset <- c(0, cumsum(table(s)))
  s <- cfg$evtM[cfg$evtM[, "EVID"] == 0, "CMT"]
  cfg$ix_endpnt <- as.integer(as.factor(s)) - 1 # to derive vecares & vecbres
  s <- cfg$evtM[cfg$evtM[, "EVID"] == 0, "ID"]
  t <- cumsum(c(0, table(s)))
  cfg$ix_idM <- cbind(t[-length(t)], t[-1] - 1) # c-index of obs records of each subject

  cfg$ares <- rep(10, cfg$nendpnt)
  cfg$bres <- rep(1, cfg$nendpnt)
  cfg$cres <- rep(1, cfg$nendpnt)
  cfg$lres <- rep(1, cfg$nendpnt)
  cfg$yj <- rep(2L, cfg$nendpnt)
  cfg$propT <- rep(0L, cfg$nendpnt)
  cfg$lambda <- rep(1.0, cfg$nendpnt)
  cfg$low <- rep(0.0, cfg$nendpnt)
  cfg$hi <- rep(1.0, cfg$nendpnt)
  cfg$ares[cfg$res.mod == 2] <- 0
  cfg$bres[cfg$res.mod == 1] <- 0
  cfg$res_offset <- cumsum(c(0L, nres))
  cfg$par.hist <- matrix(0, cfg$niter, sum(parHistThetaKeep) + sum(parHistOmegaKeep) + sum(1L - resFixed))

  cfg$DEBUG <- cfg$opt$DEBUG <- cfg$optM$DEBUG <- DEBUG
  cfg$phiMFile <- tempfile("phi-", rxode2::rxTempDir(), ".phi")
  cfg$tol <- tol
  cfg$itmax <- itmax
  cfg$type <- type
  cfg$lambdaRange <- lambdaRange
  cfg$powRange <- powRange
  cfg$odeRecalcFactor <- odeRecalcFactor
  cfg$maxOdeRecalc <- maxOdeRecalc
  cfg
}

#' Print an SAEM model fit summary
#'
#' Print an SAEM model fit summary
#'
#' @param object a saemFit object
#' @param ... others
#' @return a list
#' @export
summary.saemFit <- function(object, ...) {
  fit <- object ## Rcheck hack

  th <- fit$Plambda
  nth <- length(th)
  H <- solve(fit$Ha[1:nth, 1:nth])
  se <- sqrt(diag(H))

  m <- cbind(exp(th), th, se) # FIXME
  ## lhsVars = scan("LHS_VARS.txt", what="", quiet=TRUE)
  ## if (length(lhsVars)==nth) dimnames(m)[[1]] = lhsVars
  dimnames(m)[[2]] <- c("th", "log(th)", "se(log_th)")
  cat("THETA:\n")
  print(m)
  cat("\nOMEGA:\n")
  print(fit$Gamma2_phi1)
  if (any(fit$sig2 == 0)) {
    cat("\nSIGMA:\n")
    print(max(fit$sig2^2))
  } else {
    cat("\nARES & BRES:\n")
    print(fit$sig2)
  }

  invisible(list(theta = th, se = se, H = H, omega = fit$Gamma2_phi1, eta = fit$mpost_phi))
}

#' Print an SAEM model fit summary
#'
#' Print an SAEM model fit summary
#'
#' @param x a saemFit object
#' @param ... others
#' @return a list
#' @export
print.saemFit <- function(x, ...) {
  fit <- x ## Rcheck hack

  th <- fit$Plambda
  nth <- length(th)
  H <- solve(fit$Ha[1:nth, 1:nth])
  se <- sqrt(diag(H))

  m <- cbind(exp(th), th, se) # FIXME
  ## lhsVars = scan("LHS_VARS.txt", what="", quiet=TRUE)
  ## if (length(lhsVars)==nth) dimnames(m)[[1]] = lhsVars
  dimnames(m)[[2]] <- c("th", "log(th)", "se(log_th)")
  cat("THETA:\n")
  print(m)
  cat("\nOMEGA:\n")
  print(fit$Gamma2_phi1)
  if (any(fit$sig2 == 0)) {
    cat("\nSIGMA:\n")
    print(max(fit$sig2^2))
  } else {
    cat("\nARES & BRES:\n")
    print(fit$sig2)
  }

  invisible(list(theta = th, se = se, H = H, omega = fit$Gamma2_phi1, eta = fit$mpost_phi))
}

##' @export
ranef.saemFit <- function(object, ...) {
  object$eta
}

##' @export
fixef.saemFit <- function(object, ...) {
  object$Plambda
}

## FIXME: coef_phi0, rmcmc, coef_sa
## FIXME: Klog, rho, sa, nmc
## FIXME: N.design
## FIXME: g = gc = 1
## FIXME: ODE inits
## FIXME: Tinf for ODE
## FIXME: chk infusion poor fit

## Local Variables:
## ess-indent-level: 2
## indent-tabs-mode: nil
## End:
