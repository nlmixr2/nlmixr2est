# The parameterization impmap estimates Omega in: rxode2's chol(Omega^-1) with a
# transformed diagonal.  Omega as a function of those parameters, for a central
# difference against the analytic Jacobian.
.impOmegaAt <- function(rx, p) {
  rx$theta <- p
  rx$omega
}

.impNumJac <- function(rx, p, pairs, h = 1e-6) {
  vapply(
    seq_along(p),
    function(m) {
      .e <- replace(numeric(length(p)), m, h)
      (.impOmegaAt(rx, p + .e) - .impOmegaAt(rx, p - .e))[pairs] / (2 * h)
    },
    numeric(nrow(pairs))
  )
}

test_that(".impCovNatural() maps chol(Omega^-1) rows to the Omega elements by the delta method", {
  .om <- matrix(c(0.4, 0, 0, 0, 0.07, 0.01, 0, 0.01, 0.02), 3)
  .eta <- c("eta.ka", "eta.cl", "eta.v")
  .ini <- data.frame(
    name = c("tka", "tcl", .eta, "(eta.v,eta.cl)"),
    ntheta = c(1L, 2L, NA, NA, NA, NA),
    neta1 = c(NA, NA, 1L, 2L, 3L, 3L),
    neta2 = c(NA, NA, 1L, 2L, 3L, 2L),
    fix = FALSE
  )
  .pairs <- .foceiOmegaPairs(.om, .ini)
  .v <- crossprod(matrix(c(3, 1, 0.2, 0.5, 0.1, 0.3, 1, 2, 0.4, 0.2, 0.1, 0.6), 6, 6)) + diag(6)
  for (.x in c("sqrt", "log", "identity")) {
    .rx <- rxode2::rxSymInvCholCreate(mat = .om, diag.xform = .x)
    .p <- .rx$theta
    # what impOmegaParDeriv() forms: d(Omega)/dp = -Omega d(Omega^-1)/dp Omega
    .dOm <- lapply(.rx$d.omegaInv, function(.d) -.om %*% .d %*% .om)
    .r <- .impCovNatural(.v, 1:2, .dOm, .om, c("tka", "tcl"), .eta, .ini)
    .nm <- c("tka", "tcl", "om.eta.ka", "om.eta.cl", "cov.eta.v.eta.cl", "om.eta.v")
    expect_identical(dimnames(.r$cov), list(.nm, .nm))
    expect_equal(unname(.r$jacobian[1:2, ]), cbind(diag(2), matrix(0, 2, 4)))
    expect_equal(unname(.r$jacobian[3:6, 1:2]), matrix(0, 4, 2))
    expect_equal(unname(.r$jacobian[3:6, 3:6]), .impNumJac(.rx, .p, .pairs), tolerance = 1e-7, info = .x)
    expect_equal(unname(.r$cov), .r$jacobian %*% .v %*% t(.r$jacobian), ignore_attr = TRUE)
  }
  # a parameter count that does not match the estimated Omega elements is refused
  expect_identical(
    .impCovNatural(.v[1:5, 1:5], 1:2, .dOm[1:3], .om, c("tka", "tcl"), .eta, .ini),
    "could not be mapped to the Omega variances and covariances"
  )
})

nmTest({
  .impCovModel <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("llikObs is that of the estimates, not of the last importance-sampling covariance leg", {
    # impComputeCov() scores every subject at its fixed importance samples, at
    # perturbed parameters, after the final MAP pass set llikObs at the estimates
    for (.est in c("impmap", "imp")) {
      .ctl <- function(covMethod) {
        impmapControl(print = 0L, nIter = 5L, isample = 100L, covMethod = covMethod, calcTables = FALSE)
      }
      .imp <- .nlmixr(.impCovModel, theo_sd, .est, .ctl("imp"))
      .none <- .nlmixr(.impCovModel, theo_sd, .est, .ctl(""))
      expect_true(is.environment(.imp$env) && exists("impCovThetaN", envir = .imp$env))
      expect_identical(.imp$llikObs, .none$llikObs, label = .est)
    }
  })

  test_that("the imp covariance reports Omega on the variance-covariance scale", {
    # The finite-difference Hessian is taken over the parameters impmap estimates
    # Omega in (chol(Omega^-1), sqrt diagonal); those rows were installed as they
    # were, under om.<eta>/cov.<eta>.<eta> names.
    blk <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl + eta.v ~ c(0.3, 0.05, 0.1)
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    f <- .nlmixr(blk, theo_sd, "impmap", impmapControl(print = 0L, nIter = 20L, isample = 200L, calcTables = FALSE))
    expect_identical(f$covMethod, "imp")
    .om <- c("om.eta.ka", "om.eta.cl", "cov.eta.v.eta.cl", "om.eta.v")
    .nm <- c("tka", "tcl", "tv", "add.sd", .om)
    expect_identical(dimnames(f$cov), list(.nm, .nm))
    # the delta method, with the Jacobian of exactly the parameterization the fit
    # used: rxode2 recovers the fit's own Omega parameters from its Omega
    .rx <- rxode2::rxSymInvCholCreate(mat = f$omega, diag.xform = "sqrt")
    expect_equal(.rx$theta, f$env$impCovOmegaPar, tolerance = 1e-8)
    .pairs <- .foceiOmegaPairs(f$omega, f$ui$iniDf)
    .j <- f$env$impCovJacobian
    expect_equal(unname(.j[.om, 5:8]), .impNumJac(.rx, .rx$theta, .pairs), tolerance = 1e-6)
    expect_equal(f$cov, .j %*% f$env$impCovInternal %*% t(.j), ignore_attr = TRUE)
    expect_equal(f$cov[1:4, 1:4], f$env$impCovInternal[1:4, 1:4], ignore_attr = TRUE)
    expect_equal(f$env$impCov, f$cov)
    expect_equal(unname(f$env$impSe), unname(sqrt(diag(f$cov))))
    expect_gt(min(eigen(f$cov, symmetric = TRUE, only.values = TRUE)$values), 0)
  })

  test_that("an imp covariance that is not positive definite is reported, and its stand-in named", {
    # est = "imp" on theo_sd with the ODE one-compartment model: the Monte-Carlo
    # information is indefinite (a negative variance for tka), which was installed
    # as it was.  It is refused, and the FOCEi analytic covariance that the
    # post-fit step installs in its place says it is not the one requested.
    odeM <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .acc <- new.env(parent = emptyenv())
    .acc$w <- character(0)
    f <- withCallingHandlers(
      suppressMessages(nlmixr2(odeM, theo_sd, "imp", impmapControl(print = 0L, calcTables = FALSE))),
      warning = function(w) {
        .acc$w <- c(.acc$w, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    expect_true(any(f$runInfo == "\"imp\" covariance is not positive definite; none installed"))
    expect_true(any(.acc$w == "\"analytic (full)\" covariance installed instead of the requested \"imp\""))
    expect_identical(f$covMethod, "analytic (full)")
    expect_lt(min(eigen(f$env$impCov, symmetric = TRUE, only.values = TRUE)$values), 0)
  })
})
