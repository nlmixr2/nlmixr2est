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

test_that(".impCovNatural() maps a covariance with no estimated Omega element", {
  .v <- matrix(c(0.04, 0.01, 0.01, 0.09), 2)
  .ini <- data.frame(name = c("tka", "tcl"), ntheta = 1:2, neta1 = NA_integer_, neta2 = NA_integer_, fix = FALSE)
  for (.om in list(NULL, matrix(0, 0, 0))) {
    .r <- .impCovNatural(.v, 1:2, list(), .om, c("tka", "tcl"), character(0), .ini)
    expect_identical(.r$cov, matrix(.v, 2, dimnames = list(c("tka", "tcl"), c("tka", "tcl"))))
    expect_identical(unname(.r$jacobian), diag(2))
  }
  # every Omega element fixed: only the thetas are mapped
  .iniFix <- rbind(.ini, data.frame(name = "eta.ka", ntheta = NA_integer_, neta1 = 1L, neta2 = 1L, fix = TRUE))
  .r <- .impCovNatural(.v, 1:2, list(), matrix(0.4), c("tka", "tcl"), "eta.ka", .iniFix)
  expect_identical(.r$cov, matrix(.v, 2, dimnames = list(c("tka", "tcl"), c("tka", "tcl"))))
  # Omega elements the count does not cover are still refused
  .ini1 <- replace(.iniFix, "fix", FALSE)
  expect_identical(
    .impCovNatural(.v, 1:2, list(), matrix(0.4), c("tka", "tcl"), "eta.ka", .ini1),
    "could not be mapped to the Omega variances and covariances"
  )
})

test_that(".impCovNatural() maps only the estimated thetas and Omega elements", {
  # tcl and eta.cl fixed: the covariance has rows for tka, tv and the Omega
  # parameters of eta.ka and eta.v.  With a diagonal Omega and a square-root
  # diagonal, Omega_ii = 1 / p_i^2, so d(Omega_ii)/d(p_i) = -2 / p_i^3.
  .eta <- c("eta.ka", "eta.cl", "eta.v")
  .ini <- data.frame(
    name = c("tka", "tcl", "tv", .eta),
    ntheta = c(1L, 2L, 3L, NA, NA, NA),
    neta1 = c(NA, NA, NA, 1L, 2L, 3L),
    neta2 = c(NA, NA, NA, 1L, 2L, 3L),
    fix = c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE)
  )
  .om <- diag(c(0.4, 0.07, 0.02))
  .p <- 1 / sqrt(c(0.4, 0.02))
  .dOm <- list(diag(c(-2 / .p[1]^3, 0, 0)), diag(c(0, 0, -2 / .p[2]^3)))
  .v <- crossprod(matrix(c(3, 1, 0.2, 0.5, 0.1, 0.3, 1, 2, 0.4, 0.2, 0.1, 0.6, 1, 0.3, 0.2, 0.7), 4, 4)) + diag(4)
  .r <- .impCovNatural(.v, c(1L, 3L), .dOm, .om, c("tka", "tcl", "tv"), .eta, .ini)
  .nm <- c("tka", "tv", "om.eta.ka", "om.eta.v")
  expect_identical(dimnames(.r$cov), list(.nm, .nm))
  .j <- diag(c(1, 1, -2 / .p^3))
  expect_identical(unname(.r$jacobian), .j)
  expect_equal(unname(.r$cov), .j %*% .v %*% t(.j))
})

# A short imp/impmap control for the fits below
.impCovCtl <- function(covMethod) {
  impmapControl(print = 0L, nIter = 5L, isample = 100L, covMethod = covMethod, calcTables = FALSE)
}

# A mock fit environment for .impCovInstall(): two thetas and a 3 x 3 Omega
# with one off-diagonal element, so six estimated parameters
.impMockEnv <- function() {
  .env <- new.env(parent = emptyenv())
  .env$thetaNames <- c("tka", "tcl")
  .env$etaNames <- c("eta.ka", "eta.cl", "eta.v")
  .env$ui <- list(
    iniDf = data.frame(
      name = c("tka", "tcl", .env$etaNames, "(eta.v,eta.cl)"),
      ntheta = c(1L, 2L, NA, NA, NA, NA),
      neta1 = c(NA, NA, 1L, 2L, 3L, 3L),
      neta2 = c(NA, NA, 1L, 2L, 3L, 2L),
      fix = FALSE
    )
  )
  .env
}

.impMockOmega <- matrix(c(0.4, 0, 0, 0, 0.07, 0.01, 0, 0.01, 0.02), 3)

# An information matrix with eigenvalues ev, on fixed eigenvectors
.impMockInfo <- function(ev) {
  .q <- qr.Q(qr(matrix(c(3, 1, 0.2, 0.5, 0.1, 0.3, 1, 2, 0.4, 0.2, 0.1, 0.6), 6, 6) + diag(6)))
  .q %*% diag(ev) %*% t(.q)
}

test_that(".impCovInstall() installs a positive-definite imp covariance as \"imp\"", {
  .rx <- rxode2::rxSymInvCholCreate(mat = .impMockOmega, diag.xform = "sqrt")
  .dOm <- lapply(.rx$d.omegaInv, function(.d) -.impMockOmega %*% .d %*% .impMockOmega)
  .info <- .impMockInfo(c(40, 25, 9, 4, 2, 1))
  .env <- .impMockEnv()
  expect_no_warning(
    .ok <- .impCovInstall(.env, solve(.info), 1:2, .dOm, .impMockOmega, .rx$theta, .info)
  )
  expect_true(.ok)
  expect_identical(.env$covMethod, "imp")
  .j <- .env$impCovJacobian
  expect_equal(.env$cov, .j %*% solve(.info) %*% t(.j), ignore_attr = TRUE)
  expect_equal(.env$impCov, .env$cov)
})

test_that(".impCovInstall() repairs an imp information that is not positive definite as \"|imp|\"", {
  .rx <- rxode2::rxSymInvCholCreate(mat = .impMockOmega, diag.xform = "sqrt")
  .dOm <- lapply(.rx$d.omegaInv, function(.d) -.impMockOmega %*% .d %*% .impMockOmega)
  .ev <- c(40, 25, 9, 4, 2, -0.5)
  .info <- .impMockInfo(.ev)
  .env <- .impMockEnv()
  expect_warning(
    .ok <- .impCovInstall(.env, solve(.info), 1:2, .dOm, .impMockOmega, .rx$theta, .info),
    "\"imp\" covariance not positive definite, corrected by sqrtm(imp %*% imp) and installed as \"|imp|\"",
    fixed = TRUE
  )
  expect_true(.ok)
  expect_identical(.env$covMethod, "|imp|")
  .j <- .env$impCovJacobian
  # the raw matrix is kept, indefinite, as $impCov
  expect_equal(.env$impCov, .j %*% solve(.info) %*% t(.j), ignore_attr = TRUE)
  expect_lt(min(eigen(.env$impCov, symmetric = TRUE, only.values = TRUE)$values), 0)
  # the installed one inverts the information with its eigenvalues made positive
  .abs <- .impMockInfo(abs(.ev))
  expect_equal(.env$cov, .j %*% solve(.abs) %*% t(.j), ignore_attr = TRUE, tolerance = 1e-10)
  expect_identical(
    dimnames(.env$cov),
    rep(list(c("tka", "tcl", "om.eta.ka", "om.eta.cl", "cov.eta.v.eta.cl", "om.eta.v")), 2)
  )
  expect_gt(min(eigen(.env$cov, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that(".impCovInstall() installs nothing when the information cannot be repaired", {
  .rx <- rxode2::rxSymInvCholCreate(mat = .impMockOmega, diag.xform = "sqrt")
  .dOm <- lapply(.rx$d.omegaInv, function(.d) -.impMockOmega %*% .d %*% .impMockOmega)
  .info <- .impMockInfo(c(40, 25, 9, 4, 2, -0.5))
  .bad <- replace(.info, 1L, NA_real_)
  .zero <- .impMockInfo(c(40, 25, 9, 4, 2, 0))
  for (.i in list(NULL, .bad, .zero)) {
    .env <- .impMockEnv()
    expect_warning(
      .ok <- .impCovInstall(.env, solve(.info), 1:2, .dOm, .impMockOmega, .rx$theta, .i),
      "\"imp\" covariance is not positive definite; none installed",
      fixed = TRUE
    )
    expect_false(.ok)
    expect_null(.env$cov)
    expect_null(.env$covMethod)
    expect_true(is.matrix(.env$impCov))
  }
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
      .imp <- .nlmixr(.impCovModel, theo_sd, .est, .impCovCtl("imp"))
      .none <- .nlmixr(.impCovModel, theo_sd, .est, .impCovCtl(""))
      expect_true(is.environment(.imp$env) && exists("impCovThetaN", envir = .imp$env))
      expect_identical(.imp$llikObs, .none$llikObs, label = .est)
    }
  })

  test_that("the imp covariance reports Omega on the variance-covariance scale", {
    # The finite-difference Hessian is taken over the parameters impmap estimates
    # Omega in (chol(Omega^-1), sqrt diagonal); its om.<eta>/cov.<eta>.<eta> rows
    # are mapped to the Omega elements they are named after.
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
})
