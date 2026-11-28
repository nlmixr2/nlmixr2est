# est="its" fits; control-only tests are in test-its.R
nmTest({
  skip_on_cran()

  .oneCmt <- function() {
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

  .fitIts <- suppressMessages(suppressWarnings(
    nlmixr2(.oneCmt, nlmixr2data::theo_sd, "its", itsControl(print = 0L))
  ))
  .fitFocei <- suppressMessages(suppressWarnings(
    nlmixr2(.oneCmt, nlmixr2data::theo_sd, "focei", foceiControl(print = 0L))
  ))

  test_that("ITS matches FOCEi on a mu-referenced one-compartment model", {
    expect_identical(.fitIts$env$method, "its")
    expect_true(.fitIts$env$impConverged)
    expect_lt(.fitIts$env$impIter, 50L)
    expect_equal(.fitIts$theta, .fitFocei$theta, tolerance = 0.05)
    expect_equal(diag(.fitIts$omega), diag(.fitFocei$omega), tolerance = 0.1)
    # the reported objective is a FOCEi evaluation at the ITS estimates
    expect_identical(.fitIts$ofvType, "focei")
    expect_lt(abs(.fitIts$objf - .fitFocei$objf), 0.5)
    expect_true(!is.null(.fitIts$cov))
    # E-step moments are the MAP etas and the inverse inner Hessian
    .cv <- .fitIts$env$impCondVar
    expect_length(.cv, 12L)
    expect_true(all(vapply(.cv, function(m) all(eigen(m, only.values = TRUE)$values > 0), logical(1))))
  })

  test_that("an ITS warm start reaches the FOCEi optimum in fewer iterations", {
    .rough <- function() {
      ini({
        tka <- 0
        tcl <- 0
        tv <- 2.5
        eta.ka ~ 1
        eta.cl ~ 1
        eta.v ~ 1
        add.sd <- 2
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .ctl <- foceiControl(print = 0L, covMethod = "")
    .f0 <- suppressMessages(suppressWarnings(
      nlmixr2(.rough, nlmixr2data::theo_sd, "focei", .ctl)
    ))
    .fi <- suppressMessages(suppressWarnings(
      nlmixr2(.rough, nlmixr2data::theo_sd, "its", itsControl(print = 0L, maxIter = 10L, covMethod = ""))
    ))
    .f <- suppressMessages(suppressWarnings(nlmixr2(.fi, est = "focei", control = .ctl)))
    expect_lt(abs(.f$objf - .f0$objf), 0.05)
    expect_lt(max(.f$parHist$iter), max(.f0$parHist$iter))
  })

  test_that("interaction=FALSE reports a FOCE objective", {
    .f <- suppressMessages(suppressWarnings(
      nlmixr2(.oneCmt, nlmixr2data::theo_sd, "its", itsControl(print = 0L, interaction = FALSE, covMethod = ""))
    ))
    expect_identical(.f$ofvType, "foce")
  })

  test_that("the conditional-covariance correction removes the omega shrinkage", {
    .sim <- function() {
      ini({
        tcl <- log(4)
        tv <- log(40)
        eta.cl ~ 0.2
        eta.v ~ 0.2
        add.sd <- 0.1
      })
      model({
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        linCmt() ~ prop(add.sd)
      })
    }
    .ev <- rxode2::et(amt = 100) |>
      rxode2::et(c(2, 12)) |>
      rxode2::et(id = 1:80)
    .s <- rxode2::rxWithSeed(
      1460,
      suppressMessages(suppressWarnings(
        rxode2::rxSolve(.sim, .ev, addDosing = TRUE, returnType = "data.frame")
      ))
    )
    .dat <- .s[, c("id", "time", "amt", "evid", "sim")]
    names(.dat) <- c("ID", "TIME", "AMT", "EVID", "DV")
    .dat$DV[.dat$EVID != 0] <- 0
    .f <- suppressMessages(suppressWarnings(
      nlmixr2(.sim, .dat, "its", itsControl(print = 0L, covMethod = ""))
    ))
    .om <- diag(.f$omega)
    # the naive two-stage omega: the empirical second moment of the MAP etas
    .eta <- as.matrix(.f$eta[, c("eta.cl", "eta.v")])
    .naive <- colMeans(.eta^2)
    .cvMean <- Reduce(`+`, .f$env$impCondVar) / length(.f$env$impCondVar)
    # Omega = E[eta eta'] + mean conditional covariance at the fixed point
    expect_equal(unname(.om), unname(.naive + diag(.cvMean)), tolerance = 0.05)
    # sparse sampling shrinks the naive estimate; ITS recovers the truth
    expect_true(all(.naive < 0.8 * 0.2))
    expect_equal(unname(.om), c(0.2, 0.2), tolerance = 0.3)
  })

  test_that("a non-mu-referenced theta is estimated", {
    .nonMu <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 30
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- tv
        linCmt() ~ add(add.sd)
      })
    }
    .fi <- suppressMessages(suppressWarnings(
      nlmixr2(.nonMu, nlmixr2data::theo_sd, "its", itsControl(print = 0L, covMethod = ""))
    ))
    .ff <- suppressMessages(suppressWarnings(
      nlmixr2(.nonMu, nlmixr2data::theo_sd, "focei", foceiControl(print = 0L, covMethod = ""))
    ))
    expect_false(isTRUE(all.equal(.fi$theta[["tv"]], 30)))
    expect_equal(.fi$theta, .ff$theta, tolerance = 0.05)
    expect_lt(abs(.fi$objf - .ff$objf), 1)
  })
})
