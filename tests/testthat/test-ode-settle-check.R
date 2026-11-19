test_that("the note is one runInfo line and silent when settled", {
  expect_null(.odeSettleMsg(0.05, 0.1))
  expect_null(.odeSettleMsg(-0.1, 0.1))
  expect_identical(.odeSettleMsg(0.284, 0.1), "OFV moves 0.284 at 10x tighter ODE tolerance")
  expect_identical(.odeSettleMsg(-1.5, 0.1), "OFV moves -1.5 at 10x tighter ODE tolerance")
  expect_identical(.odeSettleMsg(NA_real_, 0.1), "OFV not finite at 10x tighter ODE tolerance")
  expect_identical(.odeSettleMsg(Inf, 0.1), "OFV not finite at 10x tighter ODE tolerance")
  expect_lt(nchar(.odeSettleMsg(-123456.789, 0.1)), 75)
})

test_that("every ODE tolerance is tightened 10x", {
  .rx <- rxode2::rxControl(rtol = 1e-4, atol = 1e-6, rtolSens = 1e-5, atolSens = 1e-7)
  .t <- .odeSettleTighten(.rx)
  for (.n in c("atol", "rtol", "atolSens", "rtolSens", "ssAtol", "ssRtol", "ssAtolSens", "ssRtolSens")) {
    expect_equal(.t[[.n]], .rx[[.n]] / 10)
  }
  expect_identical(.t$method, .rx$method)
})

test_that("foceiControl() validates odeSettleCheck/odeSettleTol", {
  .c <- foceiControl()
  expect_false(.c$odeSettleCheck)
  expect_identical(.c$odeSettleTol, 0.1)
  expect_true(foceiControl(odeSettleCheck = TRUE, odeSettleTol = 0.5)$odeSettleCheck)
  expect_error(foceiControl(odeSettleCheck = NA))
  expect_error(foceiControl(odeSettleTol = -1))
  expect_true(laplaceControl(odeSettleCheck = TRUE)$odeSettleCheck)
})

test_that("no re-solve when the check is off or no outer step ran", {
  local_mocked_bindings(.odeSettleObjf = function(fit) stop("re-solved"))
  expect_null(.odeSettleCheck(list(control = list(odeSettleCheck = FALSE, maxOuterIterations = 100L))))
  expect_null(.odeSettleCheck(list(control = list(maxOuterIterations = 100L))))
  expect_null(.odeSettleCheck(list(control = list(odeSettleCheck = TRUE, maxOuterIterations = 0L))))
  expect_null(.odeSettleCheck(list(control = NULL)))
})

nmTest({
  .oneCmt <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
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

  test_that("a one-compartment ODE fit at default tolerances is settled", {
    skip_on_cran()
    .f <- suppressMessages(nlmixr2(
      .oneCmt,
      nlmixr2data::theo_sd,
      "focei",
      foceiControl(print = 0, covMethod = "", calcTables = FALSE, odeSettleCheck = TRUE)
    ))
    expect_true(is.finite(.f$env$odeSettleObjf))
    expect_lt(abs(.f$env$odeSettleObjf - .f$objective), 0.1)
    expect_false(any(grepl("ODE tolerance", .f$runInfo)))
  })

  test_that("with the check off no extra solve occurs", {
    skip_on_cran()
    local_mocked_bindings(.odeSettleObjf = function(fit) stop("re-solved"))
    .f <- suppressMessages(nlmixr2(
      .oneCmt,
      nlmixr2data::theo_sd,
      "focei",
      foceiControl(print = 0, covMethod = "", calcTables = FALSE)
    ))
    expect_false(exists("odeSettleObjf", envir = .f$env, inherits = FALSE))
  })

  test_that("a linCmt() model is not re-solved", {
    skip_on_cran()
    .lin <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
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
    local_mocked_bindings(.odeSettleObjf = function(fit) stop("re-solved"))
    .f <- suppressMessages(nlmixr2(
      .lin,
      nlmixr2data::theo_sd,
      "focei",
      foceiControl(print = 0, covMethod = "", calcTables = FALSE, odeSettleCheck = TRUE)
    ))
    expect_false(exists("odeSettleObjf", envir = .f$env, inherits = FALSE))
  })

  test_that("a feedback ODE model at loose tolerance is flagged", {
    skip_on_cran()
    .fb <- function() {
      ini({
        tkin <- log(10)
        tkout <- log(0.5)
        tic50 <- log(2)
        tgam <- log(0.8)
        tcl <- log(1)
        tv <- log(10)
        eta.kin ~ 0.1
        eta.cl ~ 0.1
        add.sd <- 0.3
      })
      model({
        kin <- exp(tkin + eta.kin)
        kout <- exp(tkout)
        ic50 <- exp(tic50)
        gam <- exp(tgam)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        r0 <- kin / kout
        r(0) <- r0
        d / dt(center) <- -cl / v * center
        cp <- center / v
        d / dt(r) <- kin * (r0 / r)^gam * (1 - cp / (ic50 + cp)) - kout * r
        r ~ add(add.sd)
      })
    }
    .ev <- rxode2::et(amt = 100, cmt = "center") |>
      rxode2::et(c(0.5, 1, 2, 4, 6, 8, 12, 16, 24, 36, 48)) |>
      rxode2::et(id = 1:20)
    set.seed(1393)
    .eta <- data.frame(id = 1:20, eta.kin = rnorm(20, sd = 0.3), eta.cl = rnorm(20, sd = 0.3))
    .sim <- as.data.frame(rxode2::rxSolve(.fb, .ev, params = .eta, omega = NA, sigma = NA, addDosing = TRUE))
    .d <- data.frame(id = .sim$id, time = .sim$time, amt = .sim$amt, evid = .sim$evid, cmt = "center", dv = 0)
    .obs <- .d$evid == 0
    .d$cmt[.obs] <- "r"
    .d$dv[.obs] <- .sim$r[.obs] + rnorm(sum(.obs), sd = 0.3)
    .d$amt[.obs] <- 0
    .tol <- rxode2::rxControl(rtol = 1e-2, atol = 1e-2, rtolSens = 1e-2, atolSens = 1e-2)
    .f <- suppressMessages(nlmixr2(
      .fb,
      .d,
      "focei",
      foceiControl(print = 0, covMethod = "", calcTables = FALSE, odeSettleCheck = TRUE, rxControl = .tol)
    ))
    .shift <- .f$env$odeSettleObjf - .f$objective
    expect_gt(abs(.shift), 0.1)
    .msg <- grep("ODE tolerance", .f$runInfo, value = TRUE)
    expect_identical(.msg, .odeSettleMsg(.shift, 0.1))
  })
})
