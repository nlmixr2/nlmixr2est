# foceiControl(fast="auto"): the analytic outer gradient for closed-form models (#1223)
nmTest({
  .theoAlg <- function() {
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
      k <- cl / v
      cp <- WT * 4.02 * ka / (v * (ka - k)) * (exp(-k * TIME) - exp(-ka * TIME))
      cp ~ add(add.sd)
    })
  }
  .theoOde <- function() {
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
  .theoLin <- function() {
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
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  .dObs <- nlmixr2data::theo_sd[nlmixr2data::theo_sd$EVID == 0, ]
  .env <- function(data = .dObs) list2env(list(data = data))
  .resolve <- function(model, control = foceiControl(), type = "foceiControl", data = .dObs) {
    .ui <- rxode2::rxode2(model)
    suppressMessages(.foceiResolveFastAuto(control, .ui, .env(data), type))
  }

  test_that("fast defaults to \"auto\" and validates", {
    expect_identical(foceiControl()$fast, "auto")
    expect_identical(foceiControl()$outerOptTxt, "bobyqa")
    expect_error(foceiControl(fast = "yes"))
    expect_false(foceiControl(fast = FALSE)$fast)
  })

  test_that("a closed-form model resolves to the analytic gradient with nlminb", {
    .ctl <- .resolve(.theoAlg)
    expect_true(.ctl$fast)
    expect_identical(.ctl$outerOptTxt, "nlminb")
    expect_equal(.ctl$rel.tol, 10^(-.ctl$sigdig - 2))
    expect_true(.resolve(.theoAlg, type = "foceControl")$fast)
    expect_true(.resolve(.theoAlg, type = "mfoceiControl")$fast)
  })

  test_that("an explicit gradient-based outerOpt is kept", {
    .ctl <- .resolve(.theoAlg, foceiControl(outerOpt = "lbfgsb3c"))
    expect_true(.ctl$fast)
    expect_identical(.ctl$outerOptTxt, "lbfgsb3c")
  })

  test_that("everything else resolves to fast=FALSE", {
    for (.m in list(.theoOde, .theoLin)) {
      .ctl <- .resolve(.m)
      expect_false(.ctl$fast)
      expect_identical(.ctl$outerOptTxt, "bobyqa")
    }
    # derivative-free outerOpt, explicit or not
    expect_false(.resolve(.theoAlg, foceiControl(outerOpt = "bobyqa"))$fast)
    expect_false(.resolve(.theoAlg, foceiControl(outerOpt = "newuoa"))$fast)
    # methods without the auto choice, and the output pass ("")
    for (.t in c("laplaceControl", "agqControl", "foControl", "posthocControl", "impmapControl", "")) {
      expect_false(.resolve(.theoAlg, type = .t)$fast, info = .t)
    }
    # censored data is out of the analytic scope
    .dc <- .dObs
    .dc$CENS <- 0L
    .dc$CENS[1] <- 1L
    expect_false(.resolve(.theoAlg, data = .dc)$fast)
    # explicit fast=FALSE / fast=TRUE are left alone
    expect_false(.resolve(.theoAlg, foceiControl(fast = FALSE))$fast)
    .ctl <- .resolve(.theoOde, foceiControl(fast = TRUE))
    expect_true(.ctl$fast)
    expect_identical(.ctl$outerOptTxt, "trust")
  })

  test_that("a fast=\"auto\" closed-form focei fit matches the bobyqa optimum", {
    skip_on_cran()
    .ctl <- list(print = 0, covMethod = "", calcTables = FALSE)
    .fast <- suppressMessages(nlmixr2(.theoAlg, .dObs, "focei", control = .ctl))
    expect_true(.fast$control$fast)
    expect_identical(.fast$control$outerOptTxt, "nlminb")
    .slow <- suppressMessages(nlmixr2(.theoAlg, .dObs, "focei", control = c(.ctl, list(fast = FALSE))))
    expect_false(.slow$control$fast)
    expect_identical(.slow$control$outerOptTxt, "bobyqa")
    expect_lt(abs(.fast$objective - .slow$objective), 0.05)
  })
})
