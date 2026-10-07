# A zero Shi (2021) epsilon degenerated the finite-difference step search (#1174).

test_that("shiErr/hessErr must be strictly positive", {
  for (.e in list(0, -1e-6)) {
    expect_error(nlmControl(shiErr = .e), "'shiErr' must be > 0")
    expect_error(nlmControl(hessErr = .e), "'hessErr' must be > 0")
    expect_error(nlminbControl(shiErr = .e), "'shiErr' must be > 0")
    expect_error(nlminbControl(hessErr = .e), "'hessErr' must be > 0")
    expect_error(nlsControl(shiErr = .e), "'shiErr' must be > 0")
    expect_error(optimControl(shiErr = .e), "'shiErr' must be > 0")
    expect_error(trustControl(hessErr = .e), "'hessErr' must be > 0")
  }
  expect_error(nlmControl(shiErr = NA_real_))
  expect_error(nlmControl(shiErr = c(1e-6, 1e-6)))
  expect_equal(nlmControl(shiErr = 1e-6, hessErr = 2e-6)[c("shiErr", "hessErr")],
               list(shiErr = 1e-6, hessErr = 2e-6))
})

test_that("a hand-built control with shiErr = 0 falls back to the default", {
  skip_on_cran()
  mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      tlag <- log(0.2)
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      alag(depot) <- exp(tlag)
      d / dt(central) <- ka * depot - cl / v * central
      cp <- central / v
      cp ~ add(add.sd)
    })
  }
  .grad <- function(ctl) {
    x <- nlmObjectiveSetup(mod, nlmixr2data::theo_sd, control = ctl,
                           gradient = TRUE, scale = "natural")
    on.exit(.nlmFreeEnv())
    nlmLikEvalC_(x)$grad
  }
  .ctl <- nlmControl(print = 0L, eventSens = "fd")
  .ref <- .grad(.ctl)
  expect_true(all(.ref != 0))
  .ctl$shiErr <- 0
  .ctl$hessErr <- 0
  expect_equal(.grad(.ctl), .ref)
})
