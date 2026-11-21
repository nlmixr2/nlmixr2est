.atolDvOde <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    prop.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl)
    v <- exp(tv)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ prop(prop.sd)
  })
}

.atolDvLinCmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    prop.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl)
    v <- exp(tv)
    linCmt() ~ prop(prop.sd)
  })
}

.atolDvAlg <- function() {
  ini({
    a <- 1
    add.sd <- 0.7
  })
  model({
    y <- a * time
    y ~ add(add.sd)
  })
}

.atolDvTiny <- function() {
  d <- nlmixr2data::theo_sd
  d$DV <- d$DV * 1e-6
  d$AMT <- d$AMT * 1e-6
  d
}

test_that("observations near atol warn only for ODE models", {
  .ode <- rxode2::rxode2(.atolDvOde)
  .tiny <- .atolDvTiny()
  .ctl <- foceiControl()
  expect_warning(
    .preProcessAtolDv(.ode, "focei", .tiny, .ctl),
    "DV near ODE atol"
  )
  # order-1 data, a tight atol, or a disabled check do not warn
  expect_no_warning(.preProcessAtolDv(.ode, "focei", nlmixr2data::theo_sd, .ctl))
  .tight <- foceiControl(rxControl = rxode2::rxControl(atol = 1e-14))
  expect_no_warning(.preProcessAtolDv(.ode, "focei", .tiny, .tight))
  withr::with_options(
    list(nlmixr2est.atolCheck = FALSE),
    expect_no_warning(.preProcessAtolDv(.ode, "focei", .tiny, .ctl))
  )
  # simulation does not warn
  expect_no_warning(.preProcessAtolDv(.ode, "rxSolve", .tiny, .ctl))
  # linCmt() and algebraic models never warn
  expect_no_warning(.preProcessAtolDv(rxode2::rxode2(.atolDvLinCmt), "focei", .tiny, .ctl))
  expect_no_warning(.preProcessAtolDv(rxode2::rxode2(.atolDvAlg), "focei", .tiny, .ctl))
})

test_that("only uncensored observations are checked", {
  .d <- data.frame(
    ID = 1,
    TIME = c(0, 1, 2, 3),
    EVID = c(1, 0, 0, 0),
    AMT = c(1, 0, 0, 0),
    DV = c(1e-9, 0, 1e-9, 2),
    CENS = c(0, 0, 1, 0)
  )
  expect_equal(.atolDvMinObs(.d), 2)
  .d$CENS <- 0
  expect_equal(.atolDvMinObs(.d), 1e-9)
  .d$EVID <- NULL
  expect_equal(.atolDvMinObs(.d), 1e-9)
  expect_null(.atolDvMinObs(.d[1, ]))
  expect_null(.atolDvControlAtol(list()))
  expect_equal(.atolDvControlAtol(rxode2::rxControl(atol = 1e-6)), 1e-6)
})

test_that("the warning reaches $runInfo", {
  skip_on_cran()
  .fit <- suppressMessages(nlmixr2(
    .atolDvOde,
    .atolDvTiny(),
    "posthoc",
    control = foceiControl(print = 0, calcTables = FALSE)
  ))
  expect_true(any(grepl("DV near ODE atol", .fit$runInfo)))
})

test_that("atolCheck = FALSE in a control disables the check", {
  .ode <- rxode2::rxode2(.atolDvOde)
  .tiny <- .atolDvTiny()
  expect_no_warning(.preProcessAtolDv(.ode, "focei", .tiny, foceiControl(atolCheck = FALSE)))
  expect_no_warning(.preProcessAtolDv(.ode, "saem", .tiny, saemControl(atolCheck = FALSE)))
  expect_warning(.preProcessAtolDv(.ode, "saem", .tiny, saemControl()), "DV near ODE atol")
  # the focei family forwards it
  expect_false(laplaceControl(atolCheck = FALSE)$atolCheck)
  # a non-default value survives deparsing
  expect_equal(
    deparse(rxUiDeparse(foceiControl(atolCheck = FALSE), "ctl")),
    "ctl <- foceiControl(atolCheck = FALSE)"
  )
  expect_equal(
    deparse(rxUiDeparse(saemControl(atolCheck = FALSE), "ctl")),
    "ctl <- saemControl(atolCheck = FALSE)"
  )
})
