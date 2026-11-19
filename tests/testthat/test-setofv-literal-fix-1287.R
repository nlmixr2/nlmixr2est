nmTest({
  # #1287: setOfv() kept a fixed covariate term as a parameter, so a linCmt()
  # model whose time-varying covariate only enters through a fixed theta could
  # exceed the sensitivity-carry columns and stop.
  test_that("setOfv() substitutes fixed thetas before the FOCEi evaluation", {
    skip_if_not(.rxFoceiLinCmtCarryCapable())
    .d <- nlmixr2data::theo_sd
    .d$SW <- as.numeric(.d$TIME >= 6)

    .base <- function() {
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
    # the same model: the flag effect on cl is fixed at 0
    .flag0 <- .base |>
      model(cl <- exp(tcl + b.cl * SW + eta.cl)) |>
      ini(b.cl = fix(0), append = "tv")

    .ctl <- saemControl(nBurn = 10, nEm = 10, seed = 42, print = 0, covMethod = "", calcTables = FALSE)
    .fitBase <- suppressMessages(nlmixr2(.base, .d, est = "saem", control = .ctl))
    .fitFlag0 <- suppressMessages(nlmixr2(.flag0, .d, est = "saem", control = .ctl))

    # with b.cl kept as a parameter, cl/eta.cl would need 2 carry columns
    local_mocked_bindings(.rxFoceiLinCmtCarryMaxPairs = function() 1L)
    suppressMessages(setOfv(.fitBase, "focei"))
    suppressMessages(setOfv(.fitFlag0, "focei"))
    expect_true(is.finite(.fitFlag0$objf))
    expect_equal(.fitFlag0$objf, .fitBase$objf)
    expect_equal(.fitFlag0$theta[["b.cl"]], 0)
  })
})
