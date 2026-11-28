nmTest({
  # an adjoint (s) method refuses the inner model (forward sensitivities without
  # the adjoint expansion), so a fit asking for one failed; it now solves with
  # the forward method of the same code (code - 200)
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  .fit <- function(method, est = "focei") {
    suppressMessages(nlmixr(
      one.cmt,
      nlmixr2data::theo_sd,
      est,
      foceiControl(
        maxOuterIterations = 0,
        print = 0,
        covMethod = "",
        calcTables = FALSE,
        rxControl = rxode2::rxControl(method = method)
      )
    ))
  }

  test_that("a focei fit with an adjoint method solves with its forward method", {
    for (.m in c("dop853s", "ros43s", "radauiia5s")) {
      .b <- sub("s$", "", .m)
      expect_message(
        nlmixr(
          one.cmt,
          nlmixr2data::theo_sd,
          "focei",
          foceiControl(
            maxOuterIterations = 0,
            print = 0,
            covMethod = "",
            calcTables = FALSE,
            rxControl = rxode2::rxControl(method = .m)
          )
        ),
        paste0("method='", .b, "'"),
        fixed = TRUE
      )
      expect_equal(.fit(.m)$objf, .fit(.b)$objf, label = .m)
    }
  })

  test_that("the mapping keeps forward methods and maps every adjoint code to its forward code", {
    .ctl <- rxode2::rxControl(method = "liblsoda")
    expect_identical(.rxControlAdjointForward(.ctl), .ctl)
    .idx <- rxode2::odeMethodToInt(NULL)
    for (.m in c("rk4s", "ros4s", "cvodesadj", "backwardEulers", "liblsodaadj")) {
      expect_identical(
        suppressMessages(.rxControlAdjointForward(rxode2::rxControl(method = .m)))$method,
        as.integer(.idx[[.m]] - 200L),
        label = .m
      )
    }
  })
})
