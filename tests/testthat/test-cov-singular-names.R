# A failed FOCEi covariance step names the parameters at fault (#1398)
nmTest({
  skip_on_cran()

  test_that("a theta the data cannot inform is named when the covariance step fails", {
    # tz never enters the prediction, so R and S get an exact zero row
    .one <- function() {
      ini({
        tcl <- log(2)
        tz <- log(0.5)
        tv <- log(30)
        eta.v ~ 0.1
        add.sd <- 0.5
      })
      model({
        cl <- exp(tcl) + (t > 1e9) * exp(tz)
        v <- exp(tv + eta.v)
        ke <- cl / v
        d / dt(center) <- -ke * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .fit <- suppressMessages(suppressWarnings(nlmixr2(
      .one,
      nlmixr2data::theo_sd,
      "focei",
      control = foceiControl(print = 0, covFull = FALSE, calcTables = FALSE)
    )))
    expect_identical(.fit$covMethod, "failed")
    expect_true("covariance step failed: tz uninformed" %in% .fit$runInfo)
  })

  test_that("too few subjects for the S matrix is said, and a healthy fit is unchanged", {
    .two <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .d <- nlmixr2data::theo_sd
    .fit <- suppressMessages(suppressWarnings(nlmixr2(
      .two,
      .d[.d$ID %in% c(1, 2), ],
      "focei",
      control = foceiControl(print = 0, covMethod = "s", covFull = FALSE, calcTables = FALSE)
    )))
    expect_identical(.fit$covMethod, "failed")
    expect_true("covariance step failed: 2 subjects < 4 parameters" %in% .fit$runInfo)
    .fit <- suppressMessages(suppressWarnings(nlmixr2(
      .two,
      .d,
      "focei",
      control = foceiControl(print = 0, covFull = FALSE, calcTables = FALSE)
    )))
    expect_identical(.fit$covMethod, "r,s")
    expect_false(any(grepl("not PD|covariance step failed", .fit$runInfo)))
  })
})
