nmTest({
  # non-positive observed DV on a log/Box-Cox endpoint warns before fitting (#1510)
  .modFun <- function(err) {
    eval(bquote(function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        add.sd <- 0.7
        lam <- 0.5
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        .(err)
      })
    }))
  }
  .lnorm <- rxode2::rxode2(.modFun(quote(cp ~ lnorm(add.sd))))
  .boxCox <- rxode2::rxode2(.modFun(quote(cp ~ add(add.sd) + boxCox(lam))))
  .boxCox1 <- rxode2::rxode2(.modFun(quote(cp ~ add(add.sd) + boxCox(1))))
  .yj <- rxode2::rxode2(.modFun(quote(cp ~ add(add.sd) + yeoJohnson(lam))))

  .dat <- nlmixr2data::theo_sd
  .dat <- .dat[.dat$ID <= 2, ]
  .dat$DV[.dat$EVID == 0 & .dat$DV <= 0] <- 0.01
  .obs <- which(.dat$EVID == 0)
  .dat0 <- .dat
  .dat0$DV[.obs[2:4]] <- c(0, -1, NA)

  .run <- function(ui, data, est = "focei") {
    .preProcessDvDomainWarn(ui, est, data, NULL)
  }

  test_that("lnorm/boxCox endpoints warn with the count of DV <= 0", {
    expect_warning(.run(.lnorm, .dat0), "^2 DV <= 0 on lnorm endpoint 'cp'$")
    expect_warning(.run(.boxCox, .dat0), "^2 DV <= 0 on boxCox endpoint 'cp'$")
    expect_lt(nchar("100000 DV <= 0 on boxCox endpoint 'a_long_endpoint_name'"), 75)
  })

  test_that("all-positive DV, simulation and exempt transforms are silent", {
    expect_no_warning(.run(.lnorm, .dat))
    expect_no_warning(.run(.lnorm, .dat0, "rxSolve"))
    expect_no_warning(.run(.yj, .dat0))
    expect_no_warning(.run(.boxCox1, .dat0))
  })

  test_that("censored rows and dose rows are not counted", {
    .d <- .dat0
    .d$CENS <- 0L
    .d$CENS[.obs[2:3]] <- 1L
    expect_no_warning(.run(.lnorm, .d))
    .d <- .dat
    .d$DV[.d$EVID != 0] <- 0
    expect_no_warning(.run(.lnorm, .d))
  })

  test_that("multiple endpoints are matched by CMT or DVID", {
    .two <- rxode2::rxode2(function() {
      ini({
        tv <- 1
        a1 <- 0.7
        a2 <- 0.7
      })
      model({
        e1 <- tv
        e1 ~ add(a1)
        e2 <- tv
        e2 ~ lnorm(a2)
      })
    })
    .d <- data.frame(ID = 1, TIME = 1:4, DV = c(-1, 0, 2, -3), CMT = c("e1", "e1", "e2", "e2"))
    expect_warning(.run(.two, .d), "^1 DV <= 0 on lnorm endpoint 'e2'$")
    .d$CMT <- NULL
    .d$DVID <- c(1, 1, 2, 2)
    expect_warning(.run(.two, .d), "^1 DV <= 0 on lnorm endpoint 'e2'$")
    .d$DVID <- NULL
    expect_no_warning(.run(.two, .d))
  })

  test_that("the warning reaches the fit's runInfo", {
    skip_on_cran()
    fit <- suppressMessages(suppressWarnings(
      nlmixr2(.lnorm, .dat0, est = "focei", control = list(print = 0))
    ))
    expect_true(any(grepl("2 DV <= 0 on lnorm endpoint 'cp'", fit$runInfo, fixed = TRUE)))
  })
})
