# The linCmt() sensitivity carry for IOV models: the rewrite puts eta.cl
# and every occasion eta in cl's slot, the occasion etas with the multCov
# shape d(cl)/d(eta) = cl * |omega| * (occ == j); the pairs share the slot's
# tracker.  Shared fixtures in helper-lincmt-carry.R.

.carryIovLin <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    iov.cl ~ 0.04 | occ
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.carryIovOde <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    iov.cl ~ 0.04 | occ
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv + eta.v)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

# occasions and a time-varying covariate in the same slot
.carryIovWtLin <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    bwt <- 0.75
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    iov.cl ~ 0.04 | occ
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + bwt * log(wt / 70) + eta.cl + iov.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}

.carryIovWtOde <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    bwt <- 0.75
    add.sd <- 0.7
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    iov.cl ~ 0.04 | occ
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + bwt * log(wt / 70) + eta.cl + iov.cl)
    v <- exp(tv + eta.v)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

# theo_md with the occasion stepping at 144 h, or every 24 h (8 occasions),
# and wt stepping at 72 h
.carryIovDat <- function(every24 = FALSE) {
  .d <- nlmixr2data::theo_md
  .d$occ <- if (every24) 1L + pmin(floor(.d$TIME / 24), 7L) else 1L + (.d$TIME >= 144)
  .d$wt <- .d$WT + 10 * (.d$TIME >= 72)
  .d
}

# focei at the initial thetas, with the etas estimated or held at 0
.carryIovFit <- function(mod, data, iovMethod = "theta", maxInner = 1000L, carry = "auto", interp = "locf") {
  suppressWarnings(suppressMessages(nlmixr2est::nlmixr2(
    mod,
    data,
    "focei",
    control = nlmixr2est::foceiControl(
      print = 0L,
      maxOuterIterations = 0L,
      maxInnerIterations = maxInner,
      covMethod = "",
      calcTables = FALSE,
      iovMethod = iovMethod,
      linCmtSensCarry = carry,
      rxControl = rxode2::rxControl(atol = 1e-10, rtol = 1e-10, covsInterpolation = interp)
    )
  )))
}

# an occasion model written out in the IOV rewrite's form, with wt in the
# same slot: d(cl)/d(eta.o1) = cl * iovsd * (occ == 1)
.carryOccHand <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    iovsd <- 0.2
    add.sd <- 0.5
    eta.cl ~ 0.1
    eta.o1 ~ 0.1
    eta.o2 ~ 0.1
  })
  model({
    cl <- exp(tcl + 0.75 * log(wt / 70) + eta.cl + iovsd * (eta.o1 * (occ == 1) + eta.o2 * (occ == 2)))
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

.carryIovInnerTxt <- function(fit) {
  paste(rxode2::rxNorm(fit$env$innerModel), collapse = "\n")
}

test_that("occasion etas sharing a slot with eta.cl and wt match FD", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .pars <- c(
    `THETA[1]` = log(2),
    `THETA[2]` = log(20),
    `THETA[3]` = 0.2,
    `THETA[4]` = 0.5,
    `ETA[1]` = 0.3,
    `ETA[2]` = 0.4,
    `ETA[3]` = -0.5
  )
  for (.kind in c("infusion", "addl", "every")) {
    .ev <- .carryMidEv(.kind)
    .ev$occ <- 1L + (.ev$time >= 24)
    for (.i in c("locf", "nocb", "midpoint")) {
      .r <- suppressWarnings(.carryJumpFd(.carryOccHand, .pars, .ev, interp = .i))
      expect_true(grepl("rx_lcCarryS2r", .r$txt, fixed = TRUE), label = paste(.kind, .i))
      expect_lt(max(.r$err), 1e-6, label = paste(.kind, .i))
    }
  }
})

test_that("a carried linCmt() IOV fit equals its ODE twin in both IOV expansions", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .d <- .carryIovDat()
  for (.m in c("theta", "omega")) {
    # etas held at 0 (the occasion etas' g still steps at 144 h) and estimated
    for (.mi in c(0L, 1000L)) {
      .lin <- .carryIovFit(.carryIovLin, .d, .m, .mi)
      .ode <- .carryIovFit(.carryIovOde, .d, .m, .mi)
      .lab <- paste(.m, .mi)
      expect_true(grepl("rx_lcCarryS2r", .carryIovInnerTxt(.lin), fixed = TRUE), label = .lab)
      expect_equal(.lin$objf, .ode$objf, tolerance = 1e-6, label = .lab)
      expect_equal(.lin$etaObf$OBJI, .ode$etaObf$OBJI, tolerance = 1e-6, label = .lab)
      expect_equal(.lin$eta$rx.iov.cl.2, .ode$eta$rx.iov.cl.2, tolerance = 1e-5, label = .lab)
    }
  }
})

test_that("occasions and a time-varying covariate in one slot fit as the ODE twin", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  .d <- .carryIovDat()
  for (.i in c("locf", "nocb")) {
    .lin <- .carryIovFit(.carryIovWtLin, .d, interp = .i)
    .ode <- .carryIovFit(.carryIovWtOde, .d, interp = .i)
    expect_equal(.lin$objf, .ode$objf, tolerance = 1e-6, label = .i)
    expect_equal(.lin$eta$eta.cl, .ode$eta$eta.cl, tolerance = 1e-5, label = .i)
  }
})

test_that("occasion etas past the carry column limit keep the standard gradient", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  # eta.cl and 8 occasion etas: 9 pairs and cl's tracker need 10 columns
  .d <- .carryIovDat(every24 = TRUE)
  .msg <- testthat::capture_messages(
    .fit <- suppressWarnings(nlmixr2est::nlmixr2(
      .carryIovLin,
      .d,
      "focei",
      control = nlmixr2est::foceiControl(
        print = 0L,
        maxOuterIterations = 0L,
        covMethod = "",
        calcTables = FALSE,
        rxControl = rxode2::rxControl(atol = 1e-10, rtol = 1e-10)
      )
    ))
  )
  expect_true(is.finite(.fit$objf))
  expect_true(any(grepl(
    "rx.iov.cl.7, rx.iov.cl.8 keep the standard gradient (carry column limit)",
    .msg,
    fixed = TRUE
  )))
  # seven pairs are carried: eta.cl and six occasion etas
  expect_true(grepl("rx_lcCarryS6r", .carryIovInnerTxt(.fit), fixed = TRUE))
  expect_false(grepl("rx_lcCarryS7r", .carryIovInnerTxt(.fit), fixed = TRUE))
})
