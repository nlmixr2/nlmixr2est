.carryThetaMod <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    add.sd <- 0.5
  })
  model({
    cl <- exp(tcl) * (wt / 70)^0.75
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

# six subjects with a time-varying wt, observed from the model under nocb
# with additive noise
.carryThetaNlmDat <- function(ui) {
  d <- .carryFitDat(6L)
  set.seed(42)
  sim <- rxode2::rxSolve(ui, d, returnType = "data.frame", covsInterpolation = "nocb")
  d$dv <- ifelse(
    d$evid == 0,
    sim$cp[match(paste(d$id, d$time), paste(sim$id, sim$time))] +
      rnorm(nrow(d), 0, 0.3),
    0
  )
  d
}

# central difference in parameter i of the nlm objective that
# nlmObjectiveSetup() loaded, at the (scaled) parameters p
.carryThetaNlmCentral <- function(p, i, h = 1e-4) {
  pp <- p
  pm <- p
  pp[i] <- p[i] + h
  pm[i] <- p[i] - h
  (as.numeric(.nlmixrNlmFunC(pp)) - as.numeric(.nlmixrNlmFunC(pm))) / (2 * h)
}

# est = "nlm" (or "nls") with the carry on or off, under nocb; the tight
# tolerances are for an ODE twin (linCmt() is solved in closed form)
.carryThetaNlmFit <- function(ui, data, carry, est = "nlm") {
  .ctl <- if (est == "nlm") nlmixr2est::nlmControl else nlmixr2est::nlsControl
  suppressMessages(nlmixr2est::nlmixr2(
    ui,
    data,
    est = est,
    control = .ctl(
      print = 0,
      linCmtSensCarry = carry,
      rxControl = rxode2::rxControl(covsInterpolation = "nocb", atol = 1e-10, rtol = 1e-10)
    )
  ))
}

test_that("theta-side carry through a real est='nlm' fit (#1003)", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())

  ui0 <- suppressMessages(nlmixr2est::nlmixr2(.carryThetaMod))
  d <- .carryThetaNlmDat(ui0)
  fitC <- .carryThetaNlmFit(ui0, d, "auto")
  fitN <- .carryThetaNlmFit(ui0, d, "none")
  expect_true(is.finite(fitC$objf))
  # the population objective is a pure function of theta (the gradient does
  # not enter its value), so carry and naive converge to the same optimum;
  # what the carry fixes is the score itself (covariance, convergence path)
  expect_equal(fitC$objf, fitN$objf, tolerance = 1e-4)
  # and that optimum is the ODE twin's
  fitO <- .carryThetaNlmFit(rxode2::linToOde(ui0), d, "none")
  expect_equal(fitC$objf, fitO$objf, tolerance = 1e-6)
  expect_equal(unname(unlist(fitC$theta)), unname(unlist(fitO$theta)), tolerance = 1e-3)

  # value correctness: the linCmt surface equals the genuinely integrated
  # linToOde() surface at the converged point
  th <- unlist(fitC$theta)
  obs <- d[d$evid == 0, ]
  surf <- function(ui, useLin) {
    p <- c(tcl = th[["tcl"]], tv = th[["tv"]])
    s <- rxode2::rxSolve(
      ui,
      params = p,
      events = d,
      returnType = "data.frame",
      covsInterpolation = "nocb",
      useLinCmt = useLin,
      atol = 1e-10,
      rtol = 1e-10
    )
    f <- s$cp[match(paste(obs$id, obs$time), paste(s$id, s$time))]
    -2 * sum(dnorm(obs$dv, f, abs(th[["add.sd"]]), log = TRUE))
  }
  expect_equal(surf(ui0, TRUE), surf(rxode2::linToOde(ui0), FALSE), tolerance = 1e-6)

  # mechanism at the fit's own control: the built gradient model carries
  u <- rxode2::.copyUi(ui0)
  assign(
    "control",
    nlmixr2est::nlmControl(
      print = 0,
      linCmtSensCarry = "auto",
      rxControl = rxode2::rxControl(covsInterpolation = "nocb")
    ),
    envir = u
  )
  e <- suppressMessages(u$nlmEnv)
  expect_true(any(grepl("rx_lcConc_~linCmtB", strsplit(e$..nlmS, "\n")[[1]])))

  # score exactness at the initial estimates, where the score is far from 0
  # (at the optimum it is 0 and a relative error means nothing): the carried
  # score matches FD of the model's own rx_pred_ sum; the naive one is biased
  scoreErr <- function(carry) {
    u2 <- rxode2::.copyUi(ui0)
    assign(
      "control",
      nlmixr2est::nlmControl(
        print = 0,
        linCmtSensCarry = carry,
        rxControl = rxode2::rxControl(covsInterpolation = "nocb")
      ),
      envir = u2
    )
    e2 <- suppressMessages(u2$nlmEnv)
    m <- suppressWarnings(rxode2::rxode2(e2$..nlmS))
    pars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = 0.5)
    slv <- function(q) {
      rxode2::rxSolve(m, params = q, events = d, returnType = "data.frame", covsInterpolation = "nocb")
    }
    r0 <- slv(pars)
    h <- 1e-5
    a <- pars
    a[1] <- a[1] + h
    b <- pars
    b[1] <- b[1] - h
    fd <- sum((slv(a)$rx_pred_ - slv(b)$rx_pred_) / (2 * h))
    abs(sum(r0$rx__sens_rx_pred__BY_THETA_1___) - fd) / abs(fd)
  }
  expect_true(scoreErr("auto") < 1e-6)
  expect_true(scoreErr("none") > 1e-4)
})

test_that("est='nlm''s carried theta gradient is the derivative of its objective", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  # nlm's record walk after each subject's solve (nlmSolveGradId()) must
  # start the carry afresh: at the initial thetas the analytic gradient of
  # the objective equals its central differences
  ui0 <- suppressMessages(nlmixr2est::nlmixr2(.carryThetaMod))
  d <- .carryThetaNlmDat(ui0)
  ctl <- nlmixr2est::nlmControl(
    print = 0,
    linCmtSensCarry = "auto",
    rxControl = rxode2::rxControl(covsInterpolation = "nocb")
  )
  p <- suppressMessages(nlmObjectiveSetup(ui0, d, control = ctl, gradient = TRUE))
  on.exit(.nlmFreeEnv(), add = TRUE)
  f0 <- .nlmixrNlmFunC(p)
  fd <- vapply(seq_along(p), .carryThetaNlmCentral, numeric(1), p = p)
  expect_equal(attr(f0, "gradient"), fd, tolerance = 1e-6)
  # an evaluation does not depend on the ones before it
  expect_identical(as.numeric(.nlmixrNlmFunC(p)), as.numeric(f0))
})

test_that("theta-side carry through a real est='nls' fit reaches its ODE twin", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  # nls takes its Jacobian from the same record walk as nlm's gradient
  ui0 <- suppressMessages(nlmixr2est::nlmixr2(.carryThetaMod))
  d <- .carryThetaNlmDat(ui0)
  fitC <- .carryThetaNlmFit(ui0, d, "auto", est = "nls")
  fitO <- .carryThetaNlmFit(rxode2::linToOde(ui0), d, "none", est = "nls")
  expect_equal(fitC$objf, fitO$objf, tolerance = 1e-6)
  expect_equal(unname(unlist(fitC$theta)), unname(unlist(fitO$theta)), tolerance = 1e-3)
})
