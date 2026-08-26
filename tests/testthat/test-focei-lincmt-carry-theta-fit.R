test_that("theta-side carry through a real est='nlm' fit (#1003)", {
  skip_on_cran()
  skip_if_not(.rxFoceiLinCmtCarryCapable())

  modTheta <- function() {
    ini({
      tcl <- log(2)
      tv <- log(20)
      add.sd <- 0.5
    })
    model({
      cl <- exp(tcl) * (wt / 70)^0.75 # nolint: object_usage_linter.
      v <- exp(tv) # nolint: object_usage_linter.
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  d <- .carryFitDat(6L)
  set.seed(42)
  ui0 <- suppressMessages(nlmixr2est::nlmixr2(modTheta))
  sim <- rxode2::rxSolve(ui0, d, returnType = "data.frame", covsInterpolation = "nocb")
  d$dv <- ifelse(d$evid == 0,
    sim$cp[match(paste(d$id, d$time), paste(sim$id, sim$time))] +
      rnorm(nrow(d), 0, 0.3), 0
  )
  fitWith <- function(carry) {
    u <- rxode2::.copyUi(ui0)
    assign("control",
      nlmixr2est::nlmControl(
        print = 0, linCmtSensCarry = carry,
        rxControl = rxode2::rxControl(covsInterpolation = "nocb")
      ),
      envir = u
    )
    suppressMessages(nlmixr2est::nlmixr2(u, d, est = "nlm"))
  }
  fitC <- fitWith("auto")
  fitN <- fitWith("none")
  expect_true(is.finite(fitC$objf))
  # the population objective is a pure function of theta (the gradient does
  # not enter its value), so carry and naive converge to the same optimum;
  # what the carry fixes is the score itself (covariance, convergence path)
  expect_equal(fitC$objf, fitN$objf, tolerance = 1e-4)

  # value correctness: the linCmt surface equals the genuinely integrated
  # linToOde() surface at the converged point
  th <- unlist(fitC$theta)
  obs <- d[d$evid == 0, ]
  surf <- function(ui, useLin) {
    p <- c(tcl = th[["tcl"]], tv = th[["tv"]])
    s <- rxode2::rxSolve(ui,
      params = p, events = d, returnType = "data.frame",
      covsInterpolation = "nocb", useLinCmt = useLin,
      atol = 1e-10, rtol = 1e-10
    )
    f <- s$cp[match(paste(obs$id, obs$time), paste(s$id, s$time))]
    -2 * sum(dnorm(obs$dv, f, abs(th[["add.sd"]]), log = TRUE))
  }
  expect_equal(surf(ui0, TRUE), surf(rxode2::linToOde(ui0), FALSE), tolerance = 1e-6)

  # mechanism at the fit's own control: the built gradient model carries
  u <- rxode2::.copyUi(ui0)
  assign("control",
    nlmixr2est::nlmControl(
      print = 0, linCmtSensCarry = "auto",
      rxControl = rxode2::rxControl(covsInterpolation = "nocb")
    ),
    envir = u
  )
  e <- suppressMessages(u$nlmEnv)
  expect_true(any(grepl("rx_lcConc_~linCmtB", strsplit(e$..nlmS, "\n")[[1]])))

  # score exactness at the converged point: the carried score matches FD of
  # the model's own rx_pred_ sum; the naive one is biased
  scoreErr <- function(carry) {
    u2 <- rxode2::.copyUi(ui0)
    assign("control",
      nlmixr2est::nlmControl(
        print = 0, linCmtSensCarry = carry,
        rxControl = rxode2::rxControl(covsInterpolation = "nocb")
      ),
      envir = u2
    )
    e2 <- suppressMessages(u2$nlmEnv)
    m <- suppressWarnings(rxode2::rxode2(e2$..nlmS))
    pars <- c(`THETA[1]` = th[["tcl"]], `THETA[2]` = th[["tv"]], `THETA[3]` = th[["add.sd"]])
    slv <- function(q) {
      rxode2::rxSolve(m,
        params = q, events = d, returnType = "data.frame",
        covsInterpolation = "nocb"
      )
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
