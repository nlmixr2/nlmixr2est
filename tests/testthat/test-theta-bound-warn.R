.thetaBoundIniDf <- function(est, lower, upper, fix = rep(FALSE, length(est))) {
  data.frame(
    ntheta = seq_along(est),
    name = paste0("t", seq_along(est)),
    est = est,
    lower = lower,
    upper = upper,
    fix = fix
  )
}

test_that("thetas at or near a finite bound are named", {
  .d <- .thetaBoundIniDf(
    est = c(0.2000005, 0.05, 1e-6, 3, 9.99999),
    lower = c(-1, -1, 0, -Inf, 0),
    upper = c(0.2, 0.2, Inf, Inf, 10)
  )
  expect_identical(.thetaAtBound(.d, 1e-4), c("t1", "t3", "t5"))
  # one-sided bound: scaled by the initial distance from the bound
  expect_identical(.thetaAtBound(.d, 1e-4, c(t3 = 1e-3)), c("t1", "t5"))
  expect_identical(.thetaAtBound(.d, NA_real_), character(0))
  # an interior estimate, infinite bounds and fixed thetas never count
  .d$fix <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
  expect_identical(.thetaAtBound(.d, 1e-4), character(0))
  # only thetas count
  .d$ntheta[2] <- NA
  expect_identical(.thetaAtBound(.d, 1), character(0))
  expect_identical(.thetaAtBound(NULL, 1e-4), character(0))
  # nor do the IOV rewrite's own bounded thetas
  .iov <- .thetaBoundIniDf(c(0, 0), c(0, 0), c(Inf, Inf))
  .iov$backTransform <- c("nlmixr2iovExp", NA)
  expect_identical(.thetaAtBound(.iov, 1e-4), "t2")
})

test_that("the tolerance is read from the fit's control", {
  expect_identical(.thetaBoundWarnTol(list(control = list())), 1e-4)
  expect_identical(.thetaBoundWarnTol(list(control = list(thetaBoundWarn = 0.01))), 0.01)
  expect_true(is.na(.thetaBoundWarnTol(list(control = list(thetaBoundWarn = NA)))))
  expect_true(is.na(.thetaBoundWarnTol(list(control = list(thetaBoundWarn = 0)))))
  expect_identical(foceiControl()$thetaBoundWarn, 1e-4)
  expect_true(is.na(foceiControl(thetaBoundWarn = NA)$thetaBoundWarn))
  expect_error(foceiControl(thetaBoundWarn = -1))
})

test_that("the warning names the thetas on one runInfo line", {
  .n <- 20
  .fit <- list(
    control = list(thetaBoundWarn = 1e-4),
    ui = list(iniDf = .thetaBoundIniDf(rep(0, .n), rep(0, .n), rep(1, .n)))
  )
  expect_warning(.at <- .thetaBoundWarn(.fit, NULL), "^estimate at its bound: \"t1\", ")
  expect_length(.at, .n)
  .msg <- tryCatch(.thetaBoundWarn(.fit, NULL), warning = function(w) conditionMessage(w))
  expect_lt(nchar(.msg), 75L)
  expect_match(.msg, "\\+[0-9]+ more$")
  .fit$control$thetaBoundWarn <- NA
  expect_no_warning(.thetaBoundWarn(.fit, NULL))
})

test_that("a theta pinned to its bound is reported in runInfo (#1505)", {
  skip_on_cran()
  .mod <- function() {
    ini({
      tka <- c(-1, 0.1, 0.2)
      tcl <- log(2.7)
      tv <- c(0, 3.45, 10)
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  .ctl <- foceiControl(print = 0, covMethod = "", calcTables = FALSE)
  .fit <- suppressMessages(nlmixr2(.mod, theo_sd, "focei", .ctl))
  expect_equal(.fit$theta[["tka"]], 0.2, tolerance = 1e-4)
  expect_true(any(.fit$runInfo == "estimate at its bound: \"tka\""))
  # an inactive bound does not warn
  .fit <- suppressMessages(nlmixr2(.mod |> ini(tka = c(-1, 0.45, 2)), theo_sd, "focei", .ctl))
  expect_false(any(grepl("^estimate at its bound", .fit$runInfo)))
})

test_that("every method's control carries thetaBoundWarn", {
  for (.ctl in c(
    "saemControl",
    "nlmeControl",
    "vaeControl",
    "emviControl",
    "nlminbControl",
    "nlmControl",
    "bobyqaControl",
    "newuoaControl",
    "uobyqaControl",
    "optimControl",
    "nlsControl",
    "trustControl"
  )) {
    .f <- get(.ctl, mode = "function")
    expect_identical(.f()$thetaBoundWarn, 1e-4, label = .ctl)
    expect_true(is.na(.f(thetaBoundWarn = NA)$thetaBoundWarn), label = .ctl)
    expect_error(.f(thetaBoundWarn = -1), "thetaBoundWarn", label = .ctl)
  }
})

test_that("a population-only fit notes a theta at its bound (#1505)", {
  skip_on_cran()
  .mod <- function() {
    ini({
      tka <- c(-1, 0.1, 0.2)
      tcl <- log(2.7)
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  .fit <- suppressMessages(nlmixr2(.mod, theo_sd, "nlminb", nlminbControl(print = 0, covMethod = "")))
  expect_equal(.fit$theta[["tka"]], 0.2, tolerance = 1e-4)
  expect_true(any(.fit$runInfo == "estimate at its bound: \"tka\""))
  .fit <- suppressMessages(
    nlmixr2(.mod, theo_sd, "nlminb", nlminbControl(print = 0, covMethod = "", thetaBoundWarn = NA))
  )
  expect_false(any(grepl("^estimate at its bound", .fit$runInfo)))
})
