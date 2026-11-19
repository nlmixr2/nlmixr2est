.shrinkTab <- function(sd) {
  .rn <- c(
    "mean",
    "var",
    "sd",
    "skewness",
    "kurtosis",
    "var shrinkage (%)",
    "sd shrinkage (%)",
    "t statistic",
    "p-value"
  )
  .d <- as.data.frame(lapply(sd, function(v) c(rep(0, 6), v, 0, 1)))
  names(.d) <- names(sd)
  row.names(.d) <- .rn
  .d
}

test_that("only etas at or above the threshold are named", {
  .s <- .shrinkTab(c(eta.cl = 45, eta.v = 30, eta.ka = 29.9, IWRES = 80))
  expect_identical(.shrinkEtaHigh(.s, 30), c("eta.cl", "eta.v"))
  expect_identical(.shrinkEtaHigh(.s, 50), character(0))
})

test_that("a disabled threshold names nothing", {
  .s <- .shrinkTab(c(eta.cl = 95))
  expect_identical(.shrinkEtaHigh(.s, NA_real_), character(0))
  expect_identical(.shrinkEtaHigh(.s, Inf), character(0))
})

test_that("non-finite shrinkage is ignored", {
  .s <- .shrinkTab(c(eta.cl = NaN, eta.v = NA, eta.ka = Inf, eta.q = 40))
  expect_identical(.shrinkEtaHigh(.s, 30), "eta.q")
})

test_that("IOV shrinkage lists are searched per table", {
  .s <- list(
    id = .shrinkTab(c(eta.cl = 10, eta.v = 50)),
    occ = .shrinkTab(c(rx.iov.cl.1 = 80, rx.iov.cl.2 = 85, iov.cl = 60, IWRES = 90))
  )
  expect_identical(.shrinkEtaHigh(.s, 30), c("eta.v", "iov.cl"))
  # the collapsed saem IOV fit keeps only per-occasion etas
  expect_identical(.shrinkEtaHigh(.shrinkTab(c(rx.eta.cl.1 = 80, eta.ka = 40)), 30), "eta.ka")
})

test_that("the threshold is read from the fit's control", {
  expect_identical(.shrinkEtaWarnLimit(list(control = list())), 30)
  expect_identical(.shrinkEtaWarnLimit(list(control = list(shrinkEtaWarn = 20))), 20)
  expect_true(is.na(.shrinkEtaWarnLimit(list(control = list(shrinkEtaWarn = NA)))))
  expect_true(is.na(.shrinkEtaWarnLimit(list(control = list(shrinkEtaWarn = Inf)))))
  expect_identical(foceiControl()$shrinkEtaWarn, 30)
  expect_true(is.na(foceiControl(shrinkEtaWarn = NA)$shrinkEtaWarn))
  expect_error(foceiControl(shrinkEtaWarn = -1))
})

test_that("the warning names the offending etas on one runInfo line", {
  .fit <- list(
    control = list(shrinkEtaWarn = 30),
    shrink = .shrinkTab(c(eta.cl = 45, eta.v = 10, IWRES = 80))
  )
  expect_warning(.shrinkEtaWarn(.fit), "^high eta shrinkage \\(>=30%\\): eta.cl$")
  .fit$control$shrinkEtaWarn <- NA
  expect_no_warning(.shrinkEtaWarn(.fit))
  .fit$control$shrinkEtaWarn <- 50
  expect_no_warning(.shrinkEtaWarn(.fit))
  .fit$shrink <- NULL
  .fit$control$shrinkEtaWarn <- 30
  expect_no_warning(.shrinkEtaWarn(.fit))
})

test_that("a long eta list is clipped under 75 characters", {
  .etas <- sprintf("eta.parameter%02d", 1:10)
  .msg <- .shrinkEtaWarnMsg(.etas, 30)
  expect_lt(nchar(.msg), 75)
  expect_match(.msg, "\\+[0-9]+ more$")
  expect_identical(.shrinkEtaWarnMsg("eta.cl", 30), "high eta shrinkage (>=30%): eta.cl")
})

nmTest({
  test_that("a fit reports high eta shrinkage in $runInfo", {
    skip_on_cran()
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- log(c(0, 2.7, 100))
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
    .fit <- function(warn) {
      suppressMessages(nlmixr2(
        one.cmt,
        nlmixr2data::theo_sd,
        "focei",
        foceiControl(print = 0, maxOuterIterations = 0, covMethod = "", shrinkEtaWarn = warn)
      ))
    }
    .f <- .fit(0)
    .sh <- unlist(.f$shrink["sd shrinkage (%)", c("eta.ka", "eta.cl", "eta.v")])
    .msg <- grep("^high eta shrinkage", .f$runInfo, value = TRUE)
    expect_length(.msg, 1L)
    expect_identical(.msg, .shrinkEtaWarnMsg(names(.sh)[.sh >= 0], 0))
    .off <- .fit(NA)
    expect_false(any(grepl("^high eta shrinkage", .off$runInfo)))
  })
})
