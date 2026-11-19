test_that("collapsed and bound-pinned omegas are named", {
  .ini <- c(eta.cl = 0.1, eta.v = 300, eta.ka = 0.5)
  .final <- c(eta.cl = 0.08, eta.v = 0.0055, eta.ka = 0.005)
  expect_identical(
    .omegaCollapsed(.ini, .final, 1e-3),
    list(collapsed = "eta.v", bound = character(0))
  )
  expect_identical(
    .omegaCollapsed(.ini, .final, 1e-3, 100),
    list(collapsed = "eta.v", bound = "eta.ka")
  )
  # a legitimately small omega near its start is left alone
  expect_identical(
    .omegaCollapsed(c(eta.cl = 0.01), c(eta.cl = 0.008), 1e-3, 100),
    list(collapsed = character(0), bound = character(0))
  )
  expect_identical(
    .omegaCollapsed(.ini, .final, NA_real_, 100),
    list(collapsed = character(0), bound = character(0))
  )
  # etas missing from either side are ignored
  expect_identical(
    .omegaCollapsed(c(eta.q = 1), c(eta.cl = 1e-9), 1e-3, 100)$collapsed,
    character(0)
  )
})

test_that("only free, positive diagonal omegas have an initial value", {
  .iniDf <- data.frame(
    name = c("tcl", "eta.cl", "eta.v", "(eta.cl,eta.v)", "eta.ka", "eta.q"),
    neta1 = c(NA, 1, 2, 2, 3, 4),
    neta2 = c(NA, 1, 2, 1, 3, 4),
    est = c(1, 0.1, 0.2, 0.01, 0.3, 0),
    fix = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)
  )
  expect_identical(.omegaCollapseIni(.iniDf), c(eta.cl = 0.1, eta.v = 0.2))
  expect_identical(.omegaCollapseIni(NULL), setNames(numeric(0), character(0)))
})

test_that("the threshold and bound are read from the fit's control", {
  expect_identical(.omegaCollapseWarnLimit(list(control = list())), 1e-3)
  expect_identical(.omegaCollapseWarnLimit(list(control = list(omegaCollapseWarn = 0.01))), 0.01)
  expect_true(is.na(.omegaCollapseWarnLimit(list(control = list(omegaCollapseWarn = NA)))))
  expect_true(is.na(.omegaCollapseWarnLimit(list(control = list(omegaCollapseWarn = 0)))))
  expect_identical(.omegaCollapseBoundLower(list(control = list(diagOmegaBoundLower = 100))), 100)
  expect_true(is.na(.omegaCollapseBoundLower(list(control = list(diagOmegaBoundLower = 1)))))
  expect_true(is.na(.omegaCollapseBoundLower(list(control = list()))))
  expect_identical(foceiControl()$omegaCollapseWarn, 1e-3)
  expect_true(is.na(foceiControl(omegaCollapseWarn = NA)$omegaCollapseWarn))
  expect_error(foceiControl(omegaCollapseWarn = -1))
  expect_error(foceiControl(omegaCollapseWarn = 2))
})

test_that("the warning names the offending etas on one runInfo line", {
  .env <- new.env(parent = emptyenv())
  .env$ui <- list(
    iniDf = data.frame(
      name = c("eta.cl", "eta.v"),
      neta1 = c(1, 2),
      neta2 = c(1, 2),
      est = c(300, 300),
      fix = c(FALSE, FALSE)
    )
  )
  .fit <- list(
    control = list(omegaCollapseWarn = 1e-3, diagOmegaBoundLower = 100),
    omega = matrix(c(0.005, 0, 0, 3), 2, 2, dimnames = list(c("eta.cl", "eta.v"), c("eta.cl", "eta.v")))
  )
  expect_warning(
    expect_warning(.omegaCollapseWarn(.fit, .env), "omega collapsed \\(<0.001 x initial\\): eta.cl"),
    "omega at lower bound \\(initial/100\\): eta.v"
  )
  .fit$control$omegaCollapseWarn <- NA
  expect_no_warning(.omegaCollapseWarn(.fit, .env))
  .msg <- .omegaCollapseWarnMsg(paste0("eta.parameter.", 1:20), 1e-3)
  expect_lt(nchar(.msg), 75L)
  expect_match(.msg, "\\+[0-9]+ more$")
})

test_that("a random-slope omega that collapses is reported in runInfo (#1394)", {
  skip_on_cran()
  set.seed(1394)
  .d <- do.call(
    rbind,
    lapply(seq_len(30), function(i) {
      .b0 <- 250 + rnorm(1, 0, sqrt(600))
      .b1 <- 10 + rnorm(1, 0, sqrt(35))
      data.frame(ID = i, TIME = 0:9, DV = .b0 + .b1 * (0:9) + rnorm(10, 0, 25))
    })
  )
  .mod <- function() {
    ini({
      tint <- 200
      tsl <- 5
      eta.int ~ 400
      eta.sl ~ 300
      add.sd <- 100
    })
    model({
      b0 <- tint + eta.int
      b1 <- tsl + eta.sl
      y <- b0 + b1 * time
      y ~ add(add.sd)
    })
  }
  .ctl <- function(...) foceiControl(print = 0, covMethod = "", calcTables = FALSE, ...)
  # without the default omega bound the slope variance collapses toward 0
  .fit <- suppressMessages(nlmixr2(.mod, .d, "focei", .ctl(diagOmegaBoundLower = 1e8)))
  expect_lt(.fit$omega["eta.sl", "eta.sl"], 0.3)
  expect_true(any(grepl("omega collapsed (<0.001 x initial): eta.sl", .fit$runInfo, fixed = TRUE)))
  # with the default bound, a far too large start stops on the bound
  .modBig <- .mod |> ini(eta.sl ~ 30000)
  .fit <- suppressMessages(nlmixr2(.modBig, .d, "focei", .ctl()))
  expect_equal(.fit$omega["eta.sl", "eta.sl"], 300, tolerance = 0.01)
  expect_true(any(grepl("omega at lower bound (initial/100): eta.sl", .fit$runInfo, fixed = TRUE)))
  # an interior estimate does not warn
  .fit <- suppressMessages(nlmixr2(.mod, .d, "focei", .ctl()))
  expect_gt(.fit$omega["eta.sl", "eta.sl"], 10)
  expect_false(any(grepl("^omega (collapsed|at lower bound)", .fit$runInfo)))
})

test_that("a block omega is checked on the parameters FOCEi bounds", {
  .nm <- c("eta.a", "eta.b")
  .om0 <- matrix(c(1, 0.5, 0.5, 1), 2, 2, dimnames = list(.nm, .nm))
  # eta.a's bounded (inverse-Cholesky) parameter on the bound, while its
  # marginal variance stays well above 1/100 of the start
  .r <- chol(solve(.om0))
  .r[1, ] <- .r[1, ] * 10
  .om <- solve(crossprod(.r))
  dimnames(.om) <- list(.nm, .nm)
  expect_gt(.om["eta.a", "eta.a"], 0.2)
  expect_identical(.omegaCollapseAtBound(.om0, .om, 100), "eta.a")
  expect_identical(
    .omegaCollapsed(diag(.om0), diag(.om), 1e-3, 100, .om0, .om),
    list(collapsed = character(0), bound = "eta.a")
  )
  expect_identical(.omegaCollapsed(diag(.om0), diag(.om), 1e-3, 100)$bound, character(0))
  # for a diagonal omega the bounded parameter and the variance ratio agree
  .d0 <- diag(c(eta.a = 300, eta.b = 2))
  dimnames(.d0) <- list(.nm, .nm)
  .d1 <- diag(c(3, 1.9))
  dimnames(.d1) <- list(.nm, .nm)
  expect_identical(.omegaCollapseAtBound(.d0, .d1, 100), "eta.a")
  expect_identical(
    .omegaCollapsed(diag(.d0), diag(.d1), 1e-3, 100)$bound,
    .omegaCollapsed(diag(.d0), diag(.d1), 1e-3, 100, .d0, .d1)$bound
  )
  # mismatched or singular matrices fall back to the ratio rule
  expect_null(.omegaCollapseAtBound(.d0, .d1[2:1, 2:1], 100))
  expect_null(.omegaCollapseAtBound(.d0, matrix(0, 2, 2, dimnames = list(.nm, .nm)), 100))
})
