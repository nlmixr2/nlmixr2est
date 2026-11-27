.aboveRegex <- "final objective .* is above the initial estimates"

test_that(".nlmStartObjective takes the earliest finite objective", {
  .ph <- data.frame(
    iter = c(1L, 1L, 2L, 3L),
    type = c("Scaled", "Unscaled", "Scaled", "Scaled"),
    objf = c(NA, NA, 120, 100)
  )
  expect_equal(.nlmStartObjective(.ph), 120)
  .ph$objf[1] <- 150
  expect_equal(.nlmStartObjective(.ph), 150)
  expect_true(is.na(.nlmStartObjective(NULL)))
  expect_true(is.na(.nlmStartObjective(data.frame(iter = 1L, objf = NA_real_))))
})

test_that(".nlmAboveStartWarn applies the band on the -2LL scale", {
  # 1% of the doubled start (200) is 2, i.e. 1 in -LL units
  expect_warning(expect_true(.nlmAboveStartWarn(100, 101.1, 0.01)), .aboveRegex)
  expect_false(expect_silent(.nlmAboveStartWarn(100, 100.9, 0.01)))
  expect_false(expect_silent(.nlmAboveStartWarn(100, 90, 0.01)))
  # the default band when the control carries no aboveStartTol
  expect_warning(.nlmAboveStartWarn(100, 102, NULL), .aboveRegex)
  # NA/Inf disable it; an unknown start never warns
  expect_false(expect_silent(.nlmAboveStartWarn(100, 1e10, Inf)))
  expect_false(expect_silent(.nlmAboveStartWarn(100, 1e10, NA_real_)))
  expect_false(expect_silent(.nlmAboveStartWarn(NA_real_, 1e10, 0.01)))
  .w <- tryCatch(.nlmAboveStartWarn(100, 200, 0.01), warning = conditionMessage)
  expect_lt(nchar(.w), 75L)
})

test_that("nlm-family controls validate aboveStartTol", {
  for (.ctl in c(
    "nlmControl",
    "nlminbControl",
    "bobyqaControl",
    "n1qn1Control",
    "newuoaControl",
    "optimControl",
    "trustControl",
    "uobyqaControl",
    "lbfgsb3cControl"
  )) {
    .f <- get(.ctl)
    expect_equal(.f()$aboveStartTol, 0.01, info = .ctl)
    expect_equal(.f(aboveStartTol = Inf)$aboveStartTol, Inf, info = .ctl)
    expect_error(.f(aboveStartTol = -1), info = .ctl)
  }
})

nmTest({
  .one <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
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
  .runInfo <- function(fit) as.character(unlist(fit$runInfo))

  test_that("a converging nlm-family fit does not note its start", {
    .fit <- suppressMessages(nlmixr(
      .one,
      nlmixr2data::theo_sd,
      est = "nlminb",
      control = nlminbControl(print = 0, covMethod = "")
    ))
    # the recorded start is the objective at the initial estimates, above the final
    .start <- .nlmStartObjective(.fit$env$parHistData)
    expect_gt(2 * .start, .fit$objf)
    expect_false(any(grepl(.aboveRegex, .runInfo(.fit))))
  })

  test_that("a start below the final objective is noted in $runInfo", {
    # a tampered history: the start looks far better than where the fit ended
    local_mocked_bindings(.nlmStartObjective = function(parHistData) 1)
    .fit <- suppressMessages(nlmixr(
      .one,
      nlmixr2data::theo_sd,
      est = "nlminb",
      control = nlminbControl(print = 0, covMethod = "")
    ))
    expect_equal(sum(grepl(.aboveRegex, .runInfo(.fit))), 1L)
    .fit <- suppressMessages(nlmixr(
      .one,
      nlmixr2data::theo_sd,
      est = "nlminb",
      control = nlminbControl(print = 0, covMethod = "", aboveStartTol = Inf)
    ))
    expect_false(any(grepl(.aboveRegex, .runInfo(.fit))))
  })
})
