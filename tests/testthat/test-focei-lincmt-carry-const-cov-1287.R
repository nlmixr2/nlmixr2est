nmTest({
  # #1287: past the carry columns, a pair whose covariates are constant within
  # every subject (here v's lwt) has an exact standard gradient, so it is
  # dropped before the build stops.
  .d <- nlmixr2data::theo_sd
  .d$lwt <- log(.d$WT / 70)
  .d$SW <- as.numeric(.d$TIME >= 6)

  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      b.cl <- -0.2
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + 0.75 * lwt + b.cl * SW + eta.cl)
      v <- exp(tv + lwt + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  test_that("subject-constant pairs are dropped first past the cap, without a warning", {
    .pairs <- data.frame(
      slot = c(1L, 2L),
      etaName = c("eta.cl", "eta.v"),
      covs = c("lwt,SW", "lwt"),
      shape = c("mult", "mult"),
      fD = NA_character_,
      lagD = NA_character_,
      stringsAsFactors = FALSE
    )
    local_mocked_bindings(.rxFoceiLinCmtCarryMaxPairs = function() 3L)
    expect_silent(.kept <- .rxFoceiLinCmtCarryCapPairs(.pairs, c("lwt", "WT")))
    expect_identical(.kept$etaName, "eta.cl")
    # without the data's constant covariates nothing is dropped
    expect_identical(.rxFoceiLinCmtCarryCapPairs(.pairs, NULL)$etaName, c("eta.cl", "eta.v"))
  })

  test_that("a model past the cap fits as with every pair carried", {
    skip_if_not(.rxFoceiLinCmtCarryCapable())
    .ctl <- foceiControl(maxOuterIterations = 0L, covMethod = "", print = 0, calcTables = FALSE)
    .fitAll <- suppressMessages(nlmixr2(.mod, .d, est = "focei", control = .ctl))
    # cl/eta.cl and v/eta.v need 4 carry columns; with 3, v/eta.v is dropped
    local_mocked_bindings(.rxFoceiLinCmtCarryMaxPairs = function() 3L)
    .fit3 <- suppressMessages(nlmixr2(.mod, .d, est = "focei", control = .ctl))
    expect_equal(.fit3$objf, .fitAll$objf, tolerance = 1e-8)
    expect_false(any(grepl("carry columns", .fit3$runInfo, fixed = TRUE)))
  })
})
