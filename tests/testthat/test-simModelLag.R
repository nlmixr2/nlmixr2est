test_that(".simModelLaggedVars finds history-function arguments (#1173)", {
  .x <- quote({
    a <- b + lag(c0) + diff(d, 2)
    e <- first(f) - last(g) + lead(h)
    k <- lag(1 + m)
  })
  expect_setequal(.simModelLaggedVars(.x), c("c0", "d", "f", "g", "h"))
  expect_identical(.simModelLaggedVars(quote(a <- b)), character(0))
})

test_that("vpcSim() and augPred() work with lag() of a calculated variable (#1173)", {
  skip_on_cran()
  mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.f ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(central) <- ka * depot - cl / v * central
      c0 <- central / v
      cp <- (0.5 * c0 + 0.5 * lag(c0)) * exp(eta.f)
      cp ~ add(add.sd)
    })
  }
  .ui <- rxode2::rxode2(mod)
  .sim <- .getSimModel(.ui, hideIpred = FALSE)
  .txt <- deparse(.sim)
  # the lagged variable stays a real lhs; other calculated variables do not
  expect_true(any(grepl("c0 <- central/v", .txt, fixed = TRUE)))
  expect_true(any(grepl("ka ~ exp(tka)", .txt, fixed = TRUE)))
  expect_identical(attr(.sim, "lagged"), "c0")
  # a lagged parameter is not an lhs turned back on, so is not dropped
  .ui2 <- rxode2::model(.ui, cp <- (0.5 * c0 + 0.5 * lag(tka)) * exp(eta.f))
  expect_identical(attr(.getSimModel(.ui2), "lagged"), character(0))
  expect_true(any(grepl("c0 ~ central/v", deparse(.getSimModel(.ui2)), fixed = TRUE)))
  # rxode2 allows only lag()/diff() of a calculated variable
  for (.fn in c("lag", "diff")) {
    .cp <- str2lang(sprintf("cp <- (0.5 * c0 + 0.5 * %s(c0)) * exp(eta.f)", .fn))
    .ui3 <- eval(bquote(rxode2::model(.ui, .(.cp))))
    expect_error(eval(.getSimModel(.ui3)), NA)
  }

  fit <- nlmixr2(
    mod,
    nlmixr2data::theo_sd,
    est = "focei",
    control = foceiControl(print = 0L, maxOuterIterations = 0L, covMethod = "")
  )
  v <- vpcSim(fit, n = 2)
  expect_s3_class(v, "nlmixr2vpcSim")
  expect_equal(length(unique(v$sim.id)), 2L)
  expect_false(any(is.na(v$sim)))
  # the lhs kept for lag() is not an extra output column
  expect_false("c0" %in% names(v))
  vp <- vpcSim(fit, n = 2, pred = TRUE)
  expect_false(any(is.na(vp$pred)))

  ap <- augPred(fit)
  expect_true(all(c("Individual", "Population") %in% levels(ap$ind)))
  expect_false(any(is.na(ap$values)))
})
