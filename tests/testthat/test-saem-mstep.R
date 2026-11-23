test_that("saem M-step table classifies each theta", {
  one <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  ui <- rxode2::rxode2(one)
  pars <- ui$saemParamsToEstimate
  cfg <- list(
    distribution = 1L,
    nMix = 1L,
    nonMuTheta = "regress",
    nonMuThetaOptType = 2L,
    i0 = match("tv", pars) - 1L
  )
  tbl <- .saemMstepTable(ui, cfg)
  expect_equal(
    setNames(tbl$update, tbl$theta)[c("tka", "tcl", "tv")],
    c(tka = "closed-form", tcl = "closed-form", tv = "regress")
  )
  expect_equal(
    .saemMstepNote(tbl, cfg),
    "closed form: tka, tcl; numerical (newuoa): tv"
  )
  cfg$nonMuTheta <- "eta"
  tbl <- .saemMstepTable(ui, cfg)
  expect_equal(tbl$update[tbl$theta == "tv"], "stochastic")
  expect_equal(.saemMstepNote(tbl, cfg), "closed form: tka, tcl; stochastic: tv")
  # a mixture never takes the regress refinement
  cfg$nonMuTheta <- "regress"
  cfg$nMix <- 2L
  expect_equal(.saemMstepTable(ui, cfg)$update[3], "stochastic")
  # a fixed theta is held
  ui2 <- rxode2::rxode2(one)
  ui2 <- rxode2::ini(ui2, tv = fix(3.45))
  cfg$nMix <- 1L
  tbl <- .saemMstepTable(ui2, cfg)
  expect_equal(tbl$update[tbl$theta == "tv"], "fixed")
  expect_null(.saemMstepNote(tbl, cfg))
})

test_that("saem M-step note truncates long theta lists", {
  tbl <- data.frame(
    theta = c(paste0("tcf", 1:10), paste0("treg", 1:10)),
    update = rep(c("closed-form", "regress"), each = 10)
  )
  note <- .saemMstepNote(tbl, list(distribution = 1L, nonMuThetaOptType = 2L))
  expect_lte(nchar(note), 74L)
  expect_match(note, "closed form: tcf1, tcf2, ", fixed = TRUE)
  expect_match(note, "numerical (newuoa): treg1, ", fixed = TRUE)
  expect_match(note, "more$")
  tbl <- data.frame(theta = c("tka", "tv"), update = c("penalized", "fixed"))
  expect_equal(.saemMstepNote(tbl, list()), "closed form+penalty: tka")
})

test_that("saem fit reports its M-step split", {
  skip_on_cran()
  one <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  ctl <- saemControl(print = 0, nBurn = 10, nEm = 10, calcTables = FALSE, covMethod = "")
  fit <- suppressWarnings(nlmixr2(one, theo_sd, "saem", ctl))
  expect_true("closed form: tka, tcl; numerical (newuoa): tv" %in% fit$runInfo)
  expect_equal(
    setNames(fit$saemMstep$update, fit$saemMstep$theta)[c("tka", "tcl", "tv")],
    c(tka = "closed-form", tcl = "closed-form", tv = "regress")
  )
  ctl$nonMuTheta <- "eta"
  fit <- suppressWarnings(nlmixr2(one, theo_sd, "saem", ctl))
  expect_true("closed form: tka, tcl; stochastic: tv" %in% fit$runInfo)

  mu <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
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
  fit <- suppressWarnings(nlmixr2(mu, theo_sd, "saem", ctl))
  expect_true(all(fit$saemMstep$update == "closed-form"))
  expect_false(any(grepl("closed form:", fit$runInfo, fixed = TRUE)))
})
