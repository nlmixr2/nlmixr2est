nmTest({
  .simWeighted <- function(model, nSub = 40L, seed = 20261125L) {
    set.seed(seed)
    .ev <- rxode2::et(amt = 320, cmt = "depot") |>
      rxode2::et(c(0.25, 0.5, 1, 2, 3.5, 5, 7, 9, 12, 24)) |>
      rxode2::et(id = seq_len(nSub))
    .cov <- data.frame(id = seq_len(nSub), WT = round(runif(nSub, 50, 110)))
    .ev <- merge(as.data.frame(.ev), .cov, by = "id")
    .s <- rxode2::rxSolve(model, .ev, addDosing = TRUE, returnType = "data.frame")
    .s <- merge(.s[, setdiff(names(.s), "WT")], .cov, by = "id")
    .d <- data.frame(
      ID = .s$id,
      TIME = .s$time,
      EVID = .s$evid,
      AMT = ifelse(is.na(.s$amt), 0, .s$amt),
      CMT = ifelse(.s$evid == 1, 1L, 2L),
      DV = ifelse(.s$evid == 1, 0, .s$sim),
      WT = .s$WT
    )
    .d[order(.d$ID, .d$TIME, -.d$EVID), ]
  }
  .base <- function(errLine) {
    eval(parse(
      text = sprintf(
        "function() {
      ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(30); add.sd <- 2; eta.ka ~ 0.2; eta.cl ~ 0.1 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- 10 * linCmt(); %s })
    }",
        errLine
      )
    ))
  }
  .ctl <- saemControl(seed = 42L, print = 0L, nBurn = 100L, nEm = 100L)

  test_that("saem fits a time-weighted additive residual in closed form (#1463)", {
    skip_on_cran()
    mTrue <- .base("cp ~ add(add.sd * exp(-time / 8))")
    .d <- .simWeighted(mTrue)
    .f <- .nlmixr(rxode2::ini(mTrue, add.sd = 1), .d, est = "saem", control = .ctl)
    .ref <- .nlmixr(rxode2::ini(mTrue, add.sd = 1), .d, est = "focei", control = foceiControl(print = 0L))
    # no general likelihood and no temporary eta
    expect_false(any(grepl("dnorm() likelihood", .f$runInfo, fixed = TRUE)))
    expect_false(any(grepl("temporary eta", .f$runInfo, fixed = TRUE)))
    expect_false(any(grepl("^rx[.]", c(names(fixef(.f)), rownames(.f$omega)))))
    expect_equal(fixef(.f)[["add.sd"]], fixef(.ref)[["add.sd"]], tolerance = 0.1)
    expect_equal(fixef(.f)[["add.sd"]], 2, tolerance = 0.15)
    # the reported model, tables and objective use the residual as written
    expect_true(any(grepl("exp(-time/8)", vapply(.f$finalUi$lstExpr, deparse1, character(1)), fixed = TRUE)))
    expect_equal(.f$IWRES, .f$IRES / (fixef(.f)[["add.sd"]] * exp(-.f$TIME / 8)), tolerance = 1e-6)
    expect_true(sd(.f$IWRES) < 2)
    .focei <- setOfv(.f, "focei")$objective
    .eval <- .nlmixr(
      .f$finalUi,
      .d,
      est = "focei",
      control = foceiControl(print = 0L, maxOuterIterations = 0L, covMethod = "")
    )$objective
    expect_equal(.focei, .eval, tolerance = 1e-3)
  })

  test_that("a covariate weight on one of two endpoints (#1463)", {
    skip_on_cran()
    mTrue <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(30); add.sd <- 2; add.sd2 <- 1
        eta.ka ~ 0.2; eta.cl ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        cp <- 10 * linCmt()
        eff <- cp / 2
        cp ~ add(add.sd * WT / 70)
        eff ~ add(add.sd2)
      })
    }
    .ui <- rxode2::rxode2(mTrue)
    expect_length(.saemResidWeightSpecs(.ui, .saemModeledResidualCond(.ui)), 1L)
    set.seed(11)
    .ev <- rxode2::et(amt = 320, cmt = "depot") |>
      rxode2::et(c(0.5, 1, 2, 4, 8, 12, 24), cmt = "cp") |>
      rxode2::et(c(1, 4, 12), cmt = "eff") |>
      rxode2::et(id = 1:30)
    .ev <- merge(as.data.frame(.ev), data.frame(id = 1:30, WT = round(runif(30, 50, 110))), by = "id")
    .s <- rxode2::rxSolve(mTrue, .ev, addDosing = TRUE, returnType = "data.frame")
    .d <- data.frame(
      ID = .s$id,
      TIME = .s$time,
      EVID = .s$evid,
      AMT = ifelse(is.na(.s$amt), 0, .s$amt),
      CMT = .s$CMT,
      DV = ifelse(.s$evid == 1, 0, .s$sim),
      WT = .s$WT
    )
    .d <- .d[order(.d$ID, .d$TIME, -.d$EVID), ]
    .f <- .nlmixr(mTrue, .d, est = "saem", control = .ctl)
    expect_false(any(grepl("dnorm() likelihood", .f$runInfo, fixed = TRUE)))
    expect_equal(fixef(.f)[["add.sd"]], 2, tolerance = 0.2)
    expect_equal(fixef(.f)[["add.sd2"]], 1, tolerance = 0.2)
  })
})
