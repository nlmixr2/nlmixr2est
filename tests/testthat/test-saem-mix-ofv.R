nmTest({
  test_that(".saemGqMix integrates each component over its own etas only", {
    .cfg <- list(nphi1 = 4L, omegaShareSubpop = c(0L, 1L, 2L, 0L))
    .r <- .saemGqMix(list(mixProb = matrix(c(0.3, 0.7), ncol = 1)), .cfg)
    expect_equal(.r$prob, c(0.3, 0.7))
    expect_equal(.r$active, list(c(1L, 2L, 4L), c(1L, 3L, 4L)))
    # a fit stored without omegaShareSubpop integrates every eta per component
    .r <- .saemGqMix(list(mixProb = c(0.3, 0.7)), list(nphi1 = 4L))
    expect_equal(.r$active, list(1:4, 1:4))
    # no mixture: one pass over every eta, no component weights
    .r <- .saemGqMix(list(mixProb = numeric(0)), .cfg)
    expect_null(.r$prob)
    expect_equal(.r$active, list(1:4))
  })

  .mixCtl <- saemControl(print = 0, seed = 1, nBurn = 50, nEm = 50, covMethod = 0L)

  test_that("a shared-eta saem mixture's quadrature -2LL matches importance sampling (#1184)", {
    mod <- function() {
      ini({
        tka <- 0.45
        tcl1 <- log(1.5)
        tcl2 <- log(5)
        tv <- 3.45
        p1 <- 0.5
        add.sd <- 0.7
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    fit <- suppressMessages(.nlmixr(mod, theo_sd, est = "saem", control = .mixCtl))
    # the default objective is an uncalculated placeholder that setOfv() fills in
    expect_true(is.na(fit$objDf$OBJF[1]))
    suppressMessages(setOfv(fit, "foce"))
    expect_true(is.finite(fit$objDf["FOCE", "OBJF"]))

    .gq1 <- suppressMessages(calc.2LL(fit$saem, nnodes.gq = 3, nsd.gq = 1.6, fit$phiM))
    .gq2 <- suppressMessages(calc.2LL(fit$saem, nnodes.gq = 3, nsd.gq = 1.6, fit$phiM))
    # each component is solved under its own mixest, not a stale one
    expect_identical(.gq1, .gq2)

    suppressMessages(setOfv(fit, "gauss15_5"))
    suppressMessages(setOfv(fit, "imp"))
    expect_lt(abs(fit$objDf["gauss15_5", "OBJF"] - fit$objDf["IMP", "OBJF"]), 2)
    expect_equal(fit$ofvType, "IMP")
  })

  test_that("a split-eta saem mixture reports a finite FOCEi and quadrature objective (#1184)", {
    mod <- function() {
      ini({
        tka <- 0.45
        tcl1 <- log(1.5)
        tcl2 <- log(5)
        tv <- 3.45
        p1 <- 0.5
        add.sd <- 0.7
        eta.ka ~ 0.6
        eta.cl1 ~ 0.3
        eta.cl2 ~ 0.3
        eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl1), p1, exp(tcl2 + eta.cl2))
        v <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    fit <- suppressMessages(.nlmixr(mod, theo_sd, est = "saem", control = .mixCtl))
    expect_false(any(is.na(fit$etaMat)))
    expect_equal(colnames(fit$etaMat), colnames(fit$omega))
    # the compressed cfg keeps which component owns each eta
    expect_equal(attr(fit$saem, "saem.cfg")$omegaShareSubpop, c(0L, 1L, 2L, 0L))
    suppressMessages(setOfv(fit, "focei"))
    expect_true(is.finite(fit$objDf["FOCEi", "OBJF"]))
    suppressMessages(setOfv(fit, "gauss3_1.6"))
    expect_true(is.finite(fit$objDf["gauss3_1.6", "OBJF"]))
  })
})
