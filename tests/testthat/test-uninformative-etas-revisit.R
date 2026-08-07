nmTest({
  test_that("saemControl exposes revisitUninformativeEtas", {

    expect_true(saemControl()$revisitUninformativeEtas)
    expect_false(saemControl(revisitUninformativeEtas = FALSE)$revisitUninformativeEtas)
    expect_error(saemControl(revisitUninformativeEtas = "yes"))
    expect_error(saemControl(revisitUninformativeEtas = c(TRUE, TRUE)))
    expect_error(saemControl(revisitUninformativeEtas = NA))
  })
})

nmTest({
  test_that("the uninformative-eta verdict is re-decided at the end of burn-in", {

    ## Subjects 11, 12, 13 dose into the central compartment (IV), so `ka` -- and
    ## therefore eta.ka -- cannot be identified from their data at ANY theta.
    ## Subjects 21, 22, 23 dose into the depot (PO), where it can.
    .dat <- data.frame(
      ID = c(rep(11, 7), rep(12, 7), rep(13, 8), rep(21, 7), rep(22, 7), rep(23, 7)),
      TIME = c(0, 0.05, 0.25, 0.5, 1, 3, 5,   0, 0.05, 0.25, 0.5, 1, 3, 5,
               0, 0.05, 0.25, 0.5, 1, 3, 5, 8, 0, 0.25, 0.5, 1, 3, 5, 8,
               0, 0.25, 0.5, 1, 3, 5, 8,      0, 0.25, 0.5, 1, 3, 5, 8),
      DV = c(NA, 2017.85, 1323.74, 792.5, 822.72, 36.27, 3.33,
             NA, 1702, 1290.75, 1095.95, 907.6, 125.44, 14.44,
             NA, 1933.04, 1242.43, 661.22, 193.52, 1.75, NA, NA,
             NA, 706.58, 1063.14, 2257.62, 941.33, 629.69, 100,
             NA, 1462.95, 2217.76, 2739.5, 705.3, 108.47, 8.75,
             NA, 211.66, 467.23, 174.24, 153.6, 27.07, 2.81),
      AMT = c(1, rep(0, 6), 1, rep(0, 6), 1, rep(0, 7),
              5, rep(0, 6), 5, rep(0, 6), 5, rep(0, 6)),
      EVID = c(1, rep(0, 6), 1, rep(0, 6), 1, rep(0, 7),
               1, rep(0, 6), 1, rep(0, 6), 1, rep(0, 6)),
      CMT = c(rep(2, 7), rep(2, 7), rep(2, 8),
              1, rep(2, 6), 1, rep(2, 6), 1, rep(2, 6)))

    ## `tv` is the only difference between the two models below.
    .mod <- function(tvIni) {
      .f <- function() {
        ini({
          tka <- 0.45
          tcl <- -7
          tv <- -8
          eta.ka ~ 0.6
          eta.cl ~ 0.3
          eta.v ~ 0.1
          prop.sd <- 0.7
        })
        model({
          ka <- exp(tka + eta.ka)
          cl <- exp(tcl + eta.cl)
          v <- exp(tv + eta.v)
          d/dt(depot) <- -ka * depot
          d/dt(center) <- ka * depot - cl/v * center
          cp <- center/v
          cp ~ prop(prop.sd)
        })
      }
      .b <- body(.f)
      .b[[2]][[2]][[4]][[3]] <- tvIni
      body(.f) <- .b
      .f
    }

    .fit <- function(tvIni, revisit) {
      .nlmixr(.mod(tvIni)(), .dat, "saem",
              control = saemControl(nBurn = 60, nEm = 40, print = 0, seed = 42,
                                    calcTables = FALSE, logLik = FALSE, covMethod = "",
                                    revisitUninformativeEtas = revisit))
    }

    ## `tv = 20` puts v at 4.9e8, so every prediction underflows at the initial
    ## estimates.  The first test then (correctly) refuses to call anything
    ## uninformative -- there is no prediction to move -- and without a revisit that
    ## verdict stands for the whole fit, leaving eta.ka free even for the IV subjects.
    .off <- .fit(20, FALSE)
    expect_equal(unname(.off$saem$ueRevisitInfo[["ran"]]), 0L)
    expect_false(all(.off$etaObf$eta.ka[1:3] == 0))

    ## With the revisit, burn-in moves tv somewhere sensible, the predictions become
    ## real, and the test is re-run against them.
    .on <- .fit(20, TRUE)
    expect_equal(unname(.on$saem$ueRevisitInfo[["ran"]]), 1L)
    ## the mechanism must be shown to have DONE something -- an unchanged mask is
    ## indistinguishable from a revisit that never happened
    expect_gt(.on$saem$ueRevisitInfo[["froze"]], 0L)
    expect_true(all(.on$etaObf$eta.ka[1:3] == 0))

    ## From good initial estimates the first test already gets the IV subjects right,
    ## and the revisit must agree with it rather than undo it.
    .good <- .fit(-8, TRUE)
    expect_equal(unname(.good$saem$ueRevisitInfo[["ran"]]), 1L)
    expect_true(all(.good$etaObf$eta.ka[1:3] == 0))
    ## eta.ka stays estimable for the PO subjects
    expect_false(all(.good$etaObf$eta.ka[4:6] == 0))
  })
})

nmTest({
  test_that("the revisit is skipped when the mask itself is off", {

    .d <- nlmixr2data::theo_sd
    .m <- function() {
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
    .ctl <- function(...) {
      saemControl(nBurn = 20, nEm = 20, print = 0, seed = 42, calcTables = FALSE,
                  logLik = FALSE, covMethod = "", ...)
    }

    ## handleUninformativeEtas=FALSE means there is no mask to re-decide; the revisit
    ## must not run and reintroduce one.
    .noMask <- .nlmixr(.m(), .d, "saem",
                       control = .ctl(handleUninformativeEtas = FALSE))
    expect_equal(unname(.noMask$saem$ueRevisitInfo[["ran"]]), 0L)

    ## Explicitly disabled.
    .noRev <- .nlmixr(.m(), .d, "saem",
                      control = .ctl(revisitUninformativeEtas = FALSE))
    expect_equal(unname(.noRev$saem$ueRevisitInfo[["ran"]]), 0L)

    ## Enabled by default, and when it changes no verdict the fit is untouched --
    ## the probe draws no random numbers, so the RNG stream does not shift.
    .rev <- .nlmixr(.m(), .d, "saem", control = .ctl())
    expect_equal(unname(.rev$saem$ueRevisitInfo[["ran"]]), 1L)
    expect_equal(unname(.rev$saem$ueRevisitInfo[["unfroze"]]), 0L)
    expect_equal(unname(.rev$saem$ueRevisitInfo[["froze"]]), 0L)
    expect_equal(fixef(.rev), fixef(.noRev))
  })
})
