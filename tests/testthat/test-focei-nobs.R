nmTest({
  test_that("FOCEi counts only observations when EVID 2 and 9 are present (#1342)", {
    oneCompartment <- function() {
      ini({
        lcl <- log(1)
        lv <- log(10)
        eta.cl ~ 0.1
        add.sd <- 0.5
      })
      model({
        cl <- exp(lcl + eta.cl)
        v <- exp(lv)
        d/dt(central) <- -(cl / v) * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    }

    late <- data.frame(
      ID = 1L,
      TIME = c(2, 3, 5, 8),
      DV = c(NA, 9, 7, 5),
      AMT = c(100, 0, 0, 0),
      EVID = c(1L, 0L, 0L, 0L),
      CMT = 1L
    )
    # rxode2 inserts EVID 9 for subject 1. Subject 2 starts with EVID 2.
    other <- rbind(transform(late[1, ], TIME = 0, AMT = 0, EVID = 2L), late)
    dat <- rbind(late, transform(other, ID = 2L))

    fit <- suppressMessages(suppressWarnings(
      nlmixr(oneCompartment(), dat, "focei", control = foceiControl(maxOuterIterations = 0, covMethod = "", print = 0))
    ))
    expect_true(all(c(2L, 9L) %in% fit$dataSav$EVID))
    expect_equal(fit$nobs, sum(dat$EVID == 0L))
    expect_equal(stats::nobs(fit), attr(logLik(fit), "nobs"))
  })
})
