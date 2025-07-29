nmTest({
  test_that("test eta only focei models work correctly", {
    skip_if_not(utils::packageVersion("rxode2") >= "4.0.0")
    f <- function() {
      ini({
        eta.cl ~ 0.0496342898304385
        eta.v ~ 0.0150288456328789
        eta.ka ~ 0.175010734212212
      })
      model({
        base1 = D_ETA1 * (-O_ETA1 + eta.cl)
        base2 = D_ETA2 * (-O_ETA2 + eta.v)
        base3 = D_ETA3 * (-O_ETA3 + eta.ka)
        BASE_TERMS = base1 + base2 + base3
        IPRED = BASE_TERMS + OPRED
        err1 = D_EPSETA_1_1 * (-O_ETA1 + eta.cl)
        err2 = D_EPSETA_1_2 * (-O_ETA2 + eta.v)
        err3 = D_EPSETA_1_3 * (-O_ETA3 + eta.ka)
        BASE_ERROR = err1 + err2 + err3
        R2 = (BASE_ERROR + O_ResVar)
        y = IPRED
        y ~ add(R2) + var()
      })
    }

    f <- f()

    suppressMessages(expect_error(f$foceiModel, NA))

  })
})
