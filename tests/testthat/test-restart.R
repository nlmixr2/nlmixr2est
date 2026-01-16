nmTest({
  test_that("restart resid", {

    lobo  <- function() {
      ini({
        lkng <- log(0.02)
        ltau <- log(c(1, 34.1, 500))
        lec50 <- log(c(1, 50, 2000))
        kmax <- 0.01
        propSd <-  c(0, 0.3)
        addSd <- c(0, 50)
      })
      model({
        kng <- exp(lkng)
        tau <- exp(ltau)
        taulast <- tau
        ec50 <- exp(lec50)
        edrug <- kmax * cp/(ec50 + cp)
        tumor(0) <- tumor0
        d/dt(transit1) <- (edrug - transit1)/tau
        d/dt(transit2) <- (transit1 - transit2)/tau
        d/dt(transit3) <- (transit2 - transit3)/tau
        d/dt(transitlast) <- transit3/tau - transitlast/taulast
        d/dt(tumor) <- kng * tumor - transitlast * tumor
        tumor ~ prop(propSd) + add(addSd)
      })
    }

    prepfit <- readRDS(test_path("test-restart.rds"))

    expect_error(.nlmixr(lobo, data=prepfit, "focei", control=foceiControl(print=0)), NA)

  })
})
