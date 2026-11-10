nmTest({
  test_that("saem reads its predictions back the same on one and on several threads", {
    skip_on_cran()
    # the readback after each bulk solve runs over the subjects in parallel;
    # each subject's rows have a fixed place, so the fit cannot depend on the
    # thread count
    .mod <- function() {
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
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    .fit <- function(cores) {
      suppressMessages(nlmixr2(
        .mod,
        nlmixr2data::theo_sd,
        est = "saem",
        control = saemControl(
          print = 0,
          nBurn = 10,
          nEm = 10,
          calcTables = FALSE,
          seed = 99,
          rxControl = rxode2::rxControl(cores = cores)
        )
      ))
    }
    .f1 <- .fit(1L)
    for (.cores in c(2L, 4L)) {
      .fn <- .fit(.cores)
      expect_identical(.fn$theta, .f1$theta, label = paste(.cores, "threads"))
      expect_identical(.fn$omega, .f1$omega, label = paste(.cores, "threads"))
      expect_identical(.fn$eta, .f1$eta, label = paste(.cores, "threads"))
    }
  })
})
