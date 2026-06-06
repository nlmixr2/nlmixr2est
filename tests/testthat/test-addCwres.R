nmTest({
  test_that("addCwres", {
    one.compartment <- function() {
      ini({
        tka <- log(1.57)
        tcl <- log(2.72)
        tv <- log(31.5)
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    suppressMessages(
      fitNoEta <-.nlmixr(one.compartment, theo_sd, est="focei", control = list(print=0))
    )
    expect_true(inherits(fitNoEta$parHistData, "data.frame"))
    expect_error(
      suppressMessages(addCwres(fitNoEta)),
      regexp = "cannot add CWRES to a model without etas"
    )
  })

  test_that("cwres issue in 3.0.3", {
    skip_if_not(file.exists(test_path("test-cwres-3.0.3.rds")))

    data <- readRDS(test_path("test-cwres-3.0.3.rds"))

    test_model <- function() {
      ini({
        lcl  <- log(3)
        lvc   <- log(40)
        prop.err <- 0.1
        eta.cl ~ 0.1
        eta.vc  ~ 0.1
        WT_Cl <- fix(0.75)
        ClCrEff <- 1
      })
      model({
        cl<- exp(lcl + eta.cl + WT_Cl * log(Weight/81.60) + ClCrEff * log(ClCr/77.73))
        vc  <- exp(lvc + eta.vc)
        d/dt(A_cen) = - cl/vc * A_cen
        cp = A_cen/vc
        cp ~ prop(prop.err)
      })
    }

    test_run001 <- .nlmixr(test_model(), data, "saem", control = saemControlFast)

    expect_error(suppressMessages(addCwres(test_run001)), NA)
  })
})
