nmTest({
  test_that("uninformative etas", {

    data <- data.frame(ID = c(11, 11, 11, 11, 11, 11, 11, 12, 12, 12,
                              12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 21, 21, 21, 21,
                              21, 21, 21, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23),
                       TIME = c(0, 0.05, 0.25, 0.5, 1, 3, 5, 0, 0.05, 0.25, 0.5,
                                1, 3, 5, 0, 0.05, 0.25, 0.5, 1, 3, 5, 8, 0, 0.25, 0.5, 1, 3,
                                5, 8, 0, 0.25, 0.5, 1, 3, 5, 8, 0, 0.25, 0.5, 1, 3, 5, 8),
                       DV = c(NA,2017.85, 1323.74, 792.5, 822.72, 36.27, 3.33, NA, 1702, 1290.75,
                              1095.95, 907.6, 125.44, 14.44, NA, 1933.04, 1242.43, 661.22,
                              193.52, 1.75, NA, NA, NA, 706.58, 1063.14, 2257.62, 941.33, 629.69,
                              100, NA, 1462.95, 2217.76, 2739.5, 705.3, 108.47, 8.75, NA, 211.66,
                              467.23, 174.24, 153.6, 27.07, 2.81),
                       AMT = c(1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0,
                               0, 0, 5, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0),
                       EVID = c(1,0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                                1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
                       CMT = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2,
                               2, 2, 2, 2, 2),
                       DOSE = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
                       ROUTE = c(1, 1, 1, 1, 1, 1, 1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2))

    # First fit with if statement to indicate diff IV/PO
    modA <- function() {
      ini({
        tka <- 0.45   # Log Ka
        tcl <- -7     # Log Cl
        tv  <- -8     # Log V
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

    f <-.nlmixr(modA, data, "saem", control=saemControl(print=0))

    expect_true(all(f$etaObf$eta.ka[1:3] == 0))

    f2 <-.nlmixr(modA, data, "saem", control=saemControl(handleUninformativeEtas=FALSE, print=0))

    expect_false(all(f2$etaObf$eta.ka[1:3] == 0))

    # First fit with if statement to indicate diff IV/PO
    modB <- function() {
      ini({
        tka <- 0.45   # Log Ka
        tcl <- -7     # Log Cl
        tv  <- -8     # Log V
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        prop.sd <- 0.7
      })
      model({ # non mu referenced non-informative eta
        ka <- tka * exp(eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl/v * center
        cp <- center/v
        cp ~ prop(prop.sd)
      })
    }

    f <-.nlmixr(modB, data, "saem", control=saemControl(print=0))

    expect_true(all(f$etaObf$eta.ka[1:3] == 0))

    f2 <-.nlmixr(modB, data, "saem", control=saemControl(handleUninformativeEtas=FALSE, print=0))

    expect_false(all(f2$etaObf$eta.ka[1:3] == 0))

  })
})
