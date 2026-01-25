nmTest({
  test_that("focei complex event info", {

    pheno <- function() {
      ini({
        tcl <- log(0.008) # typical value of clearance
        tv <-  log(0.6)   # typical value of volume
        max_dose <- 5
        ## var(eta.cl)
        eta.cl + eta.v ~ c(1,
                           0.01, 1) ## cov(eta.cl, eta.v), var(eta.v)
        # interindividual variability on clearance and volume
        add.err <- 0.1    # residual variability
      })
      model({
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        fest <- max_dose/DOSE
        if (fest > 1) fest <- 1 # error is here
        d/dt(A1) = - ke * A1
        f(A1) <- fest
        cp = A1 / v
        cp ~ add(add.err)
      })
    }

    f <- suppressMessages(pheno())
    expect_error(f$foceiModel, NA)
  })

  test_that("Inner test", {
    ev <- eventTable() |>
      add.sampling(c(
        95.99, 119.99, 143.99, 144.25, 144.5, 144.75,
        145, 145.5, 146, 146.5, 147, 148, 150, 152, 156, 160, 164, 167.99,
        191.99, 215.99, 216.25, 216.5, 216.75, 217, 217.5, 218, 218.5, 219,
        220, 222, 224, 228, 232, 236, 240, 252, 264, 276, 288
      )) |>
      add.dosing(dose = 60000, start.time = 72, nbr.doses = 7, dosing.interval = 24)

    dv <- c(
      263.6, 164.7, 287.3, 1248.7, 1211.5, 1017.7, 1690.1, 1029.8,
      890.7, 598.4, 1009.3, 1159.8, 742.2, 724.6, 728.2, 509.7, 243.1,
      259.9, 242.2, 281.4, 1500.1, 1281.4, 1200.2, 1378.8, 1373.2,
      582.9, 960.2, 720.3, 852.6, 950.3, 654.7, 402.5, 456, 346.5,
      268.2, 134.2, 42.6, 25.9, 14.6
    )

    m1 <- function() {
      ini({
        tcl <- 1.6
        tv <- 4.5
        eta.cl ~ 0.1
        eta.v ~ 0.1
        prop.sd <- sqrt(0.1)
      })
      model({
        CL <- exp(tcl + eta.cl)
        V <- exp(tv + eta.v)
        C2 <- centr / V
        d/dt(centr) <- -CL * C2
        C2 ~ prop(prop.sd)
      })
    }

    w7 <- data.frame(ev$get.EventTable())
    w7$DV <- NA
    w7$DV[which(is.na(w7$amt))] <- dv
    w7$ID <- 1

    ETA <- matrix(c(-0.147736086922763, -0.294637022436797), ncol = 2)

    fitPi <- .nlmixr(
      m1, w7,
      est="focei",
      foceiControl(
        etaMat = ETA,
        maxOuterIterations = 0, maxInnerIterations = 0,
        covMethod = ""
      )
    )

    expect_equal(418.935, round(fitPi$objective, 3))
  })

  test_that("boundary value is not triggered by bounds on both sides of zero (#318)", {

    one.compartment <- function() {
      ini({
        tka <- c(-6, -4, 2)
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)*100
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }

    fit <- .nlmixr(one.compartment, theo_sd, est = "focei", control = list(print=0))
    # SE being present indicates that the covariance matrix was estimated
    expect_true("SE" %fin% names(fit$parFixedDf))

    # Also make sure that it correctly identifies mu-ref
    expect_equal(c("tka", "tcl", "tv", "add.sd"), row.names(fit$parFixedDf))
  })
})
