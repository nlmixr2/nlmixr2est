nmTest({

  test_that("Inner test", {
    ev <- eventTable() %>%
      add.sampling(c(
        95.99, 119.99, 143.99, 144.25, 144.5, 144.75,
        145, 145.5, 146, 146.5, 147, 148, 150, 152, 156, 160, 164, 167.99,
        191.99, 215.99, 216.25, 216.5, 216.75, 217, 217.5, 218, 218.5, 219,
        220, 222, 224, 228, 232, 236, 240, 252, 264, 276, 288
      )) %>%
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

    fitPi <- suppressMessages(nlmixr(
      m1, w7,
      est="focei",
      foceiControl(
        etaMat = ETA,
        maxOuterIterations = 0, maxInnerIterations = 0,
        covMethod = ""
      )
    ))

    expect_equal(418.935, round(fitPi$objective, 3))
  })

})
