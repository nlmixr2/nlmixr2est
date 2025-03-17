nmTest({
  test_that("character estimation works", {
    one_compartment_textcov <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        cllow <- 0.001
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      # and a model block with the error sppecification and model specification
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl + cllow*(LowID == "Yes"))
        v <- exp(tv + eta.v)
        d/dt(depot) = -ka * depot
        d/dt(center) = ka * depot - cl / v * center
        cp = center / v
        cp ~ add(add.sd)
      })
    }

    fitdata <- theo_sd
    fitdata$ID <- paste("theo", fitdata$ID)
    fitdata$LowID <- ifelse(as.numeric(theo_sd$ID) < 7, "Yes", "No")
    fit <- .nlmixr(one_compartment_textcov, fitdata, est="focei", control = foceiControl(print = 0))

    .cllow <- fit$theta["cllow"]

    expect_true(inherits(fit$ID, "factor"))

    expect_equal(levels(fit$ID),
                 c("theo 1", "theo 2", "theo 3", "theo 4", "theo 5", "theo 6", "theo 7",
                   "theo 8", "theo 9", "theo 10", "theo 11", "theo 12"))

    expect_true(inherits(fit$LowID, "factor"))
    expect_equal(levels(fit$LowID), c("No", "Yes"))

    fitdata$LowID <- factor(fitdata$LowID, c("Yes", "No"))
    fit <- .nlmixr(one_compartment_textcov, fitdata, est = "focei", control = foceiControl(print = 0))

    expect_equal(.cllow, fit$theta["cllow"])
  })
})
