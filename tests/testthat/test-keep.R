nmTest({

  # Use centralized model from helper-models.R and add keep/drop
  one.compartment.keep.drop <- function() {
    ini({
      tka <- 0.45 # Log Ka
      tcl <- 1 # Log Cl
      tv <- 3.45    # Log V
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) = -ka * depot
      d/dt(center) = ka * depot - cl / v * center
      cp = center / v
      cp ~ add(add.sd)
    })
    keep = c("WT")
    drop = c("depot")
  }

  d <- theo_sd
  d$WT2 <- d$WT + 0.5
  d$WT3 <- d$WT + 0.4

  for (est in c("fo", "foi", "foce", "focei", "saem", "posthoc")) {
    test_that(paste0("keep/drop in ", est), {
      fitF <- suppressWarnings(suppressMessages(nlmixr(one.compartment.keep.drop, d, est=est)))
      fitF2 <- suppressWarnings(suppressMessages(nlmixr(one.compartment.keep.drop, d, est=est, table=list(keep="WT2", drop="center"))))
      expect_true(any(names(fitF) == "WT"))
      expect_true(!any(names(fitF) == "WT2"))
      expect_true(!any(names(fitF) == "WT3"))
      expect_true(!any(names(fitF) == "depot"))
      expect_true(any(names(fitF) == "center"))
      expect_true(any(names(fitF) == "dosenum"))
      expect_true(!any(names(fitF) == "rxLambda"))
      expect_true(!any(names(fitF) == "rxYj"))

      expect_true(!any(names(fitF2) == "WT"))
      expect_true(any(names(fitF2) == "WT2"))
      expect_true(!any(names(fitF) == "WT3"))
      expect_true(any(names(fitF2) == "depot"))
      expect_true(!any(names(fitF2) == "center"))
      expect_true(any(names(fitF2) == "dosenum"))
      expect_true(!any(names(fitF2) == "rxLambda"))
      expect_true(!any(names(fitF2) == "rxYj"))
    })
  }

})
