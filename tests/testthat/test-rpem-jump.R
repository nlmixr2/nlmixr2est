# Dose-parameter ("jump") sensitivities: a bioavailability F scales the depot bolus,
# so its gradient is an event/dosing-parameter sensitivity.  RPEM estimates the
# structural F via its finite-difference structural M-step (rpemMstepBeta re-solves
# through the scaled dose, which is inherently jump-correct), and the FOCEI eval that
# builds the fit object computes F's SE with eventSens="jump" (set explicitly in
# .rpemBuildFit).  Both the estimate and the SE must match FOCEI.

test_that("RPEM handles a dose (bioavailability) parameter and its jump sensitivity", {
  skip_on_cran()
  skip_on_ci()  # heavy: FOCEI + RPEM fits

  sim <- rxode2::rxode2({ ka <- exp(tka + eta); cl <- exp(tcl); v <- exp(tv)
                          f(depot) <- exp(tf); cp <- linCmt() })
  set.seed(12); nsub <- 40L; obsT <- seq(0.5, 24, by = 1.5); tfTrue <- log(0.7)
  etasTrue <- rnorm(nsub, 0, sqrt(0.3))
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(sim, params = c(tka = 0.45, tcl = 1, tv = 3.45, tf = tfTrue,
                                         eta = etasTrue[i]),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, 0.1)
    d
  }))

  mod <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); tf <- log(0.5); add.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
            f(depot) <- exp(tf); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ff <- suppressMessages(nlmixr2(mod, dat, est = "focei",
                                 control = foceiControl(print = 0, calcTables = FALSE)))
  rf <- suppressMessages(nlmixr2(mod, dat, est = "rpem",
                                 control = rpemControl(nGauss = 200L, nMH = 50000L, mhBurn = 5000L,
                                                       niter = 30L, collect = 12L, seed = 100L, cores = 4L)))

  .r <- rf$parFixedDf$Estimate; names(.r) <- rownames(rf$parFixedDf)
  .f <- ff$parFixedDf$Estimate; names(.f) <- rownames(ff$parFixedDf)
  # the dose parameter (and the others) recover the truth and match FOCEI
  expect_equal(unname(.r["tf"]), tfTrue, tolerance = 0.1)
  expect_equal(unname(.r["tf"]), unname(.f["tf"]), tolerance = 0.05)
  expect_equal(unname(.r["tka"]), unname(.f["tka"]), tolerance = 0.05)
  expect_equal(unname(.r["add.sd"]), unname(.f["add.sd"]), tolerance = 0.05)
  # the dose-parameter SE (from the eval's jump sensitivities) matches FOCEI: if the
  # event sensitivities were wrong the F SE would not agree.
  .rse <- rf$parFixedDf["tf", "SE"]; .fse <- ff$parFixedDf["tf", "SE"]
  expect_true(is.finite(.rse) && .rse > 0)
  expect_equal(unname(.rse), unname(.fse), tolerance = 0.2)
})
