# The FOCEi inner-problem harness of the nonparametric engines (R/npInner.R).

nmTest({
  test_that(".npInnerSetup() solves with the control's eventSens and indTolRelax", {
    npMod <- function() {
      ini({
        tka <- log(1.5)
        tcl <- log(2.7)
        tv <- log(32)
        eta.cl ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    ui <- rxode2::assertRxUi(npMod)
    dat <- nlmixr2data::theo_sd
    N <- length(unique(dat$ID))
    for (.v in list(list(eventSens = "fd", indTolRelax = FALSE), list(eventSens = "jump", indTolRelax = TRUE))) {
      ctl <- do.call(npagControl, .v)
      expect_identical(ctl[c("eventSens", "indTolRelax")], .v)
      .env <- .npInnerSetup(ui, dat, matrix(0, N, 1L), ctl)
      .npInnerFree()
      expect_identical(.env$control[c("eventSens", "indTolRelax")], .v)
    }
  })
})
