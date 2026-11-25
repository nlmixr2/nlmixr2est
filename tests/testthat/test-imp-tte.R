nmTest({
  skip_if_not(.nlmixr2estRxode2HasTte(), "rxode2 without tte() endpoints")

  test_that("imp-family fits recover a Weibull repeated time-to-event model (#1467)", {
    skip_on_cran()
    .weib <- function() {
      ini({
        tlam <- log(0.3)
        tk <- log(1.4)
        eta.lam ~ 0.2
      })
      model({
        lam <- exp(tlam + eta.lam)
        k <- exp(tk)
        haz <- k * lam * (lam * (time + 1e-8))^(k - 1)
        ev ~ tte(haz)
      })
    }
    .s <- withr::with_seed(42, rxode2::rxSolve(.weib(), rxode2::et(c(0, 24)) |> rxode2::et(id = 1:100)))
    .d <- as.data.frame(.s)[, c("id", "time", "sim")]
    names(.d)[3] <- "dv"
    # start away from the truth so the M-step has to move the shape
    .start <- .weib |> rxode2::ini(tlam = log(0.1), tk = 0, eta.lam ~ 0.1)
    .ref <- suppressMessages(.nlmixr(.start, .d, est = "focei", foceiControl(print = 0, calcTables = FALSE)))
    for (.est in c("imp", "impmap", "qrpem")) {
      .fit <- suppressMessages(.nlmixr(.start, .d, est = .est, impmapControl(print = 0, calcTables = FALSE)))
      expect_equal(
        unname(.fit$theta[c("tlam", "tk")]),
        unname(.ref$theta[c("tlam", "tk")]),
        tolerance = 0.1,
        label = .est
      )
      expect_equal(unname(.fit$omega[1, 1]), unname(.ref$omega[1, 1]), tolerance = 0.5, label = .est)
    }
  })
})
