## Regression for the neonate-style data path: .vaeDataPrep must not crash when
## the data has neither an AMT nor an EVID column (dose-free observations). The
## former `ifelse(is.na(d$AMT) | d$AMT == 0, ...)` yielded a length-0 vector when
## d$AMT was NULL, giving "replacement has 0 rows, data has N". Fast: no fit.

nmTest({
  test_that(".vaeDataPrep handles data with no AMT and no EVID columns", {
    theo <- function() {
      ini({
        lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
        add.err <- 0.7
      })
      model({
        ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V
        cp ~ add(add.err)
      })
    }
    ui <- rxode2::assertRxUi(theo)

    d <- nlmixr2data::theo_sd
    d <- d[d$EVID == 0, , drop = FALSE]   # observations only
    d$AMT <- NULL
    d$EVID <- NULL
    expect_false("AMT" %in% names(d))
    expect_false("EVID" %in% names(d))

    prep <- .vaeDataPrep(ui, d)
    expect_equal(prep$N, length(unique(d$ID)))
    expect_equal(prep$Nobs, nrow(d))       # every row treated as an observation
    expect_true(all(vapply(prep$subj, function(s) all(s$ev$EVID == 0L), logical(1))))
  })
})
