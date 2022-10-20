nmTest({
  test_that("non-mu simple test", {

    # Data' purpose illustrates the error and my data set
    df <- data.frame(
      ID = c(rep(1, 6), rep(2, 6)),
      TIME = c(0.00, 12.11, 18.41, 23.89, 36.00, 43.51, 0.00, 12.00, 20.00, 24.00, 36.80, 45.00),
      AMT = c(1000, 1000, NA, 1000, 1000, NA, 1000, 2000, NA, 1000, 1000, NA),
      DUR = c(2.5, 2.5, NA, 2.5, 2.5, NA, 2.5, 2.5, NA, 2.5, 2.5, NA),
      DV = c(NA, NA, 3.0, NA, NA, 9.6, NA, NA, 7.0, NA, NA, 2.8),
      WT = c(rep(55, 6), rep(48, 6))
    ) |>
      dplyr::mutate(EVID = ifelse(is.na(DV), 1, 0))

    fun <- function() {
      ini({
        tvCl <- c(0, 4, Inf)
        tvVc <- c(0, 48, Inf)
        eta.Vc ~ 0.62
        prop.sd <- 0.051529
      })
      model({
        Cl <- tvCl
        Vc <- tvVc*(WT/70)*exp(eta.Vc)
        # dynamical system
        linCmt() ~ prop(prop.sd)
      })
    }

    fit <- nlmixr2(fun, df, list(print=0), est="posthoc")
    
    expect_error(fit$dataMergeInner, NA)
    
    tmp <- fit$dataMergeInner
    
    # Should have llikObs
    expect_true("llikObs" %in% names(tmp))
    
    expect_true(all(names(fit$etaSE) == c("ID", "eta.Vc")))
    
    expect_true(all(names(fit$etaRSE) == c("ID", "rse(eta.Vc)%")))

  })
})
