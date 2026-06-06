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

    fit <-.nlmixr(fun, df, list(print=0), est="posthoc")

    expect_error(fit$dataMergeInner, NA)
    expect_error(fit$fitMergeInner, NA)

    tmp <- fit$dataMergeInner

    # Should have llikObs
    expect_true("nlmixrLlikObs" %in% names(tmp))

    expect_true(all(names(fit$etaSE) == c("ID", "eta.Vc")))

    expect_true(all(names(fit$etaRSE) == c("ID", "rse(eta.Vc)%")))

  })

  test_that("another merge issue", {
    dat <- xgxr::case1_pkpd |>
      dplyr::rename(DV=LIDV) |>
      dplyr::filter(CMT %in% 1:2) |>
      dplyr::filter(TRTACT != "Placebo")

    doses <- unique(dat$DOSE)
    nid <- 3 # 7 ids per dose group
    dat2 <- do.call("rbind",
                    lapply(doses, function(x) {
                      ids <- dat |>
                        dplyr::filter(DOSE == x) |>
                        dplyr::reframe(ids=unique(ID)) |>
                        dplyr::pull()
                      ids <- ids[seq(1, nid)]
                      dat |>
                        dplyr::filter(ID %in% ids)
                    }))

    # Use centralized model from helper-models.R
    cmt2 <- two.compartment

    cmt2fit.logn <-
      .nlmixr(
        cmt2, dat2, "posthoc",
        control=list(print=0),
        table=tableControl(cwres=TRUE, npde=TRUE)
      )

    expect_error(cmt2fit.logn$dataMergeLeft, NA)
    expect_error(cmt2fit.logn$fitMergeLeft, NA)
    expect_true(any(names(cmt2fit.logn$dataMergeLeft) == "nlmixrLlikObs"))

    # Now force an error

    .llikObs <- c(cmt2fit.logn$env$llikObs, 10)
    assign("llikObs", .llikObs, envir=cmt2fit.logn$env)

    expect_warning(cmt2fit.logn$dataMergeLeft)
    expect_warning(cmt2fit.logn$fitMergeLeft)

    .dat <- suppressWarnings(cmt2fit.logn$dataMergeLeft)
    expect_false(any(names(.dat) == "nlmixrLlikObs"))

  })
})
