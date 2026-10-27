nmTest({
  test_that("predict", {
    # Use centralized fit from helper-fits.R
    ## The fit was performed with focei, maxOuterIterations = 0L, and print = 0

    fit <- one.compartment.fit.focei

    # given no solving options both entry points inherit the fit's own rxControl
    # (rtol=1e-3); pass it explicitly here so the pairing is visible
    md <- suppressMessages(do.call("predict", c(list(fit, theo_md), fit$control$rxControl)))

    md2 <- .nlmixr(fit, theo_md, "predict")

    expect_equal(as.data.frame(md), as.data.frame(md2), tolerance = 1e-4)

    # with no solving options both entry points use the fit's own rxControl
    expect_equal(
      as.data.frame(suppressMessages(predict(fit, theo_md))),
      as.data.frame(.nlmixr(fit, theo_md, "predict")),
      tolerance = 1e-4
    )

    # so comparing against rxControl() defaults needs them passed on both sides
    md <- suppressMessages(
      do.call("predict", c(list(fit, theo_md), rxode2::rxControl()))
    )

    md2 <- .nlmixr(fit, theo_md, "predict", control = rxode2::rxControl())

    expect_equal(as.data.frame(md), as.data.frame(md2), tolerance = 1e-4)

    ipred <- suppressMessages(predict(fit, theo_sd, level = "individual"))
    expect_equal(ipred$ipredSim, fit$IPRED, tolerance = 1e-4)

    # Test explicit population level (default)
    ppred <- suppressMessages(predict(fit, theo_sd, level = "population"))
    expect_true("pred" %in% names(ppred))

    # Test numeric level=0 (population)
    ppred0 <- suppressMessages(predict(fit, theo_sd, level = 0))
    expect_true("pred" %in% names(ppred0))

    # Test numeric level=1 (individual)
    ipred1 <- suppressMessages(predict(fit, theo_sd, level = 1))
    expect_equal(ipred1$ipredSim, fit$IPRED, tolerance = 1e-4)

    # Test alias level="pred" (population)
    ppredAlias <- suppressMessages(predict(fit, theo_sd, level = "pred"))
    expect_true("pred" %in% names(ppredAlias))

    # Test alias level="ppred" (population)
    ppredAlias2 <- suppressMessages(predict(fit, theo_sd, level = "ppred"))
    expect_true("pred" %in% names(ppredAlias2))

    # Test alias level="ipred" (individual)
    ipredAlias <- suppressMessages(predict(fit, theo_sd, level = "ipred"))
    expect_equal(ipredAlias$ipredSim, fit$IPRED, tolerance = 1e-4)

    # Test invalid numeric level throws error
    expect_error(
      suppressMessages(predict(fit, theo_sd, level = 2)),
      "level numeric must be 0 \\(population\\) or 1 \\(individual\\)"
    )

    # Test individual predictions with new data (theo_md)
    ipredNewData <- suppressMessages(predict(fit, theo_md, level = "individual"))
    expect_true("ipredSim" %in% names(ipredNewData))
  })
})

# individual predictions of `data`, named by "id time"; tight tolerances, so
# a solve of fewer subjects agrees with one of all of them
.predIndKeyed <- function(fit, data) {
  .p <- as.data.frame(suppressMessages(predict(fit, data, level = "individual", atol = 1e-12, rtol = 1e-12)))
  stats::setNames(.p$ipredSim, paste(.p$id, .p$time))
}

nmTest({
  test_that("individual predictions take each subject's own etas, in any order or subset", {
    for (.est in c("focei", "saem")) {
      .fit <- get(paste0("one.compartment.fit.", .est))
      .ref <- .predIndKeyed(.fit, theo_sd)
      expect_length(.ref, sum(theo_sd$EVID == 0))
      # subjects in descending order
      .desc <- theo_sd[order(-theo_sd$ID, seq_len(nrow(theo_sd))), ]
      expect_equal(.predIndKeyed(.fit, .desc)[names(.ref)], .ref, tolerance = 1e-12, label = .est)
      # two subjects only
      .sub <- .predIndKeyed(.fit, theo_sd[theo_sd$ID %in% c(3, 5), ])
      expect_length(.sub, sum(theo_sd$EVID == 0 & theo_sd$ID %in% c(3, 5)))
      expect_equal(.sub, .ref[names(.sub)], tolerance = 1e-12, label = .est)
      # subject 3 with its dose only: nothing to predict for it, no warning
      .dd <- theo_sd[!(theo_sd$ID == 3 & theo_sd$EVID == 0), ]
      .w <- testthat::capture_warnings(.dose <- .predIndKeyed(.fit, .dd))
      expect_identical(.w, character(0), label = .est)
      expect_false(any(startsWith(names(.dose), "3 ")))
      expect_equal(.dose, .ref[names(.dose)], tolerance = 1e-12, label = .est)
    }
  })

  test_that("a subject the fit did not estimate is predicted with etas of 0, with a warning", {
    for (.est in c("focei", "saem")) {
      .fit <- get(paste0("one.compartment.fit.", .est))
      .ref <- .predIndKeyed(.fit, theo_sd)
      .dn <- theo_sd
      .dn$ID[.dn$ID == 12] <- 99
      .w <- testthat::capture_warnings(.new <- .predIndKeyed(.fit, .dn))
      expect_identical(.w, "1 subject(s) not in the fit get etas of 0", label = .est)
      .pop <- as.data.frame(suppressMessages(predict(.fit, .dn, level = "population", atol = 1e-12, rtol = 1e-12)))
      .pop <- stats::setNames(.pop$pred, paste(.pop$id, .pop$time))
      .k99 <- startsWith(names(.new), "99 ")
      expect_equal(sum(.k99), sum(theo_sd$EVID == 0 & theo_sd$ID == 12))
      expect_equal(.new[.k99], .pop[names(.new)[.k99]], tolerance = 1e-8, label = .est)
      # the other subjects keep their own etas
      expect_equal(.new[!.k99], .ref[names(.new)[!.k99]], tolerance = 1e-12, label = .est)
    }
  })

  test_that("a fit with a dose-only subject predicts its other subjects as IPRED", {
    .dd <- theo_sd[!(theo_sd$ID == 3 & theo_sd$EVID == 0), ]
    .fit <- suppressWarnings(suppressMessages(nlmixr2(
      one.compartment,
      .dd,
      "focei",
      control = foceiControl(
        print = 0,
        maxOuterIterations = 0L,
        covMethod = "",
        rxControl = rxode2::rxControl(atol = 1e-10, rtol = 1e-10)
      )
    )))
    expect_identical(as.character(.fit$eta$ID), as.character(c(1:2, 4:12)))
    # subject 3 has no observations, so it is not predicted and gives no warning
    .w <- testthat::capture_warnings(.p <- .predIndKeyed(.fit, .dd))
    expect_identical(.w, character(0))
    .tab <- as.data.frame(.fit)
    .tab <- .tab[!is.na(.tab$IPRED), ]
    .key <- paste(.tab$ID, .tab$TIME)
    expect_identical(sort(names(.p)), sort(.key))
    expect_equal(unname(.p[.key]), .tab$IPRED, tolerance = 1e-8)
  })

  test_that(".ipredObsRows() finds the observations of data, and of an event table", {
    .d <- data.frame(ID = 1, TIME = c(0, 0, 1, 2), EVID = c(1, 0, 0, 2), AMT = c(1, 0, 0, 0), DV = c(NA, 1, NA, NA))
    expect_identical(.ipredObsRows(.d), c(FALSE, TRUE, FALSE, FALSE))
    # lower-case names, and no DV column: every observation record counts
    .e <- data.frame(id = 1, time = c(0, 1, 2), evid = c(1L, 0L, 0L), amt = c(1, 0, 0))
    expect_identical(.ipredObsRows(.e), c(FALSE, TRUE, TRUE))
  })
})
