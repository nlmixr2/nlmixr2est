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
    # an observation record counts whatever its DV, as a prediction grid's
    .d <- data.frame(ID = 1, TIME = c(0, 0, 1, 2), EVID = c(1, 0, 0, 2), AMT = c(1, 0, 0, 0), DV = c(NA, 1, NA, NA))
    expect_identical(.ipredObsRows(.d), c(FALSE, TRUE, TRUE, FALSE))
    # lower-case names, and no DV column
    .e <- data.frame(id = 1, time = c(0, 1, 2), evid = c(1L, 0L, 0L), amt = c(1, 0, 0))
    expect_identical(.ipredObsRows(.e), c(FALSE, TRUE, TRUE))
  })

  test_that("a prediction grid with no DV takes each subject's etas, and warns for a new one", {
    .fit <- one.compartment.fit.focei
    .ref <- .predIndKeyed(.fit, theo_sd)
    .g <- theo_sd
    .g$DV[.g$EVID == 0] <- NA
    expect_identical(testthat::capture_warnings(.p <- .predIndKeyed(.fit, .g)), character(0))
    expect_equal(.p, .ref, tolerance = 1e-12)
    .g$ID[.g$ID == 12] <- 99
    expect_identical(
      testthat::capture_warnings(.predIndKeyed(.fit, .g)),
      "1 subject(s) not in the fit get etas of 0"
    )
  })

  test_that("newdata without an ID column is predicted for each fitted subject", {
    .et <- rxode2::et(amt = 320) |> rxode2::et(c(0.5, 1, 2, 4, 8, 12, 24))
    .g <- as.data.frame(.et)
    .n <- nrow(.g)
    # each fitted subject's curve, from the grid with its ID
    .fit <- one.compartment.fit.focei
    .all <- .g[rep(seq_len(.n), 12), ]
    .all$ID <- rep(1:12, each = .n)
    .ref <- .predIndKeyed(.fit, .all)
    .p <- as.data.frame(suppressMessages(predict(.fit, .et, level = "individual", atol = 1e-12, rtol = 1e-12)))
    # one solve per fitted subject, in the fit's order
    expect_identical(as.integer(.p$sim.id), rep(1:12, each = 7L))
    expect_equal(.p$ipredSim, unname(.ref[paste(.p$sim.id, .p$time)]), tolerance = 1e-12)
    # a one-subject fit predicts its subject's curve
    .one <- suppressWarnings(suppressMessages(nlmixr2(
      one.compartment,
      theo_sd[theo_sd$ID == 4, ],
      "focei",
      control = foceiControl(print = 0, maxOuterIterations = 0L, covMethod = "")
    )))
    .g$ID <- 4
    .ref1 <- .predIndKeyed(.one, .g)
    .p1 <- .predIndKeyed(.one, .et)
    expect_length(.p1, 7L)
    expect_equal(unname(.p1), unname(.ref1), tolerance = 1e-12)
    expect_true(all(.one$eta$eta.ka != 0))
  })

  test_that("individual predictions match subjects whose labels are not their codes", {
    # IDs 1e5 to 1.2e6, read back as integers in another order
    .d <- theo_sd
    .d$ID <- .d$ID * 1e5
    .fit <- suppressWarnings(suppressMessages(nlmixr2(
      one.compartment,
      .d,
      "focei",
      control = foceiControl(
        print = 0,
        maxOuterIterations = 0L,
        covMethod = "",
        rxControl = rxode2::rxControl(atol = 1e-10, rtol = 1e-10)
      )
    )))
    .tab <- as.data.frame(.fit)
    .want <- stats::setNames(.tab$IPRED, paste(as.numeric(as.character(.tab$ID)), .tab$TIME))
    .n <- .d[order(-.d$ID, seq_len(nrow(.d))), ]
    .n$ID <- as.integer(.n$ID)
    expect_identical(testthat::capture_warnings(.p <- .predIndKeyed(.fit, .n)), character(0))
    names(.p) <- paste(as.numeric(sub(" .*", "", names(.p))), sub("^[^ ]* ", "", names(.p)))
    expect_identical(sort(names(.p)), sort(names(.want)))
    expect_equal(.p[names(.want)], .want, tolerance = 1e-8)
  })

  test_that("saem's individual predictions are its IPRED", {
    .fit <- one.compartment.fit.saem
    .tab <- as.data.frame(.fit)
    .p <- as.data.frame(suppressMessages(predict(.fit, theo_sd, level = "individual")))
    expect_identical(nrow(.p), nrow(.tab))
    expect_equal(.p$ipredSim, .tab$IPRED, tolerance = 1e-6)
  })
})

test_that(".ipredIdLabels() matches a numeric ID to the fit's label by value", {
  .fit <- list(eta = data.frame(ID = factor(c("100000", "2"))))
  expect_identical(.ipredIdLabels(.fit, c(1e5, 2, 3, 2.5e5)), c("100000", "2", "3", "250000"))
  expect_identical(.ipredIdLabels(.fit, c(100000L, 2L)), c("100000", "2"))
  # a label written in scientific notation matches its value
  .fit$eta$ID <- factor(c("1e+05", "2"))
  expect_identical(.ipredIdLabels(.fit, 100000L), "1e+05")
  # an IOV-only fit's labels are its occasion tables'
  .iov <- list(env = list(iov = list(occ = data.frame(ID = factor(c("300000", "300000")), occ = 1:2))))
  expect_identical(.ipredIdLabels(.iov, c(3e5, 4e5)), c("300000", "400000"))
  # text and factor IDs match by their text
  expect_identical(.ipredIdLabels(.fit, c("1e+05", "100000")), c("1e+05", "100000"))
  expect_identical(.ipredIdLabels(.fit, factor("2")), "2")
})

test_that(".ipredHasId() looks for a subject column in the data frame", {
  expect_true(.ipredHasId(data.frame(Id = 1, time = 0)))
  expect_false(.ipredHasId(data.frame(time = 0)))
  # an event table of one subject has none
  expect_false(.ipredHasId(rxode2::et(amt = 1) |> rxode2::et(1)))
  expect_true(.ipredHasId(rxode2::et(amt = 1, id = 1:2) |> rxode2::et(1)))
})
