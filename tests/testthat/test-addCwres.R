nmTest({
  test_that("addCwres", {
    one.compartment <- function() {
      ini({
        tka <- log(1.57)
        tcl <- log(2.72)
        tv <- log(31.5)
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    suppressMessages(
      fitNoEta <-.nlmixr(one.compartment, theo_sd, est="focei", control = list(print=0))
    )
    expect_true(inherits(fitNoEta$parHistData, "data.frame"))
    expect_error(
      suppressMessages(addCwres(fitNoEta)),
      regexp = "cannot add CWRES to a model without etas"
    )
  })

  test_that("cwres issue in 3.0.3", {
    skip_if_not(file.exists(test_path("test-cwres-3.0.3.rds")))

    data <- readRDS(test_path("test-cwres-3.0.3.rds"))

    test_model <- function() {
      ini({
        lcl  <- log(3)
        lvc   <- log(40)
        prop.err <- 0.1
        eta.cl ~ 0.1
        eta.vc  ~ 0.1
        WT_Cl <- fix(0.75)
        ClCrEff <- 1
      })
      model({
        cl<- exp(lcl + eta.cl + WT_Cl * log(Weight/81.60) + ClCrEff * log(ClCr/77.73))
        vc  <- exp(lvc + eta.vc)
        d/dt(A_cen) = - cl/vc * A_cen
        cp = A_cen/vc
        cp ~ prop(prop.err)
      })
    }

    test_run001 <- .nlmixr(test_model(), data, "saem", control = saemControlFast)

    expect_error(suppressMessages(addCwres(test_run001)), NA)
  })

  test_that(".foceiObjfWithoutCwres() flags the fits addCwres() cannot fix later", {
    .objDf <- function(rn, objf) {
      .d <- data.frame(OBJF = objf, AIC = objf, BIC = objf,
                       "Log-likelihood" = -objf / 2, check.names = FALSE)
      row.names(.d) <- rn
      .d
    }
    # these are the rows addCwres() itself would add, so it could not add them twice
    expect_true(.foceiObjfWithoutCwres(list(objDf = .objDf("FOCEi", 100))))
    expect_true(.foceiObjfWithoutCwres(list(objDf = .objDf("lFOCEi", 100))))
    expect_true(.foceiObjfWithoutCwres(list(objDf = .objDf("FOCE", 100))))
    # an uncalculated objective function is replaced by addCwres(), not appended
    expect_false(.foceiObjfWithoutCwres(list(objDf = .objDf("FOCEi", NA_real_))))
    # a quadrature/FO row does not collide with the row addCwres() adds
    expect_false(.foceiObjfWithoutCwres(list(objDf = .objDf("Laplace", 100))))
    expect_false(.foceiObjfWithoutCwres(list(objDf = .objDf("FO", 100))))
    expect_false(.foceiObjfWithoutCwres(list(objDf = NULL)))
  })

  test_that("addCwres() works on a focei fit whose table skipped CWRES", {
    one.compartment <- function() {
      ini({
        tka <- log(1.57)
        tcl <- log(2.72)
        tv <- log(31.5)
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    fit <- .nlmixr(one.compartment, theo_sd, est = "focei",
                   control = foceiControlFast,
                   table = tableControl(cwres = FALSE))
    # the fit reports the very objective function addCwres() adds, and has no CWRES
    expect_true("FOCEi" %in% rownames(fit$objDf))
    expect_false("CWRES" %in% names(fit))

    fit2 <- suppressMessages(addCwres(fit, updateObject = FALSE))
    expect_true("CWRES" %in% names(fit2))
    # the objective function row is not duplicated (and not lost)
    expect_equal(sum(rownames(fit2$objDf) == "FOCEi"), 1L)
    expect_equal(fit2$objDf["FOCEi", "OBJF"], fit$objDf["FOCEi", "OBJF"])
  })

  test_that("addCwres() still works after setOfv(fit, 'focei')", {
    # setOfv() adds the focei objective function row WITHOUT the residual
    # columns (calcTables=FALSE), which used to leave addCwres() no way to add
    # them afterwards.  A fresh fit, not a shared fixture: setOfv() writes into
    # the fit environment in place.
    one.compartment <- function() {
      ini({
        tka <- log(1.57)
        tcl <- log(2.72)
        tv <- log(31.5)
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }
    fit <- .nlmixr(one.compartment, theo_sd, est = "saem",
                   control = saemControlFast)
    expect_false("CWRES" %in% names(fit))

    suppressMessages(setOfv(fit, "focei"))
    # now a CALCULATED focei objective function row with no CWRES
    expect_true("FOCEi" %in% rownames(fit$objDf))
    .objf <- fit$objDf["FOCEi", "OBJF"]
    expect_false(is.na(.objf))
    expect_false("CWRES" %in% names(fit))

    fit2 <- suppressMessages(addCwres(fit, updateObject = FALSE))
    expect_true("CWRES" %in% names(fit2))
    expect_equal(sum(rownames(fit2$objDf) == "FOCEi"), 1L)
    expect_equal(fit2$objDf["FOCEi", "OBJF"], .objf)
  })
})
