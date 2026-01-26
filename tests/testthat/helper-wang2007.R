# Helper functions for Wang2007 test suite
# These functions support the split test files for parallel execution

################################################################################
# Data Setup Functions
################################################################################

# Function to setup Wang2007 base dataset
getWang2007BaseData <- function() {
  dat <- Wang2007
  dat$DV <- dat$Y
  dat
}

# Function to setup Wang2007 dataset with dosing (dat2)
getWang2007DoseData <- function() {
  dat <- Wang2007
  dat$DV <- dat$Y
  dat2 <- dat[dat$Time == 0, ]
  dat2$EVID <- 101
  dat2$AMT <- 10
  dat2 <- rbind(dat2, data.frame(dat, EVID = 0, AMT = 0))
  dat2 <- dat2[(order(dat2$ID, -dat2$EVID, dat2$Time)), ]
  dat2
}

# Function to setup log-transformed data
getWang2007LogData <- function() {
  dat <- Wang2007
  dat$DV <- dat$Y
  datl <- dat
  datl$DV <- log(datl$DV)
  datl
}

# Function to setup log-transformed data with dosing
getWang2007LogDoseData <- function() {
  dat2 <- getWang2007DoseData()
  datl2 <- dat2
  datl2$DV <- log(datl2$DV)
  datl2
}

################################################################################
# Model Definition Functions
################################################################################

# Base PK model without ODE
getWang2007BaseModel <- function() {
  f <- function() {
    ini({
      tke <- 0.5
      eta.ke ~ 0.04
      prop.sd <- sqrt(0.1)
    })
    model({
      ke <- tke * exp(eta.ke)
      ipre <- 10 * exp(-ke * t)
      f2 <- ipre / (ipre + 5)
      f3 <- f2 * 3
      lipre <- log(ipre)
      ipre ~ prop(prop.sd)
    })
  }
  return(.nlmixr(f))
}

# Base PK model with ODE
getWang2007OdeModel <- function() {
  fo <- function() {
    ini({
      tke <- 0.5
      eta.ke ~ 0.04
      prop.sd <- sqrt(0.1)
    })
    model({
      ke <- tke * exp(eta.ke)
      d/dt(ipre) <- -ke * ipre
      f2 <- ipre / (ipre + 5)
      f3 <- f2 * 3
      lipre <- log(ipre)
      ipre ~ prop(prop.sd)
    })
  }
  return(.nlmixr(fo))
}

################################################################################
# Core Testing Function
################################################################################

testWang2007ErrorModel <- function(type, fun, val = rep(NA_real_, 10), addProp = 2, log=FALSE) {
  ## message(type)
  valName <- as.character(substitute(val))

  # Get base models
  f <- getWang2007BaseModel()
  fo <- getWang2007OdeModel()

  .f <- fun(f)
  .fo <- fun(fo)
  .dode <- getWang2007DoseData()
  .d <- getWang2007BaseData()
  if (log) {
    .d <- getWang2007LogData()
    .dode <- getWang2007LogDoseData()
  }

  fit1 <- .nlmixr(.fo, .dode, "focei",
                  control = foceiControl(
                    maxOuterIterations = 0, covMethod = "",
                    addProp = paste0("combined", addProp)))

  fit2 <- .nlmixr(.f, .d, "focei",
                  control = foceiControl(
                    maxOuterIterations = 0, covMethod = "",
                    addProp = paste0("combined", addProp)))

  fit3 <- .nlmixr(.fo, .dode, "foce",
                  control = foceiControl(
                    maxOuterIterations = 0, covMethod = "",
                    addProp = paste0("combined", addProp)))

  fit4 <- .nlmixr(.f, .d, "foce",
                  control = foceiControl(
                    maxOuterIterations = 0, covMethod = "",
                    addProp = paste0("combined", addProp)))

  fit5 <- .nlmixr(.fo, .dode, "fo",
                  control = foceiControl(
                    maxOuterIterations = 0, covMethod = "",
                    addProp = paste0("combined", addProp)))
  setOfv(fit5, "fo")

  fit6 <- .nlmixr(.f, .dode, "fo", control = foceiControl(
    maxOuterIterations = 0, covMethod = "",
    addProp = paste0("combined", addProp)))
  setOfv(fit6, "fo")

  fit7 <- .nlmixr(.fo, .dode, "agq",
                  control = agqControl(
                    maxOuterIterations = 0, covMethod = "",
                    addProp = paste0("combined", addProp)))

  fit8 <- .nlmixr(.f, .d, "agq",
                  control = agqControl(
                    maxOuterIterations = 0, covMethod = "",
                    addProp = paste0("combined", addProp)))

  fit9 <- .nlmixr(.fo, .dode, "laplace",
                  control = laplaceControl(
                    maxOuterIterations = 0, covMethod = "",
                    addProp = paste0("combined", addProp)))

  fit10 <- .nlmixr(.f, .d, "laplace",
                  control = laplaceControl(
                    maxOuterIterations = 0, covMethod = "",
                    addProp = paste0("combined", addProp)))

  .n <- paste(type, c("focei ode", "focei", "foce ode", "foce", "fo ode", "fo",
                      "agq ode", "agq", "laplace ode", "laplace"),
              paste0("combined", addProp))
  ret <- c(fit1$objective, fit2$objective, fit3$objective, fit4$objective,
           fit5$objective, fit6$objective, fit7$objective, fit8$objective,
           fit9$objective, fit10$objective)
  ret <- setNames(ret, .n)
  if (length(val) < length(.n)) {
    val <- c(val, rep(NA_real_, length(.n) - length(val)))
  }
  val <- setNames(val, .n)

  ## if (!identical(round(val, 3), round(ret, 3))) {
  ##   print("\n")
  ##   t <- try(print(str2lang(paste0(valName, " <- ", deparse1(round(setNames(ret, NULL), 3))))), silent=TRUE)
  ##   if (inherits(t, "try-error")) {
  ##     print(type)
  ##     try(print(str2lang(deparse1(round(setNames(ret, NULL), 3)))))
  ##   }
  ## }

  ## Now test
  if (!all(is.na(val))) {
    test_that(
      type,{
        expect_equal(setNames(ret, NULL), setNames(val, NULL), tolerance=1e-3)
      })
    test_that(paste0(type, " print"),{
      utils::capture.output(suppressMessages({
        withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit1), NA))
        withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit1), NA))
        withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit2), NA))
        withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit2), NA))
        withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit3), NA))
        withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit3), NA))
        withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit4), NA))
        withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit4), NA))
        withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit5), NA))
        withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit5), NA))
        withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit6), NA))
        withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit6), NA))
        withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit7), NA))
        withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit7), NA))
        withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit8), NA))
        withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit8), NA))
        withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit9), NA))
        withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit9), NA))
        withr::with_options(list(cli.unicode=TRUE),expect_error(print(fit10), NA))
        withr::with_options(list(cli.unicode=FALSE),expect_error(print(fit10), NA))
      }))
    })
  }
  invisible(ret)
}
