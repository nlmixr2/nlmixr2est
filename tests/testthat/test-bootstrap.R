skip_on_cran()

test_that("sampling should return different datasets at each call", {
  a <- digest::digest(nlmixr2est:::sampling(theo_sd))
  b <- digest::digest(nlmixr2est:::sampling(theo_sd))
  expect_false(isTRUE(all.equal(a, b)))
})

test_that("resuming the fit should not return the same datasets as before", {

  one.cmt <- function() {
    ini({
      tka <- 0.45 ; label("Log Ka")
      tcl <- 1 ; label("Log Cl")
      tv <- 3.45 ; label("log V")
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  fit <- suppressMessages(suppressWarnings(nlmixr(
    one.cmt,
    theo_sd,
    est = "focei",
    control = list(print = 0),
    table = list(npde = TRUE, cwres = TRUE))))

  fit1 <- suppressMessages(nlmixr2est:::bootstrapFit(fit, nboot = 2, restart = TRUE))
  fit2 <- suppressMessages(nlmixr2est:::bootstrapFit(fit, nboot = 4, restart = FALSE))

  output_dir <-
    paste0("nlmixr2BootstrapCache_", "fit", "_", fit$bootstrapMd5)

  fnameBootDataPattern <- paste0("boot_data",
                                 "_", "[0-9]+", ".rds",
                                 sep = "")

  files <- list.files(paste0("./", output_dir), pattern = fnameBootDataPattern, full.names=TRUE)

  ## print(output_dir)

  fitdata <- lapply(files, function(x) {
    readRDS(x)
  })

  a <- digest::digest(fitdata[[1]])
  b <- digest::digest(fitdata[[3]])
  expect_false(isTRUE(all.equal(a, b)))

  a <- digest::digest(fitdata[[2]])
  b <- digest::digest(fitdata[[4]])
  expect_false(isTRUE(all.equal(a, b)))

  unlink(output_dir, recursive = TRUE, force = TRUE)
})

test_that("different confidence levels should result in different bands", {

  one.cmt <- function() {
    ini({
      tka <- 0.45 ; label("Log Ka")
      tcl <- 1 ; label("Log Cl")
      tv <- 3.45 ; label("log V")
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  fit <- suppressMessages(suppressWarnings(nlmixr2(
    one.cmt,
    theo_sd,
    est = "focei",
    control = list(print = 0),
    table = list(npde = TRUE, cwres = TRUE))))

  fitlist <- suppressMessages(nlmixr2est:::modelBootstrap(fit, nboot = 4, restart = TRUE)[[1]])
  bootSummary1 <- nlmixr2est:::getBootstrapSummary(fitlist, ci = 0.95)
  bootSummary2 <- nlmixr2est:::getBootstrapSummary(fitlist, ci = 0.75)

  a <- digest::digest(bootSummary1$parFixedDf$confLower)
  b <- digest::digest(bootSummary2$parFixedDf$confLower)

  expect_false(isTRUE(all.equal(a, b)))

  a <- digest::digest(bootSummary1$parFixedDf$confUpper)
  b <- digest::digest(bootSummary2$parFixedDf$confUpper)
  expect_false(isTRUE(all.equal(a, b)))

})

test_that("expected columns in fit$parFixedDf object should match", {

  one.cmt <- function() {
    ini({
      tka <- 0.45 ; label("Log Ka")
      tcl <- 1 ; label("Log Cl")
      tv <- 3.45 ; label("log V")
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  fit <- suppressMessages(suppressWarnings(nlmixr2(
    one.cmt,
    theo_sd,
    est = "focei",
    control = list(print = 0),
    table = list(npde = TRUE, cwres = TRUE))))

  colsBefore <- colnames(fit$parFixedDf)
  fitlist <- suppressMessages(nlmixr2est:::modelBootstrap(fit, nboot = 4, restart = TRUE)[[1]])

  bootSummary <- nlmixr2est:::getBootstrapSummary(fitlist, ci = 0.95)

  colsAfter <- colnames(fit$parFixedDf)

  expect_equal(colsAfter, colsBefore)

  lapply(
    list.files("./", pattern = "nlmixr2BootstrapCache_.*"),
    function(x) {
      unlink(x, recursive = TRUE, force = TRUE)
    })
})


test_that("saem bootstrap", {

  one.cmt <- function() {
    ini({
      tka <- 0.45 ; label("Log Ka")
      tcl <- 1 ; label("Log Cl")
      tv <- 3.45 ; label("log V")
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  fit <- suppressMessages(suppressWarnings(nlmixr(
    one.cmt,
    theo_sd,
    est = "saem",
    control = list(print = 0),
    table = list(npde = TRUE, cwres = TRUE))))

  expect_error(fit1 <- suppressMessages(nlmixr2est:::bootstrapFit(fit, nboot = 2, restart = TRUE)), NA)

  output_dir <-
    paste0("nlmixr2BootstrapCache_", "fit", "_", fit$bootstrapMd5)


})
