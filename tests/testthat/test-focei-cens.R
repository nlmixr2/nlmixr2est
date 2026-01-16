nmTest({
  # M2 all observations

  dat <- Wang2007
  dat$DV <- dat$Y # Add the required DV data item


  f <- function() {
    ini({
      tvK <- 0.5 # Typical Value of K
      bsvK ~ 0.04 # Between Subject Variance of K
      prop.sd <- sqrt(0.1)
    })
    model({
      ke <- tvK * exp(bsvK)
      v <- 1
      ipre <- 10 * exp(-ke * t)
      ipre ~ prop(prop.sd)
    })
  }

  ct <- function(model, censInfo) {
    expect_equal(as.character(model$censInformation), censInfo)
  }

  dat2 <- dat
  dat2$limit <- 0

  dat3 <- dat
  dat3$limit <- 3

  dat4 <- dat
  dat4$limit <- 12

  f.foce <- suppressMessages(suppressWarnings(nlmixr(f, dat, "posthoc", control = list(interaction = FALSE))))
  f.focei <- suppressWarnings(suppressMessages(nlmixr(f, dat, "posthoc")))

  test_that("censoring information is correct", {
    ct(f.foce, "No censoring")
    ct(f.focei, "No censoring")
  })




  test_that("censoring changes results - focei", {
    f.focei2 <- suppressWarnings(suppressMessages(nlmixr(f, dat2, "posthoc")))
    expect_false(isTRUE(all.equal(f.focei$objf, f.focei2$objf)))
    ct(f.focei2, "M2 censoring")
  })

  test_that("censoring changes results - saem", {
    skip_on_os("linux") # SAEM on linux CI seems unstable for this test, though runs fine locally
    f.saem2 <- suppressWarnings(suppressMessages(nlmixr(f, dat2, "saem")))
    ct(f.saem2, "M2 censoring")
  })


  test_that("Limit affects values", {

    f.focei3 <- suppressMessages(suppressWarnings(nlmixr(f, dat3, "posthoc")))
    ct(f.focei3, "M2 censoring")


    f.focei4 <- suppressMessages(suppressWarnings(nlmixr(f, dat4, "posthoc")))

    ct(f.focei4, "M2 censoring")

    f.foce2 <- suppressMessages(suppressWarnings(nlmixr(f, dat2, "posthoc", control = list(interaction = FALSE))))
    expect_false(isTRUE(all.equal(f.foce$objf, f.foce2$objf)))

    ct(f.foce2, "M2 censoring")

    f.foce3 <- suppressMessages(suppressWarnings(nlmixr(f, dat3, "posthoc", control = list(interaction = FALSE))))
    expect_false(isTRUE(all.equal(f.foce2$objf, f.foce3$objf)))

    ct(f.foce3, "M2 censoring")

  })

  test_that("M3/M4 -- Missing, assume LLOQ=3 at t=1.5", {
    skip_if(Sys.getenv("R_ARCH") == "/i386", "windows32")

    datL <- rbind(dat[, names(dat) != "Y"], data.frame(ID = 1:10, Time = 1.5, DV = 3))
    datL$cens <- ifelse(datL$Time == 1.5, 1, 0)
    datL <- datL[order(datL$ID, datL$Time), ]

    datL4 <- datL
    datL4$limit <- 0

    f.foceiL <- suppressMessages(suppressWarnings(nlmixr(f, datL, "posthoc")))
    expect_false(isTRUE(all.equal(f.focei$objf, f.foceiL$objf)))
    ct(f.foceiL, "M3 censoring")

    datL2o3 <- datL
    datL2o3$limit <- NA
    datL2o3$limit[1] <- 0

    assign("curdat", datL2o3, env=globalenv())
    assign("f", f, env=globalenv())

    f.foceiL2o3 <- suppressMessages(suppressWarnings(nlmixr(f, datL2o3, "posthoc")))
    ct(f.foceiL2o3, "M2 and M3 censoring")

    f.foceiL <- suppressMessages(suppressWarnings(nlmixr(f, datL, "posthoc")))
    expect_false(isTRUE(all.equal(f.focei$objf, f.foceiL$objf)))
    ct(f.foceiL, "M3 censoring")


    f.foceiL4 <- suppressMessages(suppressWarnings(nlmixr(f, datL4, "posthoc")))
    expect_false(isTRUE(all.equal(f.focei$objf, f.foceiL4$objf)))
    expect_false(isTRUE(all.equal(f.foceiL$objf, f.foceiL4$objf)))
    ct(f.foceiL4, "M2 and M4 censoring")

    datL4only <- datL4
    datL4only$limit <- ifelse(datL4only$cens == 0, NA, 0)

    f.foceiL4only <- suppressMessages(suppressWarnings(nlmixr(f, datL4only, "posthoc")))
    ct(f.foceiL4only, "M4 censoring")

    w <- which(datL4only$cens == 1)
    datL4and3 <- datL4only
    datL4and3$limit[w[1]] <- NA

    f.foceiL4and3 <- suppressMessages(suppressWarnings(nlmixr(f, datL4and3, "posthoc")))
    ct(f.foceiL4and3, "M3 and M4 censoring")

    datL4o3o2 <- datL4and3
    datL4o3o2$limit[1] <- 0

    f.foceiL4o3o2 <- suppressMessages(suppressWarnings(nlmixr(f, datL4o3o2, "posthoc")))
    ct(f.foceiL4o3o2, "M2, M3 and M4 censoring")

    datL <- rbind(dat[, names(dat) != "Y"], data.frame(ID = 1:10, Time = 1.5, DV = 3))
    datL$cens <- ifelse(datL$Time == 1.5, 1, 0)
    datL <- datL[order(datL$ID, datL$Time), ]

    datL4 <- datL
    datL4$limit <- 0

    f.foceiL <- suppressMessages(suppressWarnings(nlmixr(f, datL, "posthoc")))
    expect_false(isTRUE(all.equal(f.focei$objf, f.foceiL$objf)))
    ct(f.foceiL, "M3 censoring")

    f.foceiL4 <- suppressMessages(suppressWarnings(nlmixr(f, datL4, "posthoc")))
    expect_false(isTRUE(all.equal(f.focei$objf, f.foceiL4$objf)))
    expect_false(isTRUE(all.equal(f.foceiL$objf, f.foceiL4$objf)))
    ct(f.foceiL, "M3 censoring")

    ## foce

    datL <- rbind(dat[, names(dat) != "Y"], data.frame(ID = 1:10, Time = 1.5, DV = 3))
    datL$cens <- ifelse(datL$Time == 1.5, 1, 0)
    datL <- datL[order(datL$ID, datL$Time), ]

    datL4 <- datL
    datL4$limit <- 0

    f.foceL <- suppressMessages(suppressWarnings(nlmixr(f, datL, "posthoc", control = list(interaction = FALSE))))
    expect_false(isTRUE(all.equal(f.foce$objf, f.foceL$objf)))
    ct(f.foceiL4, "M2 and M4 censoring")


    f.foceL4 <- suppressMessages(suppressWarnings(nlmixr(f, datL4, "posthoc", control = list(interaction = FALSE))))
    expect_false(isTRUE(all.equal(f.foce$objf, f.foceL4$objf)))
    expect_false(isTRUE(all.equal(f.foceL$objf, f.foceL4$objf)))
    ct(f.foceiL4, "M2 and M4 censoring")

    upperDat <- datL4
    names(upperDat)[4] <- "CENS"
    names(upperDat)[5] <- "LIMIT"

    f.foceL4u <- suppressMessages(suppressWarnings(nlmixr(f, datL4, "posthoc", control = list(interaction = FALSE))))

    expect_equal(names(f.foceL4u), names(f.foceL4))

  })

})
