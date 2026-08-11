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
    # focei/foce append the censored 2nd-derivative type " (laplace)"/" (gauss)" to the
    # censoring text; strip it here so these checks test the censoring METHOD (M2/M3/M4).
    expect_equal(sub(" \\((laplace|gauss)\\)$", "", as.character(model$censInformation)), censInfo)
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

  test_that("censInformation notes the censored 2nd-derivative type (laplace/gauss)", {
    fg <- suppressWarnings(suppressMessages(nlmixr(f, dat2, "posthoc")))  # gauss is the default
    fl <- suppressWarnings(suppressMessages(nlmixr(f, dat2, "posthoc", control = list(censOption = "laplace"))))
    expect_match(as.character(fg$censInformation), "\\(gauss\\)$")
    expect_match(as.character(fl$censInformation), "\\(laplace\\)$")
    expect_equal(as.character(f.focei$censInformation), "No censoring")   # no suffix when uncensored
  })




  test_that("censoring changes results - focei", {
    f.focei2 <- suppressWarnings(suppressMessages(nlmixr(f, dat2, "posthoc")))
    expect_false(isTRUE(all.equal(f.focei$objf, f.focei2$objf)))
    ct(f.focei2, "M2 censoring")
  })

  test_that("censoring changes results - saem", {
    f.saem2 <- suppressWarnings(suppressMessages(nlmixr(f, dat2, "saem")))
    ct(f.saem2, "M2 censoring")
    # censOption is inert for SAEM (no Laplace inner Hessian) -> censoring text stays PLAIN
    expect_equal(as.character(f.saem2$censInformation), "M2 censoring")
    expect_no_match(as.character(f.saem2$censInformation), "\\((laplace|gauss)\\)")
  })

  # M3 for the SAEM family (#876).  The chain scores a censored row with
  # doCensNormal1 and so does the FOCEi inner that f-SAEM's IMH acceptance uses,
  # so the two agree and fsaem needs no censoring gate -- these fits are what
  # says so.  M2 alone would not: with limit=0 the M2 tail probability is ~1, so
  # its term is ~0 and a wrong SIGN or SCALE on it does not show.  M3 has no such
  # cover, and when the chain fed doCensNormal1 the loss (not the
  # log-likelihood), the SD (not the variance), and the raw DV, these fits came
  # back at ke = -23 with a residual SD in the thousands.
  .m3 <- rbind(dat[, names(dat) != "Y"],
               data.frame(ID = 1:10, Time = 1.5, DV = 3))
  .m3$cens <- ifelse(.m3$Time == 1.5, 1, 0)
  .m3 <- .m3[order(.m3$ID, .m3$Time), ]
  .m3ctl <- saemControl(print = 0, nBurn = 200, nEm = 200, seed = 42)

  test_that("M3 censoring - saem estimates stay on the model", {
    f.saemM3 <- suppressWarnings(suppressMessages(nlmixr(f, .m3, "saem", control = .m3ctl)))
    f.foceiM3 <- suppressWarnings(suppressMessages(nlmixr(f, .m3, "focei")))
    ct(f.saemM3, "M3 censoring")
    # both fit the same model to the same data, so the elimination rate agrees to
    # well within 20%; the flipped-sign loss put it on the wrong side of zero
    expect_equal(fixef(f.saemM3)[["tvK"]], fixef(f.foceiM3)[["tvK"]], tolerance = 0.2)
    # a proportional residual SD is a fraction, not a count
    expect_lt(fixef(f.saemM3)[["prop.sd"]], 1)
    expect_lt(f.saemM3$omega[1, 1], 1)

    # M4 is its own branch (a limit turns the tail into a difference of two), so
    # it needs its own case
    .m4 <- .m3
    .m4$limit <- ifelse(.m4$cens == 1, 0, NA)
    f.saemM4 <- suppressWarnings(suppressMessages(nlmixr(f, .m4, "saem", control = .m3ctl)))
    f.foceiM4 <- suppressWarnings(suppressMessages(nlmixr(f, .m4, "focei")))
    ct(f.saemM4, "M4 censoring")
    expect_equal(fixef(f.saemM4)[["tvK"]], fixef(f.foceiM4)[["tvK"]], tolerance = 0.2)
    expect_lt(fixef(f.saemM4)[["prop.sd"]], 1)
  })

  test_that("M3 censoring - fsaem lands where saem does", {
    f.saemM3 <- suppressWarnings(suppressMessages(nlmixr(f, .m3, "saem", control = .m3ctl)))
    f.fsaemM3 <- suppressWarnings(suppressMessages(nlmixr(f, .m3, "fsaem", control = .m3ctl)))
    ct(f.fsaemM3, "M3 censoring")
    # the fast kernel really ran on the censored data (not degraded to saem)
    expect_gt(f.fsaemM3$fsaemDiag$nStep, 0)
    expect_gt(f.fsaemM3$fsaemDiag$nAcc, 0)
    # the IMH acceptance scores a censored row the way the chain does, so
    # preconditioning the chain with it must not move where the fit lands
    expect_equal(fixef(f.fsaemM3), fixef(f.saemM3), tolerance = 0.05)
    expect_equal(f.fsaemM3$omega[1, 1], f.saemM3$omega[1, 1], tolerance = 0.5)
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
