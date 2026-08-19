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

  test_that("saem M3/M4 residual SD matches focei (#916)", {
    # #916: the SAEM M-step's residual SSR counted a censored (M3/M4) row's
    # recorded LOQ/limit as if it had been measured, biasing the residual SD
    # low vs focei.  Data augmentation (simulate the censored DV from the
    # truncated normal implied by the current fit each E-step) fixes this --
    # pin BOTH tvK and prop.sd against focei instead of only bounding prop.sd
    # from above.
    datL <- rbind(dat[, names(dat) != "Y"], data.frame(ID = 1:10, Time = 1.5, DV = 3))
    datL$cens <- ifelse(datL$Time == 1.5, 1, 0)
    datL <- datL[order(datL$ID, datL$Time), ]

    f.foceiL <- suppressMessages(suppressWarnings(nlmixr(f, datL, "focei")))
    f.saemL <- suppressMessages(suppressWarnings(nlmixr(f, datL, "saem")))
    ct(f.saemL, "M3 censoring")
    expect_equal(as.numeric(f.saemL$theta[["tvK"]]), as.numeric(f.foceiL$theta[["tvK"]]),
                 tolerance = 0.1)
    expect_equal(as.numeric(f.saemL$theta[["prop.sd"]]), as.numeric(f.foceiL$theta[["prop.sd"]]),
                 tolerance = 0.15)

    datL4 <- datL
    datL4$limit <- 0
    f.foceiL4 <- suppressMessages(suppressWarnings(nlmixr(f, datL4, "focei")))
    f.saemL4 <- suppressMessages(suppressWarnings(nlmixr(f, datL4, "saem")))
    ct(f.saemL4, "M2 and M4 censoring")
    expect_equal(as.numeric(f.saemL4$theta[["tvK"]]), as.numeric(f.foceiL4$theta[["tvK"]]),
                 tolerance = 0.1)
    expect_equal(as.numeric(f.saemL4$theta[["prop.sd"]]), as.numeric(f.foceiL4$theta[["prop.sd"]]),
                 tolerance = 0.15)
  })

  test_that("saem mixture model with censored data fits (#916 coverage)", {
    # augmentCensY()/applyCensLoss() are also threaded through the nMix>1
    # ("parallel"/soft-EM) mixture branch (cens_mix/limit_mix), a different
    # vector-slicing path from the plain (non-mixture) case above -- fit a
    # small 2-component mixture on the same M3 data and check it converges to
    # a sane (finite, reproducible) answer instead of pinning exact values,
    # since mixture fits are prone to label-switching.
    datL <- rbind(dat[, names(dat) != "Y"], data.frame(ID = 1:10, Time = 1.5, DV = 3))
    datL$cens <- ifelse(datL$Time == 1.5, 1, 0)
    datL <- datL[order(datL$ID, datL$Time), ]

    fMix <- function() {
      ini({
        tvK1 <- 0.3
        tvK2 <- 0.8
        p1 <- 0.5
        bsvK ~ 0.04
        prop.sd <- sqrt(0.1)
      })
      model({
        ke <- mix(tvK1, p1, tvK2) * exp(bsvK)
        v <- 1
        ipre <- 10 * exp(-ke * t)
        ipre ~ prop(prop.sd)
      })
    }

    f.mix <- suppressMessages(suppressWarnings(
      nlmixr(fMix, datL, "saem",
             control = saemControl(nBurn = 10, nEm = 10, calcTables = FALSE, print = 0))
    ))
    expect_true(is.finite(f.mix$objf))
    expect_true(all(is.finite(unlist(f.mix$theta))))
    expect_true(all(as.numeric(f.mix$theta[c("tvK1", "tvK2")]) > 0.05 &
                       as.numeric(f.mix$theta[c("tvK1", "tvK2")]) < 3))
    expect_true(as.numeric(f.mix$theta[["prop.sd"]]) > 0.01 &&
                  as.numeric(f.mix$theta[["prop.sd"]]) < 2)
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
