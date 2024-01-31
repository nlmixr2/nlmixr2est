test_that("Standard theo linCmt()", {

  one.cmt <- function() {
    ini({
      tka <- 0.45 ; label("Ka")
      tcl <- log(c(0, 2.7, 100)) ; label("Log Cl")
      tv <- 3.45; label("log V")
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

  skip_if_not(rxode2parse::.linCmtSens())

  f <- nlmixr(one.cmt)

  expect_equal(f$saemModel0,
               quote(rxModelVars({
                 ka <- exp(tka)
                 cl <- exp(tcl)
                 v <- exp(tv)
                 rx_pred_ <- linCmt()
               })))

  expect_equal(f$saemParamsToEstimate,
               c("tka", "tcl", "tv"))

  expect_equal(f$saemParams,
               "params(tka,tcl,tv)")

  expect_equal(f$saemInParsAndMuRefCovariates,
               list(inPars = character(0), covars = character(0)))

  expect_equal(f$saemFixed,
               c(tka = FALSE, tcl = FALSE, tv = FALSE))

  expect_equal(f$saemEtaTrans,
               1:3)

  expect_equal(f$saemOmegaTrans,
               1:3)

  expect_equal(f$saemModelOmega, diag(3))

  expect_equal(f$saemModelOmegaFixed,
               matrix(rep(0, 3 * 3), 3, 3))

  expect_equal(f$saemModelOmegaFixedValues,
               structure(c(0.6, 0, 0, 0, 0.3, 0, 0, 0, 0.1), .Dim = c(3L, 3L)))

  expect_equal(f$saemLow, -Inf)

  expect_equal(f$saemHi, Inf)

  expect_equal(f$saemPropT, 0) # not proportional on transformed scale

  expect_equal(f$saemYj, 2) # normal translation

  expect_equal(f$saemResMod, 1) # additive

  # keep all omegas
  expect_equal(f$saemParHistOmegaKeep,
               c(eta.ka = 1L, eta.cl = 1L, eta.v = 1L))

  expect_equal(f$saemParHistNames,
               c("tka", "tcl", "tv", "V(eta.ka)", "V(eta.cl)", "V(eta.v)", "add.sd"))

  # y= f + (a + b*f^c)*varepsilon # combined1
  # y= f + sqrt(a^2 + b^2*f^(2*c))*varepsilon # combined2

  # Ares
  expect_equal(f$saemAres, 0.7)

  # Bres
  expect_equal(f$saemBres, 1)

  # Cres
  expect_equal(f$saemCres, 1)

  # lambda for boxCox or yeoJohnson
  expect_equal(f$saemLres, 1)

  expect_equal(f$saemLogEta, c(tka = TRUE, tcl = TRUE, tv = TRUE))

  expect_equal(f$saemInitTheta,
               structure(c(1.56831218549017, 2.7, 31.5003923087479), .Names = c("", "", "")))

  expect_equal(f$saemInitOmega,
               c(tka = 0.6, tcl = 0.3, tv = 0.1))

  expect_equal(f$saemThetaDataFrame,
               structure(list(lower = c(-Inf, -Inf, -Inf, -Inf), theta = c(0.45, 0.993251773010283, 3.45, 0.7), fixed = c(FALSE, FALSE, FALSE, FALSE), upper = c(Inf, Inf, Inf, Inf)), class = "data.frame", row.names = c("tka", "tcl", "tv", "add.sd")))

  expect_equal(f$saemModelPredReplaceLst,
               c(tka = "THETA[1] + ETA[1]", tcl = "THETA[2] + ETA[2]", tv = "THETA[3] + ETA[3]", add.sd = "THETA[4]"))

})

test_that("non mu-ref theo linCmt() with fixed components", {

  one.cmt <- function() {
    ini({
      tka <- exp(0.45) ; label("Ka")
      tcl <- log(c(0, 2.7, 100)) ; label("Log Cl")
      tv <- fix(3.45); label("log V")
      eta.ka ~ 0.6
      eta.cl ~ fix(0.3)
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- tka * exp(eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  f <- nlmixr(one.cmt)

  expect_equal(f$saemModel0,
               quote(rxModelVars({
                 ka <- tka * exp(eta.ka)
                 cl <- exp(tcl)
                 v <- exp(tv)
                 rx_pred_ <- linCmt()
               })))

  expect_equal(f$saemParamsToEstimate,
               c("tka", "tcl", "tv", "eta.ka"))

  expect_equal(f$saemParams,
               "params(tka,tcl,tv,eta.ka)")

  expect_equal(f$saemInParsAndMuRefCovariates,
               list(inPars = character(0), covars = character(0)))

  expect_equal(f$saemFixed,
               c(tka = FALSE, tcl = FALSE, tv = TRUE, eta.ka=TRUE))

  expect_equal(f$saemEtaTrans,
               c(4L, 2L, 3L))

  expect_equal(f$saemOmegaTrans,
               c(3L, 1L, 2L))

  expect_equal(f$saemModelOmega,
               structure(c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1), .Dim = c(4L, 4L)))

  expect_equal(f$saemModelOmegaFixed,
               structure(c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(4L, 4L)))

  expect_equal(f$saemModelOmegaFixedValues,
               structure(c(0, 0, 0, 0, 0, 0.3, 0, 0, 0, 0, 0.1, 0, 0, 0, 0, 0.6), .Dim = c(4L, 4L)))

  expect_equal(f$saemLow, -Inf)

  expect_equal(f$saemHi, Inf)

  expect_equal(f$saemPropT, 0) # not proportional on transformed scale

  expect_equal(f$saemYj, 2) # normal translation

  expect_equal(f$saemResMod, 1) # additive

  # Don't keep the eta.cl since it is fixed
  expect_equal(f$saemParHistOmegaKeep,
               c(eta.cl = 0L, eta.v = 1L, eta.ka = 1L))

  # fixed volume and fixed eta.cl
  expect_equal(f$saemParHistNames,
               c("tka", "tcl", "V(eta.v)", "V(eta.ka)", "add.sd"))

  # y= f + (a + b*f^c)*varepsilon # combined1
  # y= f + sqrt(a^2 + b^2*f^(2*c))*varepsilon # combined2

  # Ares
  expect_equal(f$saemAres, 0.7)

  # Bres
  expect_equal(f$saemBres, 1)

  # Cres
  expect_equal(f$saemCres, 1)

  # lambda for boxCox or yeoJohnson
  expect_equal(f$saemLres, 1)

  expect_equal(f$saemLogEta, c(tka = FALSE, tcl = TRUE, tv = TRUE, eta.ka=TRUE))

  expect_equal(f$saemInitTheta,
               c(1.56831218549017, 2.7, FIXED=31.5003923087479, FIXED=1))

  expect_equal(f$saemInitOmega,
               c(tka = 1.0, tcl = 0.3, tv = 0.1, eta.ka=0.6))

  expect_equal(f$saemThetaDataFrame,
               structure(list(lower = c(-Inf, -Inf, -Inf, -Inf),
                              theta = c(exp(0.45), 0.993251773010283, 3.45, 0.7),
                              fixed = c(FALSE, FALSE, TRUE, FALSE),
                              upper = c(Inf, Inf, Inf, Inf)),
                         class = "data.frame",
                         row.names = c("tka", "tcl", "tv", "add.sd")))

  expect_equal(f$saemModelPredReplaceLst,
               c(tka = "THETA[1]", tcl = "THETA[2] + ETA[2]", tv = "THETA[3] + ETA[3]", eta.ka = "ETA[1]", add.sd = "THETA[4]"))

})

test_that("theo wt cov parsing", {

  one.cmt <- function() {
    ini({
      tka <- 0.45 ; label("Ka")
      tcl <- log(c(0, 2.7, 100)) ; label("Log Cl")
      tv <- 3.45; label("log V")
      cl.wt <- 0
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + wt * cl.wt)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }

  f <- one.cmt()

  expect_equal(f$muRefCovariateDataFrame,
               structure(list(theta = "tcl", covariate = "wt", covariateParameter = "cl.wt"), row.names = c(NA, -1L), class = "data.frame"))

  expect_equal(f$saemMuRefCovariateDataFrame,
               structure(list(theta = "tcl", covariate = "wt", covariateParameter = "cl.wt"), row.names = c(NA, -1L), class = "data.frame"))

  expect_equal(f$saemModelPredReplaceLst,
               c(tka = "THETA[1] + ETA[1]",
                 tcl = "THETA[2] + ETA[2] + wt * THETA[4]",
                 tv = "THETA[3] + ETA[3]",
                 cl.wt = "THETA[4]",
                 add.sd = "THETA[5]"))

})

test_that("nimo parsing", {

  nimo <- function() {
    ini({
      tcl <- log(0.001)
      tv1 <- log(1.45)
      tQ <- log(0.004)
      tv2 <- log(44)
      tkss <- log(12)
      tkint <- log(0.3)
      tksyn <- log(1)
      tkdeg <- log(7)
      eta.cl  ~ 2
      eta.v1  ~ 2
      eta.kss ~ 2
      add.err <- 10
    })
    model({
      cl <- exp(tcl + eta.cl)
      v1 <- exp(tv1 + eta.v1)
      Q  <- exp(tQ)
      v2 <- exp(tv2)
      kss <- exp(tkss + eta.kss)
      kint <- exp(tkint)
      ksyn <- exp(tksyn)
      kdeg <- exp(tkdeg)

      k <- cl/v1
      k12 <- Q/v1
      k21 <- Q/v2

      eff(0) <- ksyn/kdeg ##initializing compartment

      ## Concentration is calculated
      conc = 0.5*(central/v1-eff-kss)+0.5*sqrt((central/v1-eff-kss)**2+4*kss*central/v1)

      d/dt(central)  = -(k+k12)*conc*v1+k21*peripheral-kint*eff*conc*v1/(kss+conc)
      d/dt(peripheral) = k12*conc*v1-k21*peripheral  ##Free Drug second compartment amount
      d/dt(eff) = ksyn - kdeg*eff - (kint-kdeg)*conc*eff/(kss+conc)

      IPRED=log(conc)

      IPRED ~ add(add.err)
    })
  }

  f <- nlmixr(nimo)

  expect_equal(
    f$saemModel0,
    quote(rxModelVars({
      cl <- exp(tcl)
      v1 <- exp(tv1)
      Q <- exp(tQ)
      v2 <- exp(tv2)
      kss <- exp(tkss)
      kint <- exp(tkint)
      ksyn <- exp(tksyn)
      kdeg <- exp(tkdeg)
      k <- cl/v1
      k12 <- Q/v1
      k21 <- Q/v2
      eff(0) <- ksyn/kdeg
      conc = 0.5 * (central/v1 - eff - kss) + 0.5 * sqrt((central/v1 -
                                                            eff - kss)^2 + 4 * kss * central/v1)
      d/dt(central) = -(k + k12) * conc * v1 + k21 * peripheral -
        kint * eff * conc * v1/(kss + conc)
      d/dt(peripheral) = k12 * conc * v1 - k21 * peripheral
      d/dt(eff) = ksyn - kdeg * eff - (kint - kdeg) * conc * eff/(kss + conc)
      IPRED = log(conc)
      rx_pred_ <- IPRED
    })))

  expect_equal(
    f$saemParamsToEstimate,
    c("tcl", "tv1", "tQ", "tv2", "tkss", "tkint", "tksyn", "tkdeg")
  )

  expect_equal(
    f$saemParams,
    "params(tcl,tv1,tQ,tv2,tkss,tkint,tksyn,tkdeg)"
  )

  expect_equal(
    f$saemInParsAndMuRefCovariates,
    list(inPars = character(0), covars = character(0))
  )

  expect_equal(
    f$saemFixed,
    c(tcl = FALSE, tv1 = FALSE, tQ = FALSE, tv2 = FALSE, tkss = FALSE, tkint = FALSE, tksyn = FALSE, tkdeg = FALSE)
  )

  # eta.cl, eta.v1, eta.kss
  expect_equal(
    f$saemEtaTrans,
    c(1L, 2L, 5L)
  )

  expect_equal(f$saemOmegaTrans, 1:3)

  expect_equal(
    f$saemModelOmegaFixedValues,
    structure(c(2, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(8L, 8L))
  )

  expect_equal(f$saemLow, -Inf)

  expect_equal(f$saemHi, Inf)

  expect_equal(f$saemPropT, 0) # not proportional on transformed scale

  expect_equal(f$saemYj, 2) # normal translation

  expect_equal(f$saemResMod, 1) # additive

  # Don't keep the eta.cl since it is fixed
  expect_equal(f$saemParHistOmegaKeep,
               c(eta.cl = 1L, eta.v1 = 1L, eta.kss = 1L))

  # fixed volume and fixed eta.cl
  expect_equal(f$saemParHistNames,
               c("tcl", "tv1", "tQ", "tv2", "tkss", "tkint", "tksyn", "tkdeg",
                 "V(eta.cl)", "V(eta.v1)", "V(eta.kss)", "add.err"))

  # y= f + (a + b*f^c)*varepsilon # combined1
  # y= f + sqrt(a^2 + b^2*f^(2*c))*varepsilon # combined2

  # Ares
  expect_equal(f$saemAres, 10)

  # Bres
  expect_equal(f$saemBres, 1)

  # Cres
  expect_equal(f$saemCres, 1)

  # lambda for boxCox or yeoJohnson
  expect_equal(f$saemLres, 1)

  expect_equal(
    f$saemLogEta,
    c(tcl = TRUE, tv1 = TRUE, tQ = TRUE, tv2 = TRUE, tkss = TRUE, tkint = TRUE, tksyn = TRUE, tkdeg = TRUE)
  )

  expect_equal(
    f$saemInitTheta,
    structure(c(0.001, 1.45, 0.004, 44, 12, 0.3, 1, 7), .Names = c("", "", "", "", "", "", "", ""))
  )

  expect_equal(
    f$saemInitOmega,
    c(tcl = 2, tv1 = 2, tQ = 1, tv2 = 1, tkss = 2, tkint = 1, tksyn = 1, tkdeg = 1)
  )

  expect_equal(
    f$saemThetaDataFrame,
    structure(
      list(
        lower = c(-Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf, -Inf),
        theta = c(-6.90775527898214, 0.371563556432483, -5.52146091786225, 3.78418963391826, 2.484906649788, -1.20397280432594, 0, 1.94591014905531, 10),
        fixed = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
        upper = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf)),
      class = "data.frame",
      row.names = c("tcl", "tv1", "tQ", "tv2", "tkss", "tkint", "tksyn", "tkdeg", "add.err"))
  )

})
