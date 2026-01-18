nmTest({
  test_that("rmEta()", {

    mod <- function ()  {
      description <- "One compartment PK model with linear clearance"
      ini({
        lka <- 0.45
        lcl <- 1
        lvc <- 3.45
        propSd <- c(0, 0.5)
        etaKa ~ 0.1
      })
      model({
        ka <- exp(lka + etaKa)
        cl <- exp(lcl)
        vc <- exp(lvc)
        Cc <- linCmt()
        Cc ~ prop(propSd)
      })
    }

    mod1 <- mod |>
      rmEta("etaKa")

    expect_true(all(is.na(mod1$iniDf$neta1)))

    mod1 <- mod |>
      rmEta(etaKa)

    expect_true(all(is.na(mod1$iniDf$neta1)))

    mod <- function ()  {
      description <- "One compartment PK model with linear clearance"
      ini({
        lka <- 0.45
        lcl <- 1
        lvc <- 3.45
        propSd <- c(0, 0.5)
        etaKa ~ 0.1
        etaCl ~ 0.2
        etaVc ~ 0.3
      })
      model({
        ka <- exp(lka + etaKa)
        cl <- exp(lcl + etaCl)
        vc <- exp(lvc + etaVc)
        Cc <- linCmt()
        Cc ~ prop(propSd)
      })
    }

    mod1 <- mod |>
      rmEta(c("etaKa", "etaCl"))

    expect_true(rxode2::testExists(mod1, "etaVc"))
    expect_true(!rxode2::testExists(mod1, "etaKa"))
    expect_true(!rxode2::testExists(mod1, "etaCl"))

    expect_error(mod |>
                   rmEta("etaX"), regexp = "'etaX' not in the model")


    mod <- function ()  {
      description <- "One compartment PK model with linear clearance"
      ini({
        lka <- 0.45
        lcl <- 1
        lvc <- 3.45
        propSd <- c(0, 0.5)
        etaKa ~ 0.1
        etaCl ~ 0.2
        etaVc ~ 0.3
      })
      model({
        ka <- exp(lka + etaKa)
        cl <- exp(lcl + etaCl)
        vc <- lvc*exp(etaVc)
        Cc <- linCmt()
        Cc ~ prop(propSd)
      })
    }

    mod1 <- mod |>
      rmEta("etaVc")

    expect_true(!rxode2::testExists(mod1, "etaVc"))
    expect_true(rxode2::testExists(mod1, "etaKa"))
    expect_true(rxode2::testExists(mod1, "etaCl"))



  })
})
