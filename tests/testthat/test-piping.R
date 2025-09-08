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

    mod1 <- mod %>% rmEta("etaKa")

    expect_true(all(is.na(mod1$iniDf$neta1)))

    mod1 <- mod %>% rmEta(etaKa)

    expect_true(all(is.na(mod1$iniDf$neta1)))

  })
})
