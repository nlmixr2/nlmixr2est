skip_on_cran()
test_that("test augPred", {

  PKdata <- nlmixr2data::warfarin %>%
    dplyr::filter(dvid == "cp") %>%
    dplyr::select(-dvid) %>%
    dplyr::mutate(sex = ifelse(sex == "male", 1, 0))

  One.comp.KA.solved <- function() {
      ini({
        # Where initial conditions/variables are specified
        lka  <- log(1.15)  #log ka (1/h)
        lcl  <- log(0.135) #log Cl (L/h)
        lv   <- log(8)     #log V (L)
        prop.err <- 0.15   #proportional error (SD/mean)
        add.err  <- 0.6    #additive error (mg/L)
        eta.ka ~ 0.5   #IIV ka
        eta.cl ~ 0.1   #IIV cl
        eta.v  ~ 0.1   #IIV v
      })
      model({
        # Where the model is specified
        cl <- exp(lcl + eta.cl)
        v  <- exp(lv + eta.v)
        ka <- exp(lka + eta.ka)
        linCmt() ~ prop(prop.err) + add(add.err)
      })
    }


  fitOne.comp.KA.solved_S <-
    nlmixr(
      One.comp.KA.solved,    #the model definition
      PKdata,                #the data set
      est = "saem",          #the estimation algorithm (SAEM)
      #the SAEM minimisation options:
      saemControl(nBurn = 200, #200 SAEM burn-in iterations (the default)
                  nEm   = 300, #300 EM iterations (the default)
                  print = 50), #print every 50th iteration
      #only print every 50th estimation step (default=1 which gives endless output)
      tableControl(cwres = TRUE,npde=TRUE) #calculates NONMEM-style conditional weighted residuals and npde for diagnostics
    )

  expect_error(augPred(fitOne.comp.KA.solved_S), NA)

  # Data' purpose illustrates the error and my data set
  df <- tibble::tibble(
    ID = c(rep(1, 6), rep(2, 6)),
    TIME = c(0.00, 12.11, 18.41, 23.89, 36.00, 43.51, 0.00, 12.00, 20.00, 24.00, 36.80, 45.00),
    AMT = c(1000, 1000, NA, 1000, 1000, NA, 1000, 2000, NA, 1000, 1000, NA),
    DUR = c(2.5, 2.5, NA, 2.5, 2.5, NA, 2.5, 2.5, NA, 2.5, 2.5, NA),
    DV = c(NA, NA, 3.0, NA, NA, 9.6, NA, NA, 7.0, NA, NA, 2.8),
    WT = c(rep(55, 6), rep(48, 6))
  ) %>%
    dplyr::mutate(EVID = ifelse(is.na(DV), 1, 0))

  fun <- function() {
    ini({
      tvCl <- c(0, 4, Inf)
      tvVc <- c(0, 48, Inf)

      eta.Vc ~ 0.62
      prop.sd <- 0.051529

    })
    model({
      Cl <- tvCl
      Vc <- tvVc*(WT/70)*exp(eta.Vc)

      # dynamical system
      linCmt() ~ prop(prop.sd)
    })
  }

  fit <- nlmixr2(fun, df, list(print=0), est="posthoc")

  expect_error(augPred(fit), NA)

})
