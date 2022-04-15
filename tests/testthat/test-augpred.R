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

})
