nmTest({
  test_that("integerish error focei", {

    PKdata <- nlmixr2data::warfarin |>
      dplyr::filter(dvid == "cp") |>
      dplyr::select(-dvid) |>
      dplyr::mutate(sex = ifelse(sex == "male", 1, 0))



      One.comp.KA.ODE <- function() {
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
          ## ODE example
          d/dt(depot)   = -ka * depot
          d/dt(central) =  ka * depot - (cl/v) * central
          ## Concentration is calculated
          cp = central/v
          ## And is assumed to follow proportional and additive error
          cp ~ prop(prop.err) + add(add.err)
        })
      }

      ## estimate parameters using nlmixr/FOCEI:

      expect_error({
        .nlmixr(
          One.comp.KA.ODE,          #the model definition
          PKdata,                   #the data set
          est = "focei",            #the estimation algorithm (FOCEi)
          #FOCEi options:
          foceiControl(print = 5))  #only print every 5th estimation step
      }, NA)

  })
})
