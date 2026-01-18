nmTest({
  test_that("test lag with warfarin", {
    KA1Lode <- function() {
      ini({
        # Where initial conditions/variables are specified
        ltlag <- log(0.2)  #log Tlag (h)
        lka  <- log(1.15)  #log ka (/h)
        lcl  <- log(0.135) #log Cl (L/h)
        lv   <- log(8)     #log V (L)
        prop.err <- 0.15   #proportional error (SD/mean)
        add.err  <- 0.6    #additive error (mg/L)
        eta.tlag ~ 0.5 #IIV tlag
        eta.ka ~ 0.5   #IIV ka
        eta.cl ~ 0.1   #IIV cl
        eta.v  ~ 0.1   #IIV v
      })
      model({
        # Where the model is specified
        tlag <- exp(ltlag + eta.tlag)
        ka <- exp(lka + eta.ka)
        cl <- exp(lcl + eta.cl)
        v  <- exp(lv + eta.v)
        ## ODE example
        d/dt(gut)= - ka*gut
        d/dt(central)= ka*gut - (cl/v)*central
        lag(gut)=tlag
        cp=central/v
        ## where residual error is assumed to follow proportional and additive error
        cp ~ prop(prop.err) + add(add.err)
      })
    }

    d <- nlmixr2data::warfarin |>
      dplyr::filter(dvid=="cp")

    f <- .nlmixr(KA1Lode, d, "focei")

    expect_true(f$objf < 500)
  })
})
