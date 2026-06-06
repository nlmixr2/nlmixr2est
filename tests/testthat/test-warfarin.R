nmTest({
  test_that("warfarin residuals do not go to zero", {

    PKdata <- nlmixr2data::warfarin
    PKdata <- PKdata[PKdata$dvid == "cp", ]

    One.comp.transit.allo <- function() {
      ini({
        # Where initial conditions/variables are specified
        lktr <- log(1.15)  #log k transit (/h)
        lcl  <- log(0.15)  #log Cl (L/hr)
        lv   <- log(7)     #log V (L)
        ALLC <- fix(0.75)  #allometric exponent cl
        ALLV <- fix(1.00)  #allometric exponent v
        prop.err <- 0.15   #proportional error (SD/mean)
        add.err <- 0.6     #additive error (mg/L)
        eta.ktr ~ 0.5
        eta.cl ~ 0.1
        eta.v ~ 0.1
      })
      model({
        #Allometric scaling on weight
        cl <- exp(lcl + eta.cl + ALLC * log(wt / 70))
        v  <- exp(lv + eta.v + ALLV * log(wt / 70))
        ktr <- exp(lktr + eta.ktr)
        # RxODE-style differential equations are supported
        d/dt(depot)   = -ktr * depot
        d/dt(central) =  ktr * trans - (cl/v) * central
        d/dt(trans)   =  ktr * depot - ktr * trans
        ## Concentration is calculated
        cp = central/v
        # And is assumed to follow proportional and additive error
        cp ~ prop(prop.err) + add(add.err)
      })
    }

    run007F <-
      suppressMessages(nlmixr(One.comp.transit.allo,
                              PKdata,
                              est = "focei",
                              foceiControl(print = 0)))

    expect_true(all(run007F$theta[c("prop.err", "add.err")] > 0.0001))

    One.comp.transit.allo <- function() {
      ini({
        # Where initial conditions/variables are specified
        lktr <- log(1.15)  #log k transit (/h)
        lcl  <- log(0.15)  #log Cl (L/hr)
        lv   <- log(7)     #log V (L)
        prop.err <- 0.15   #proportional error (SD/mean)
        add.err <- 0.6     #additive error (mg/L)
        eta.ktr ~ 0.5
        eta.cl ~ 0.1
        eta.v ~ 0.1
      })
      model({
        #Allometric scaling on weight
        cl <- exp(lcl + eta.cl + 0.75 * log(wt / 70))
        v  <- exp(lv + eta.v + 1.00 * log(wt / 70))
        ktr <- exp(lktr + eta.ktr)
        # RxODE-style differential equations are supported
        d/dt(depot)   = -ktr * depot
        d/dt(central) =  ktr * trans - (cl/v) * central
        d/dt(trans)   =  ktr * depot - ktr * trans
        ## Concentration is calculated
        cp = central/v
        # And is assumed to follow proportional and additive error
        cp ~ prop(prop.err) + add(add.err)
      })
    }

    run007F2 <-
      suppressMessages(nlmixr(One.comp.transit.allo,
                              PKdata,
                              est = "focei",
                              foceiControl(print = 0)))

    expect_equal(run007F2$objf, run007F$objf)


  })
})
