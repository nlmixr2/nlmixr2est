nmTest({
  test_that("omega boundary", {

    one.compartment.IV.MM.model <- function(){
      ini({
        lVM <- 7
        label("log Vmax (mg/hr)")
        lKM <- 5.7
        label("log KM (mg/L)")
        lVc <- 4.5
        label("log Vc (L)")
        prop.err <- 0.3
        label("Proportional error")
        eta.Vc ~ 0.15
        label("IIV Vc")
        eta.VM ~ 0.15
        label("IIV Vmax")
        eta.KM ~ 0.15
        label("IIV KM")
      })
      model({
        Vc <- exp(lVc + eta.Vc)
        VM <- exp(lVM + eta.VM)
        KM <- exp(lKM + eta.KM)
        d/dt(centr)  <- -(VM*centr/Vc)/(KM+centr/Vc)
        cp <- centr / Vc
        cp ~ prop(prop.err)
      })
    }


    datr <- nlmixr2data::Infusion_1CPTMM

    dataX<-datr[datr$SD==0, ]

    .nlmixr <- function(...) {
      suppressWarnings(suppressMessages(nlmixr2(...)))
    }

    expect_true(inherits(.nlmixr(one.compartment.IV.MM.model,dataX,est="focei",
                                 control=list(print=0)), "nlmixr2FitCore"))

  })
})
