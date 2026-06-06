test_that("est='simulation'", {

  f <- function(){
    ini({
      lCl <- 1.6      #log Cl (L/hr)
      lVc <- log(90)  #log Vc (L)
      lKA <- 0.1      #log Ka (1/hr)
      prop.err <- c(0, 0.2, 1)
      eta.Cl ~ 0.1   # BSV Cl
      eta.Vc ~ 0.1   # BSV Vc
      eta.KA ~ 0.1   # BSV Ka
    })
    model({
      Cl <- exp(lCl + eta.Cl)
      Vc = exp(lVc + eta.Vc)
      KA <- exp(lKA + eta.KA)
      d/dt(depot) <- -KA*depot
      d/dt(central) <- KA*depot - Cl/Vc*central
      cp <- central/Vc
      cp ~ prop(prop.err)
    })
  }

  d_sim_prep <-
    data.frame(
      AMT=c(10, rep(0, 10)),
      CMT=c("depot", rep("central", 10)),
      TIME=c(0, 0:9),
      ID=1
    )

  expect_error(.nlmixr(f, data=d_sim_prep, est="simulate"), NA)

})
