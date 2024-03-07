test_that("compressed ui in nlmixr2 fit doesn't interefere with solve (#422)", {
  nimo <- function() {
    ini({
      ## Note that the UI can take expressions
      ## Also note that these initial estimates should be provided on the log-scale
      tcl <- log(0.001)
      tv1 <- log(1.45)
      tQ <- log(0.004)
      tv2 <- log(44)
      tkss <- log(12)
      tkint <- log(0.3)
      tksyn <- log(1)
      tkdeg <- log(7)
      ## Initial estimates should be high for SAEM ETAs
      eta.cl  ~ 2
      eta.v1  ~ 2
      eta.kss ~ 2
      ##  Also true for additive error (also ignored in SAEM)
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

  fit <- nlmixr(nimo, nlmixr2data::nimoData, est="saem", control=saemControl(print=0, nBurn=2, nEm=2))
  expect_error(rxSolve(fit, nlmixr2data::nimoData), NA)
})
