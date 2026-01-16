test_that("fix parameter saem (#445)", {

  theo_sd2 <- nlmixr2data::theo_sd

  theo_sd2$lwt<-log(theo_sd2$WT/70)

  # The basic model consists of an ini block that has initial estimates
  one.compartment <- function() {
    ini({
      tka <- log(1.57); label("Ka")
      tcl <- log(2.72); label("Cl")
      tv <- log(31.5); label("V")
      covwt<- fix(0.01)
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + covwt*lwt)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  fit0 <- .nlmixr(one.compartment, theo_sd2, est="saem",
                  saemControl(print = 0, seed = 1234, nBurn = 1, nEm = 1,
                              calcTables = FALSE))

  fit1 <- .nlmixr(one.compartment |>
                    ini(covwt=fix(100)), theo_sd2, est="saem",
                  saemControl(print = 0, seed = 1234, nBurn = 1, nEm = 1,
                              calcTables = FALSE))

  theta0 <- fit0$theta
  theta0 <- theta0[names(theta0) != "covwt"]

  theta1 <- fit1$theta
  theta1 <- theta1[names(theta1) != "covwt"]

  expect_true(!all(theta1 == theta0))
})
