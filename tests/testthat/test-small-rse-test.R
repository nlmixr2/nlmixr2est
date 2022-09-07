test_that("rse for focei does not inflate", {

  u002 <- function() {
    ini({ # Where initial conditions/variables are specified
      # '<-' or '=' defines population parameters
      # Simple numeric expressions are supported
      lCl <- 1.6      #log Cl (L/hr)
      lVc <- 4.5      #log V (L)
      # Bounds may be specified by c(lower, est, upper), like NONMEM:
      # Residuals errors are assumed to be population parameters
      prop.err <- c(0, 0.3, 1)
      # Between subject variability estimates are specified by '~'
      # Semicolons are optional
      eta.Vc ~ 0.1   #IIV V
      eta.Cl ~ 0.1   #IIV Cl
    })
    model({ # Where the model is specified
      # The model uses the ini-defined variable names
      Vc <- exp(lVc + eta.Vc)
      Cl <- exp(lCl + eta.Cl)
      # RxODE-style differential equations are supported
      d / dt(centr) <- -(Cl / Vc) * centr
      ## Concentration is calculated
      cp <- centr / Vc
      # And is assumed to follow proportional error estimated by prop.err
      cp ~ prop(prop.err)
    })
  }
  datr <- nlmixr2data::Bolus_1CPT
  dat <- datr[datr$SD == 0, ]
  dat <- dat[, names(dat) != "SS"]

  f <- nlmixr2(u002, dat, "focei", control=foceiControl(print=0))

  expect_true(f$parFixedDf["lVc", "SE"] < 0.029)
  
  
})
