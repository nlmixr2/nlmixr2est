test_that("estimated parameters that cannot reach the objective are found (#1382)", {
  .unused <- function(f) {
    .unusedEstimatedParams(rxode2::rxUiDecompress(rxode2::rxode2(f)))
  }
  dead <- function() {
    ini({
      tka <- 0.5
      tcl <- -3
      tv <- 3
      tdead <- 1
      eta.cl ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      junk <- tdead * 2
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  expect_equal(.unused(dead), "tdead")

  used <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      tf <- 0.1
      tlag <- -1
      tc <- 2
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      fd <- exp(tf)
      if (tc > 1) {
        ka <- ka * 1.1
      }
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      f(depot) <- fd
      alag(depot) <- exp(tlag)
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  expect_equal(.unused(used), character(0))

  # a theta and an eta behind an intermediate nothing reads; fixed ones are skipped
  chain <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      tfoo <- 1
      tfix <- fix(2)
      eta.ka ~ 0.6
      eta.z ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl)
      v <- exp(tv)
      z <- tfoo + eta.z
      z2 <- z * tfix
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  expect_equal(.unused(chain), c("tfoo", "eta.z"))

  # a parameter used only in a general log-likelihood endpoint is used
  llMod <- function() {
    ini({
      tcl <- 1
      tv <- 3.45
      tsd <- 0.1
    })
    model({
      cl <- exp(tcl)
      v <- exp(tv)
      cp <- linCmt()
      ll(err) ~ -log(tsd) - 0.5 * ((DV - cp) / tsd)^2
    })
  }
  expect_equal(.unused(llMod), character(0))

  # linCmt() as the endpoint itself
  linEnd <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      linCmt() ~ add(add.sd)
    })
  }
  expect_equal(.unused(linEnd), character(0))

  # a state nothing reads does not keep its parameters
  spare <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      tdead <- 0.2
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      d / dt(spare) <- tdead
      cp <- center / v
      cp ~ add(add.sd)
    })
  }
  expect_equal(.unused(spare), "tdead")

  # an unread linCmt() keeps nothing
  linJunk <- function() {
    ini({
      tcl <- 1
      tv <- 3.45
      tb <- 1
      add.sd <- 0.7
    })
    model({
      cl <- exp(tcl)
      v <- exp(tv)
      junk <- linCmt()
      cp <- tb + 0 * t
      cp ~ add(add.sd)
    })
  }
  expect_equal(.unused(linJunk), c("tcl", "tv"))

  # nor does it in a model with its own ODE; a linCmt() state endpoint does
  linOde <- function() {
    ini({
      tcl <- 1
      tv <- 3.45
      tk <- 0.1
      add.sd <- 0.7
    })
    model({
      cl <- exp(tcl)
      v <- exp(tv)
      junk <- linCmt()
      d / dt(eff) <- -exp(tk) * eff
      eff ~ add(add.sd)
    })
  }
  expect_equal(.unused(linOde), c("tcl", "tv"))
  linCentral <- function() {
    ini({
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      cl <- exp(tcl)
      v <- exp(tv)
      junk <- linCmt()
      central ~ add(add.sd)
    })
  }
  expect_equal(.unused(linCentral), character(0))
})

test_that("a parameter with a prior affects the objective", {
  pri <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      tdead <- 1
      add.sd <- 0.7
      prior(tdead) ~ dnorm(2, 1)
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      junk <- tdead * 2
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  .ui <- rxode2::rxUiDecompress(rxode2::rxode2(pri))
  skip_if_not("prior" %in% names(.ui$iniDf), "this lotri has no prior support")
  expect_equal(.unusedEstimatedParams(.ui), character(0))
  expect_equal(.unusedEstimatedParams(.ui, list(priorMethod = "none")), "tdead")
})

test_that("unused parameter hook warns by name", {
  dead <- function() {
    ini({
      tka <- 0.5
      tcl <- -3
      tv <- 3
      tdead <- 1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      junk <- tdead * 2
      cp <- linCmt()
      cp ~ add(add.sd)
    })
  }
  .ui <- rxode2::rxode2(dead)
  expect_warning(
    .preProcessUnusedParams(.ui, "focei", NULL, foceiControl()),
    "'tdead' does not affect the objective function"
  )
  expect_no_warning(.preProcessUnusedParams(.ui, "rxSolve", NULL, NULL))
  # a registered likelihood contributor may read any parameter
  .Call("_nlmixr2est_registerTestContrib", PACKAGE = "nlmixr2est")
  on.exit(.Call("_nlmixr2est_removeTestContrib", PACKAGE = "nlmixr2est"), add = TRUE)
  expect_no_warning(.preProcessUnusedParams(.ui, "focei", NULL, foceiControl()))
})
