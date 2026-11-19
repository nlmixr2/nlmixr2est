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
})
