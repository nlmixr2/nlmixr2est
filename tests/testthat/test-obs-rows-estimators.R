# est = "vae", "emvi", "fbvi" and "saem" with a record that is not an observation
# at the time of one.  An EVID 0 record with a missing DV or MDV = 1 is not an
# observation; once translated it is the same event as an EVID 2 record, so the
# three fits must agree, and stay near the fit without the record (the extra
# solver stop moved the estimates 3.0e-3 relative for vae, 2.9e-5 for emvi and
# 8.4e-8 for saem).

nmTest({
  .ovMod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  .ovBase <- nlmixr2data::theo_sd
  .ovBase$MDV <- as.integer(.ovBase$EVID != 0)
  .ovT <- .ovBase$TIME[.ovBase$ID == 1 & .ovBase$EVID == 0][3]
  .ovAdd <- function(evid, dv, mdv) {
    .r <- data.frame(
      ID = 1L,
      TIME = .ovT,
      DV = dv,
      AMT = 0,
      EVID = evid,
      CMT = 2L,
      WT = .ovBase$WT[.ovBase$ID == 1][1],
      MDV = mdv
    )
    .d <- rbind(.ovBase, .r)
    .d[order(.d$ID, .d$TIME, -.d$EVID), ]
  }
  .ovData <- list(
    base = .ovBase,
    evid2 = .ovAdd(2L, NA_real_, 1L),
    dvNa = .ovAdd(0L, NA_real_, 0L),
    mdv1 = .ovAdd(0L, 5, 1L)
  )

  .ovExpect <- function(est, control, tol) {
    .fits <- lapply(.ovData, function(.d) {
      suppressMessages(suppressWarnings(nlmixr2(.ovMod, .d, est = est, control = control)))
    })
    for (.n in names(.fits)) {
      expect_true(all(is.finite(.fits[[.n]]$theta)))
      expect_true(all(is.finite(.fits[[.n]]$omega)))
      expect_equal(.fits[[.n]]$theta, .fits$base$theta, tolerance = tol)
    }
    for (.n in c("dvNa", "mdv1")) {
      expect_equal(.fits[[.n]]$theta, .fits$evid2$theta)
      expect_equal(.fits[[.n]]$omega, .fits$evid2$omega)
    }
  }

  test_that("est = \"vae\" does not count a missing-DV or MDV = 1 record as an observation", {
    skip_on_cran()
    .ovExpect(
      "vae",
      vaeControl(
        itersBurnIn = 40L,
        klWarmup = 20L,
        gammaIter = 60L,
        iters = 80L,
        hiddenDim = 25L,
        seed = 1L,
        covariateSelection = FALSE,
        print = 0L,
        covMethod = ""
      ),
      tol = 1e-2
    )
  })

  test_that("est = \"emvi\" and \"fbvi\" do not either", {
    skip_on_cran()
    .ovExpect("emvi", emviControl(iters = 150L, seed = 7L, print = 0L), tol = 1e-3)
    .ovExpect("fbvi", fbviControl(iters = 150L, seed = 7L, print = 0L), tol = 1e-3)
  })

  test_that("est = \"saem\" does not either", {
    skip_on_cran()
    .ovExpect("saem", saemControl(nBurn = 50L, nEm = 50L, print = 0L, covMethod = ""), tol = 1e-5)
  })
})
