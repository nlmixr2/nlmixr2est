# est = "vae", "emvi", "fbvi" and "saem" with a record that is not an observation
# at the time of one.  An EVID 0 record with a missing DV or MDV = 1 is not an
# observation; once translated it is the same event as an EVID 2 record, so the
# three fits must be the same fit, and near the one without the record.  The
# record adds a solver stop, which moved the estimates at most 3.0e-3 relative for
# vae, 2.9e-5 for emvi and 8.4e-8 for saem; the tolerances are about 3x, 30x and
# 100x that.  The vae test pins the missing-DV record that ended its fit on missing
# values; emvi, fbvi and saem are pinned as reading their data the same way.  Also
# here: a subject with no observations (vae and emvi stopped on a subject-count
# mismatch) and data out of time order (vae paired the observations in data order,
# moving its estimates by 1.2e-2; emvi and fbvi solve the translated data).

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
      expect_equal(.fits[[.n]]$theta, .fits$base$theta, tolerance = tol)
      expect_equal(.fits[[.n]]$omega, .fits$base$omega, tolerance = tol)
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

  .ovVaeControl <- function() {
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
    )
  }

  test_that("a subject with no observations is dropped, with a warning, as FOCEi does", {
    skip_on_cran()
    .none <- .ovBase
    .none$DV[.none$ID == 12 & .none$EVID == 0] <- NA
    .without <- .ovBase[.ovBase$ID != 12, ]
    for (.c in list(list("vae", .ovVaeControl()), list("emvi", emviControl(iters = 150L, seed = 7L, print = 0L)))) {
      .f <- suppressMessages(nlmixr2(.ovMod, .none, est = .c[[1]], control = .c[[2]]))
      .w <- suppressMessages(suppressWarnings(nlmixr2(.ovMod, .without, est = .c[[1]], control = .c[[2]])))
      expect_true("IDs without observations dropped: 12" %in% .f$runInfo, label = .c[[1]])
      expect_equal(.f$theta, .w$theta, label = .c[[1]])
      expect_equal(.f$omega, .w$omega, label = .c[[1]])
    }
  })

  test_that("est = \"vae\", \"emvi\" and \"fbvi\" give the same fit for data out of time order", {
    skip_on_cran()
    # each subject's observations in reverse time order
    .u <- do.call(
      rbind,
      lapply(split(.ovBase, .ovBase$ID), function(.s) {
        rbind(.s[.s$EVID != 0, ], .s[rev(which(.s$EVID == 0)), ])
      })
    )
    .ctl <- list(
      vae = .ovVaeControl(),
      emvi = emviControl(iters = 150L, seed = 7L, print = 0L),
      fbvi = fbviControl(iters = 150L, seed = 7L, print = 0L)
    )
    for (.est in names(.ctl)) {
      .s <- suppressMessages(suppressWarnings(nlmixr2(.ovMod, .ovBase, est = .est, control = .ctl[[.est]])))
      .r <- suppressMessages(suppressWarnings(nlmixr2(.ovMod, .u, est = .est, control = .ctl[[.est]])))
      expect_equal(.r$theta, .s$theta, label = .est)
      expect_equal(.r$omega, .s$omega, label = .est)
    }
  })
})
