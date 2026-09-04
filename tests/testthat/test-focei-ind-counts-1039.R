nmTest({
  # Issue #1039 -- every per-subject block in the FOCEi setup (gVid, ga/gc, gB,
  # gcH*, llikObsFull) is sized and strided from rxode2's per-subject event
  # counts, read through the rxode2 pointer table.  A build whose rxode2 and
  # nlmixr2est disagree on the solve layout returns those counts from the wrong
  # bytes, and the setup then sizes megabytes of per-subject storage from
  # garbage.  Both checks below need more than one subject: subject 0 sits at
  # the start of the subject array, so it reads correctly under any stride.

  test_that("rxode2's per-subject event counts split the dataset (#1039)", {
    # unequal per-subject record counts, so a wrong stride cannot pass by
    # reading a neighbouring (identical) subject
    nObsI <- c(3, 5, 2, 7, 4, 6)
    d <- do.call(rbind, lapply(seq_along(nObsI), function(id) {
      rbind(data.frame(ID = id, TIME = 0, DV = NA_real_, AMT = 320, EVID = 1,
                       CMT = 1),
            data.frame(ID = id, TIME = seq(0.5, 12, length.out = nObsI[id]),
                       DV = 5, AMT = 0, EVID = 0, CMT = 1))
    }))

    m <- rxode2::rxode2({
      ka <- 1
      cl <- 3
      v <- 30
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
    })
    suppressMessages(rxode2::rxSolve(m, d))

    cnt <- foceiIndEventCounts_()
    expect_equal(nrow(cnt), length(nObsI))
    expect_equal(as.integer(cnt[, "nAllTimes"]), as.integer(nObsI) + 1L)
    expect_equal(as.integer(cnt[, "nDoses"]), rep(1L, length(nObsI)))
    expect_equal(as.integer(cnt[, "nEvid2"]), rep(0L, length(nObsI)))
    # this is the invariant foceiCheckIndCounts() enforces before it sizes
    # anything: the per-subject counts have to add back up to the dataset
    expect_equal(sum(cnt[, "nAllTimes"]), nrow(d))
  })

  test_that("a multi-subject focei fit sizes its per-subject blocks correctly", {
    # the report's own reproduction: on a mismatched build this failed at
    # .foceiFitInternal() with "dataset too large", an R_Calloc failure, or a
    # segfault, depending on what the mis-read bytes happened to hold
    oneCompartment <- function() {
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

    fit <- suppressMessages(suppressWarnings(
      nlmixr(oneCompartment(), nlmixr2data::theo_sd, "focei",
             control = foceiControl(maxOuterIterations = 0, covMethod = "",
                                    print = 0))))
    expect_true(length(unique(fit$ID)) > 1)
    expect_true(is.finite(fit$objf))
    expect_true(all(is.finite(fit$IPRED)))
    # llikObsFull is the one per-subject block handed back to R in record
    # order, so its per-subject offsets come straight from these counts: a
    # subject given the wrong offset leaves its records unwritten and moves the
    # dose records (the NA entries) off the dose rows of the event table
    ll <- fit$llikObs
    expect_equal(length(ll), nrow(fit$dataSav))
    expect_equal(which(is.na(ll)), which(fit$dataSav$EVID != 0))
    expect_true(all(is.finite(ll[!is.na(ll)])))
  })
})
