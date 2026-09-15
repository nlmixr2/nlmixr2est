nmTest({
  .src <- file.path("..", "..", "src")

  test_that("sampling seeds are sequential, never hashed", {
    skip_if(!dir.exists(.src), "source tree not available (installed package)")
    for (.f in file.path(.src, c("saem.cpp", "npb.cpp"))) {
      skip_if(!file.exists(.f))
      .txt <- paste(readLines(.f, warn = FALSE), collapse = "\n")
      expect_false(grepl("2654435761", .txt, fixed = TRUE), info = .f)
      expect_true(grepl("nmSeqSeedStart(", .txt, fixed = TRUE), info = .f)
    }
    # FOCEi's etaRestart draws take sequential seeds too
    .inner <- paste(readLines(file.path(.src, "inner.cpp"), warn = FALSE), collapse = "\n")
    expect_false(grepl("2654435761", .inner, fixed = TRUE))
    expect_true(grepl("nmSeqSeed(", .inner, fixed = TRUE))
  })

  test_that("saem seeds are sequential by iteration, step and individual", {
    # (nu1B, nphi1, nphi0, nMix, nM, nmc, ntotal): no mode 1B, mode 1B, and a
    # 3-component mixture with no phi0 block
    for (.a in list(c(0L, 3L, 2L, 1L, 4L, 2L, 5L), c(3L, 3L, 2L, 1L, 4L, 2L, 5L),
                    c(0L, 2L, 0L, 3L, 6L, 3L, 4L))) {
      .s <- saemSeedLayoutTest_(c(2L, 2L, 2L), .a[1], .a[2], .a[3], .a[4], .a[5],
                                .a[6], .a[7], 4L)
      # in draw order every seed is the next one: distinct, dense, in order
      expect_equal(.s, seq(0, length(.s) - 1))
      # an iteration's seeds do not depend on how many iterations run
      .s2 <- saemSeedLayoutTest_(c(2L, 2L, 2L), .a[1], .a[2], .a[3], .a[4], .a[5],
                                 .a[6], .a[7], 2L)
      expect_equal(.s[seq_along(.s2)], .s2)
    }
  })

  test_that("a seeded saem fit does not depend on the thread count", {
    skip_on_cran()
    one.compartment <- function() {
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
        linCmt() ~ add(add.sd)
      })
    }
    .fit <- function(threads) {
      .old <- rxode2::getRxThreads(verbose = FALSE)
      on.exit(rxode2::setRxThreads(.old))
      rxode2::setRxThreads(threads)
      if (threads > 1L) skip_if(rxode2::getRxThreads(verbose = FALSE) < threads)
      suppressMessages(nlmixr2(one.compartment, theo_sd, est = "saem",
                               control = saemControl(print = 0, nBurn = 10, nEm = 10,
                                                     seed = 42L, calcTables = FALSE,
                                                     covMethod = "")))
    }
    .f1 <- .fit(1L)
    .f2 <- .fit(2L)
    # the setup solve advances rxode2's seed sequence by the thread count; the
    # kernel restarts it, so the draws, and the fit, are identical
    expect_equal(.f1$objf, .f2$objf)
    expect_equal(.f1$theta, .f2$theta)
  })
})
