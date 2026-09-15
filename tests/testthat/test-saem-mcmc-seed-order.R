nmTest({
  # SAEM's sampling draws -- do_mcmc proposals and augmentCensY()'s
  # rxTruncNorm() -- read the per-thread threefry engine (rxNormEng/rxUnifEng,
  # truncNorm.h).  rxode2's par_*() loops re-seed that SAME engine once per
  # subject, setSeedEng1(getRxSeed1() + id), so when a solve returns the engine
  # holds the seed of whichever subject the thread happened to solve LAST --
  # which rx->ordId, the solve order, decides.  Any draw taken after a solve
  # without restoring the sampling block's own seed therefore depends on the
  # solve order rather than on anything keyed to a subject.
  #
  # Two of the eighteen user_fn() call sites wrapped themselves in nmRngGuard();
  # the rest did not.  The restore is taken inside user_function() so that no
  # caller has to remember, and so a nineteenth call site cannot reintroduce it.
  #
  # Measured before the restore was moved: a 131-subject SAEM fit gave four
  # different objective values over eight seeded runs once rx->ordId was a
  # non-identity permutation, and a different value again for each of five solve
  # orders.  After: 8/8 identical, and identical across all five orders.

  .src <- file.path("..", "..", "src")

  test_that("user_function() restores the sampling seed on every exit", {
    skip_if(!dir.exists(.src), "source tree not available (installed package)")
    .f <- file.path(.src, "saem.cpp")
    skip_if(!file.exists(.f))
    .l <- readLines(.f, warn = FALSE)
    .start <- grep("^mat user_function\\(", .l)
    expect_equal(length(.start), 1L)
    # the guard has to be a destructor in the function's own scope, so that an
    # early return or a throw still restores
    .head <- .l[seq(.start, min(.start + 25L, length(.l)))]
    expect_true(any(grepl("nmRestoreMcmcSeed()", .head, fixed = TRUE)))
    expect_true(any(grepl("~_SaemMcmcSeedGuard", .head, fixed = TRUE)))
  })

  test_that("sampling seeds are sequential, never hashed", {
    skip_if(!dir.exists(.src), "source tree not available (installed package)")
    for (.f in file.path(.src, c("saem.cpp", "npb.cpp"))) {
      skip_if(!file.exists(.f))
      .txt <- paste(readLines(.f, warn = FALSE), collapse = "\n")
      expect_false(grepl("2654435761", .txt, fixed = TRUE), info = .f)
      expect_true(grepl("nmSeqSeedReserve(", .txt, fixed = TRUE), info = .f)
      expect_true(grepl("nmSeqSeedStart(", .txt, fixed = TRUE), info = .f)
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

  test_that("saem seeds are sequential by iteration, step and individual", {
    # (nphi1, nphi0, nMix, nM, nmc, ntotal): no mixture, and a 3-component
    # mixture with no phi0 block
    for (.a in list(c(3L, 2L, 1L, 4L, 2L, 5L), c(2L, 0L, 3L, 6L, 3L, 4L))) {
      .s <- saemSeedLayoutTest_(c(2L, 2L, 2L), .a[1], .a[2], .a[3], .a[4], .a[5],
                                .a[6], 4L)
      # in draw order every seed is the next one: distinct, dense, in order
      expect_equal(.s, seq(0, length(.s) - 1))
      # an iteration's seeds do not depend on how many iterations run
      .s2 <- saemSeedLayoutTest_(c(2L, 2L, 2L), .a[1], .a[2], .a[3], .a[4], .a[5],
                                 .a[6], 2L)
      expect_equal(.s[seq_along(.s2)], .s2)
    }
  })
})
