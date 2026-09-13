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

  test_that("the sampling-block seeds are keyed, not taken from the engine", {
    skip_if(!dir.exists(.src), "source tree not available (installed package)")
    .f <- file.path(.src, "saem.cpp")
    skip_if(!file.exists(.f))
    .txt <- paste(readLines(.f, warn = FALSE), collapse = "\n")
    # every nmSetSeedEng1() argument is folded from the fit's own seed plus
    # iteration/chain/endpoint indices -- never from a running counter, a thread
    # id, or a subject the solve order picked
    expect_true(grepl("_saemSeedDoMcmc", .txt, fixed = TRUE))
    expect_true(grepl("_saemSeedCensAug", .txt, fixed = TRUE))
    expect_true(grepl("s = s * 2654435761u + (uint32_t)kiter;", .txt, fixed = TRUE))
  })
})
