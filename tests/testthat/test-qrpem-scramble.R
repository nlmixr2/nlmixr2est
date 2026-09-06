# impmapControl(qrScramble=) -- Owen / linear-matrix scrambling of the QRPEM
# Sobol point set.  A Cranley-Patterson shift randomizes the set but leaves the
# correlation structure between high-order dimensions intact; scrambling
# permutes the digits and breaks it.
#
# Mostly pure-kernel tests through the impQrPoints_ hook, so no fits: cheap
# enough for the essential (non-slow) subset.  The load-bearing assertion is
# that qrScramble="none" is bit-identical to the unscrambled path, since every
# existing qr fit depends on it.
nmTest({
  test_that("qrScramble validation and round-trip", {
    expect_identical(impmapControl()$qrScramble, "none")
    expect_identical(impmapControl(qrScramble = "owen")$qrScramble, "owen")
    expect_identical(impmapControl(qrScramble = "lms")$qrScramble, "lms")
    expect_identical(qrpemControl(qrScramble = "owen")$qrScramble, "owen")
    expect_error(impmapControl(qrScramble = "tezuka"))

    .ctl <- impmapControl(qrScramble = "owen")
    expect_identical(do.call(impmapControl, .ctl)$qrScramble, "owen")

    expect_true("qrScramble" %in% .impmapIsControlNames)
    expect_true("qrScramble" %in% .npInertImpCtl)
  })

  test_that("qrScramble='none' is bit-identical to the unscrambled point set", {
    # the default argument and the explicit token must both be the historical
    # path, byte for byte -- every existing qr=TRUE fit rides on this
    .Z0 <- impQrPoints_(256L, 3L, NULL)
    expect_identical(impQrPoints_(256L, 3L, NULL, "none"), .Z0)
    expect_identical(impQrPoints_(256L, 3L, NULL, "none", 999L), .Z0)
    # the shift path is likewise untouched
    .sh <- c(0.371, 0.842, 0.117)
    expect_identical(impQrPoints_(256L, 3L, .sh, "none"),
                     impQrPoints_(256L, 3L, .sh))
  })

  test_that("scrambled points stay a valid stratified N(0,1) set", {
    .n <- 1024L
    .bins <- function(Z, nb = 16L) {
      .U <- pnorm(Z)
      max(vapply(seq_len(ncol(.U)), function(j) {
        max(abs(table(cut(.U[, j], breaks = seq(0, 1, by = 1 / nb))) - .n / nb))
      }, numeric(1)))
    }
    for (.s in c("owen", "lms")) {
      .Z <- impQrPoints_(.n, 4L, NULL, .s, 42L)
      expect_true(all(is.finite(.Z)))
      expect_true(all(dim(.Z) == c(.n, 4L)))
      # a nested digit permutation / linear matrix scramble is
      # stratification-preserving, so the equidistribution survives
      expect_lte(.bins(.Z), 2)
      # standard normal to quasi-random accuracy
      expect_lt(max(abs(colMeans(.Z))), 0.01)
      expect_equal(unname(apply(.Z, 2, sd)), rep(1, 4L), tolerance = 0.02)
    }
  })

  test_that("scrambling preserves the 2-D net property", {
    # A 1-D marginal check CANNOT see a wrong-orientation linear scramble: over
    # the first 2^m points the low digits are a linear function of the top m, so
    # an upper-triangular L still permutes them and stays 1-D stratified.  It
    # does destroy the (t,m,s)-net elementary-interval balance, which is where
    # the QMC accuracy actually lives -- an earlier draft of impLmsScramble had
    # exactly that bug and left 512-767 empty boxes below.
    .net <- function(Z) {
      .U <- pnorm(Z)
      .m <- 10L
      .empty <- 0L
      for (.d1 in 0:.m) {
        .d2 <- .m - .d1
        .b1 <- pmin(floor(.U[, 1] * 2^.d1), 2^.d1 - 1)
        .b2 <- pmin(floor(.U[, 2] * 2^.d2), 2^.d2 - 1)
        .tb <- table(factor(.b1 * 2^.d2 + .b2, levels = 0:(2^.m - 1)))
        .empty <- max(.empty, sum(.tb == 0))
      }
      .empty
    }
    # the raw sequence leaves 30 (boost skips the zero point); a correct
    # scramble leaves 1.  A broken one leaves hundreds.
    expect_lte(.net(impQrPoints_(1024L, 2L, NULL)), 40L)
    for (.s in c("owen", "lms")) {
      for (.seed in c(42L, 1L, 12345L)) {
        expect_lte(.net(impQrPoints_(1024L, 2L, NULL, .s, .seed)), 4L)
      }
    }
  })

  test_that("scrambling is seeded, reproducible and distinct per method", {
    .a <- impQrPoints_(512L, 3L, NULL, "owen", 42L)
    .b <- impQrPoints_(512L, 3L, NULL, "owen", 42L)
    .c <- impQrPoints_(512L, 3L, NULL, "owen", 43L)
    .d <- impQrPoints_(512L, 3L, NULL, "lms", 42L)
    .n0 <- impQrPoints_(512L, 3L, NULL)
    # same seed reproduces exactly (this is what makes a fit reproducible)
    expect_identical(.a, .b)
    # a different seed gives a genuinely different randomization ...
    expect_false(identical(.a, .c))
    # ... with the same first two moments
    expect_equal(colMeans(.a), colMeans(.c), tolerance = 0.02)
    # the two methods are not the same map, and neither is the identity
    expect_false(identical(.a, .d))
    expect_false(identical(.a, .n0))
    expect_false(identical(.d, .n0))
  })

  test_that("a scramble replaces the shift rather than composing with it", {
    expect_error(impQrPoints_(64L, 2L, c(0.1, 0.2), "owen", 42L), "shift")
    expect_error(impQrPoints_(64L, 2L, c(0.1, 0.2), "lms", 42L), "shift")
  })

  test_that("qrScramble reaches the kernel and keeps a fit thread-independent", {
    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- 1; tv <- 3.45
        eta.ka ~ 0.6; eta.cl ~ 0.3
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }
    .dat <- nlmixr2data::theo_sd
    .thr0 <- rxode2::getRxThreads()
    on.exit(rxode2::setRxThreads(.thr0), add = TRUE)
    .run <- function(scr, nthr) {
      rxode2::setRxThreads(nthr)
      suppressWarnings(suppressMessages(
        nlmixr2(one.cmt, .dat, "qrpem",
                qrpemControl(print = 0L, nIter = 3L, isample = 200L,
                             qrScramble = scr, covMethod = "",
                             calcTables = FALSE))))
    }
    .none <- .run("none", 1L)
    .owen1 <- .run("owen", 1L)
    .owen2 <- .run("owen", 2L)

    # the resolved mode is reported back
    expect_identical(.none$env$impQrScramble, "none")
    expect_identical(.owen1$env$impQrScramble, "owen")
    # scrambling actually changed the sampler
    expect_false(isTRUE(all.equal(.owen1$env$impObjTrace, .none$env$impObjTrace,
                                  tolerance = 1e-8)))
    # ... and the fit is still bit-identical across thread counts, because the
    # scramble key is arithmetic in (seed, subject, iteration, dimension) rather
    # than drawn from the RNG
    expect_equal(.owen2$env$impObjTrace, .owen1$env$impObjTrace, tolerance = 1e-12)
    expect_equal(.owen2$objf, .owen1$objf, tolerance = 1e-12)
  })
})
