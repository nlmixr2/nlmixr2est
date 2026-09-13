# etaDistInit() converts between a declared family and its log-normal surrogate
# TWICE -- forward to seed the surrogate, back to read its answer.  Both
# conversions used the ARITHMETIC mean/variance, and for a heavy-tailed family
# that is catastrophically wrong: a log-normal and a gamma sharing a log-variance
# differ by ~70x in relative variance at gamma shape 0.5.
#
# The failure was silent end to end.  The surrogate fit converges, the solver
# reports success (it DID match the moments it was handed), and every returned
# value is finite -- so nothing downstream could tell.  Measured on Bauer's g4
# it seeded relative variance 68 against a truth of 2, and the fit it warm-started
# ran away to CL 2312 (MARE 9132%) where the cold start reached 25.9%.
#
# No fitting here: these are closed-form identities, so the test is exact and fast.
nmTest({

  test_that("a declared family reports its log-scale moments", {
    # For X ~ gamma(shape a, rate b):  E[log X] = digamma(a) - log(b)
    #                                  Var[log X] = trigamma(a)
    .a <- 0.5
    .m <- exp(1.63)
    .b <- .a / .m
    .d <- sprintf("dgamma(shape = %.17g, rate = %.17g)", .a, .b)
    .mv <- .etaDistMoments(.d, list(), .etaDistGh(60L))
    expect_false(is.null(.mv))
    # arithmetic pair, unchanged behaviour
    expect_equal(.mv[["mean"]], .m, tolerance = 1e-3)
    expect_equal(.mv[["var"]], .m^2 / .a, tolerance = 1e-2)
    # log pair, the one the surrogate needs
    expect_equal(.mv[["meanlog"]], digamma(.a) - log(.b), tolerance = 1e-3)
    expect_equal(.mv[["varlog"]], trigamma(.a), tolerance = 1e-2)
  })

  test_that("matching a surrogate's log moments recovers the declared family", {
    # the surrogate fitted on g4: meanlog 0.4755, varlog 4.2386
    .mu <- 0.4755
    .w <- 4.2386
    .dist <- "dgamma(shape = 1/exp(lclrv), rate = 1/(exp(lclrv) * exp(lclm)))"
    .start <- c(lclm = 1.5, lclrv = 0.7)
    .gh <- .etaDistGh(60L)

    .sol <- .etaDistSolveThetas(
      .dist, c("lclm", "lclrv"), .start,
      c(mean = exp(.mu + 0.5 * .w), var = (exp(.w) - 1) * exp(2 * .mu + .w),
        meanlog = .mu, varlog = .w), .gh)
    expect_false(is.null(.sol))
    # truth for g4 is CL 5.104, relative variance 2.0
    expect_equal(exp(.sol[["lclm"]]), 5.104, tolerance = 0.15)
    expect_equal(exp(.sol[["lclrv"]]), 2.0, tolerance = 0.20)

    # and the arithmetic route, which is what shipped, is nowhere near:
    # exp(varlog) - 1 is the relative variance of the LOG-NORMAL, not of the
    # gamma that has to reproduce it
    expect_gt(exp(.w) - 1, 50)          # what the old target asked for
    expect_lt(exp(.sol[["lclrv"]]), 5)  # what the family actually needs
  })

  test_that("a warm start that is no improvement is refused", {
    # the guard compares the proposal against the model's OWN starting values on
    # the surrogate's criterion, so it needs no threshold for what counts as an
    # implausible parameter -- which matters, since dist() takes any family and a
    # relative variance of 68 can be honest for some of them
    .dist <- "dgamma(shape = 1/exp(lclrv), rate = 1/(exp(lclrv) * exp(lclm)))"
    .gh <- .etaDistGh(60L)
    # mirrors the guard in etaDistInit(): Inf when the family cannot report
    # log-scale moments at all, which is itself a rejection
    .fitTo <- function(v, mu, w) {
      .mm <- .etaDistMoments(.dist, as.list(v), .gh)
      if (is.null(.mm) || !all(is.finite(.mm[c("meanlog", "varlog")])) ||
            .mm[["varlog"]] <= 0) return(Inf)
      (.mm[["meanlog"]] - mu)^2 + (log(.mm[["varlog"]] / w))^2
    }
    # g4's surrogate, against g4's own ini(): the solved values must be closer
    .ini <- c(lclm = 1.5, lclrv = 0.7)
    .good <- c(lclm = 1.6195, lclrv = 0.6028)
    expect_lt(.fitTo(.good, 0.4755, 4.2386), .fitTo(.ini, 0.4755, 4.2386))
    # What the arithmetic route produced is FARTHER away, so the guard drops it.
    # It is in fact so far away that the implied family (gamma shape 0.015) has
    # no finite log-scale moments at all -- its quantiles reach zero -- which the
    # guard scores as Inf and rejects for the same reason.
    .bad <- c(lclm = 2.5939, lclrv = 4.2263)
    expect_gt(.fitTo(.bad, 0.4755, 4.2386), .fitTo(.ini, 0.4755, 4.2386))
  })

})
