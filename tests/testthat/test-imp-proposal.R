# impmapControl(proposal=) -- the importance-sampling proposal family.
#
# The algebra is the part most likely to be silently wrong, and it is testable
# without a fit: every family's closed-form normalizer is pinned against NUMERIC
# QUADRATURE rather than against itself, through the impPropKernel_ hook.  Two
# invariants the rest of the kernel is built on get their own tests, because
# breaking either corrupts xi (and therefore both gamma controllers) rather than
# producing a visible error:
#   * the peak-normalized kernel is EXACTLY 0 at quad = 0 -- qCenter omits the
#     kernel term for that reason;
#   * the Ci correction is gamma-free -- that is what lets Ci factor at all.
nmTest({
  .K <- function(type, quad, gamma = 1, p = 2, df = 0,
                 cs = c(1, 9), ws = c(0.9, 0.1)) {
    impPropKernel_(type, df, cs, ws, quad, gamma, p)
  }

  test_that("proposal validation, defaults and round-trip", {
    expect_identical(impmapControl()$proposal, "auto")
    expect_identical(impmapControl(proposal = "laplace")$proposal, "laplace")
    expect_identical(impmapControl(df = 8, proposal = "t")$proposal, "t")
    expect_identical(qrpemControl(proposal = "mixture")$proposal, "mixture")

    # df keeps its meaning and "auto" preserves its historical double duty
    expect_identical(impmapControl(df = 0)$proposal, "auto")
    expect_identical(impmapControl(df = 8)$proposal, "auto")

    # contradictions are refused rather than silently resolved
    expect_error(impmapControl(proposal = "normal", df = 8), "normal")
    expect_error(impmapControl(proposal = "t", df = 0), "df")
    expect_error(impmapControl(proposal = "laplace", df = 8), "df")
    expect_error(impmapControl(proposal = "mixture", df = 3), "df")
    expect_error(impmapControl(proposal = "cauchy"))

    # mixture parameters
    expect_error(impmapControl(proposal = "mixture", propMixScale = 1), "propMixScale")
    expect_error(impmapControl(propMixScale = c(1, 4), propMixWeight = c(1, 1, 1)),
                 "propMixWeight")
    expect_error(impmapControl(propMixScale = c(2, 9)), "start at 1")
    expect_error(impmapControl(propMixScale = c(1, 9, 4),
                               propMixWeight = c(0.5, 0.3, 0.2)), "increasing")
    expect_error(impmapControl(propMixScale = c(1, -3)), "propMixScale")
    # weights are normalized, not required to sum to 1
    expect_equal(impmapControl(propMixWeight = c(3, 1))$propMixWeight, c(0.75, 0.25))

    .ctl <- impmapControl(proposal = "laplace")
    expect_identical(do.call(impmapControl, .ctl)$proposal, "laplace")

    expect_true(all(c("proposal", "propMixScale", "propMixWeight") %in%
                      .impmapIsControlNames))
    expect_true(all(c("proposal", "propMixScale", "propMixWeight") %in%
                      .npInertImpCtl))
  })

  test_that("the kernel is exactly zero at the peak for every family", {
    # NOT approximately: qCenter (src/imp.cpp) drops the kernel term outright,
    # so any family whose kernel is nonzero at d = 0 silently biases xi
    for (.p in 1:5) {
      for (.g in c(0.5, 1, 2.7)) {
        expect_identical(unname(.K("normal", 0, .g, .p)[1]), 0)
        expect_identical(unname(.K("t", 0, .g, .p, df = 7)[1]), 0)
        expect_identical(unname(.K("laplace", 0, .g, .p)[1]), 0)
        expect_identical(unname(.K("mixture", 0, .g, .p)[1]), 0)
      }
    }
  })

  test_that("the Ci correction is gamma-free for every family", {
    for (.ty in c("normal", "t", "laplace", "mixture")) {
      .cc <- vapply(c(0.3, 1, 5),
                    function(g) unname(.K(.ty, 1.3, g, 4, df = 6)[2]), numeric(1))
      expect_identical(diff(range(.cc)), 0)
    }
  })

  test_that("every family's density integrates to one (quadrature, not self)", {
    # g(x) = peak * K(x) with peak = gaussPeak * exp(-corr); with Sigma = I this
    # must integrate to 1 over R^p.  This is what validates the closed-form
    # normalizers independently of the code that produced them.
    .intg <- function(ty, p, df = 0) {
      .corr <- unname(.K(ty, 1, 1, p, df)[2])
      .logPeak <- -0.5 * p * log(2 * pi) - .corr
      .Sp <- 2 * pi^(p / 2) / gamma(p / 2)
      .f <- function(r) {
        vapply(r, function(rr) {
          rr^(p - 1) * exp(-unname(.K(ty, rr^2, 1, p, df)[1]))
        }, numeric(1))
      }
      exp(.logPeak) * .Sp *
        integrate(.f, 0, Inf, rel.tol = 1e-10, subdivisions = 2000L)$value
    }
    for (.p in 1:5) {
      expect_equal(.intg("normal", .p), 1, tolerance = 1e-8)
      expect_equal(.intg("t", .p, df = 6), 1, tolerance = 1e-8)
      expect_equal(.intg("laplace", .p), 1, tolerance = 1e-8)
      expect_equal(.intg("mixture", .p), 1, tolerance = 1e-8)
    }
  })

  test_that("the t correction reduces to the Gaussian as df grows", {
    expect_equal(unname(.K("t", 2, 1, 3, df = 1e7)[2]), 0, tolerance = 1e-5)
    # NEGATIVE at low df, and that is correct: a small-df t concentrates more
    # mass at its center than the Gaussian with the same scale matrix, so its
    # peak density is HIGHER and log(gaussPeak) - log(tPeak) < 0.  The
    # quadrature test above is what proves the value itself.
    expect_lt(unname(.K("t", 2, 1, 3, df = 4)[2]), 0)
    # monotone toward the Gaussian as the tail lightens
    expect_gt(unname(.K("t", 2, 1, 3, df = 30)[2]),
              unname(.K("t", 2, 1, 3, df = 4)[2]))
    # and the normal family's correction is exactly zero, which is what makes
    # the xi normalization a no-op there
    expect_identical(unname(.K("normal", 2, 1, 3)[2]), 0)
  })

  test_that("the Laplace scale is covariance matched", {
    # Cov = (p+1) * S, so S = Sigma/(p+1) keeps gamma/iscaleMin/iscaleMax
    # meaning what they mean for the normal proposal
    for (.p in 1:5) {
      expect_equal(unname(.K("laplace", 1, 1, .p)[3]), 1 / (.p + 1),
                   tolerance = 1e-12)
    }
  })

  test_that("the mixture kernel uses the mixture density, not a component", {
    # the single most likely wrong-but-plausible implementation.  A single
    # Gaussian component would give quad/(2 c gamma); the mixture is strictly
    # LESS than the narrow component's kernel at large quad, because the wide
    # component carries the tail.
    .p <- 3
    .narrow <- 1 / (2 * 1) * 20        # component 1 (c = 1) at quad = 20
    .mix <- unname(.K("mixture", 20, 1, .p, cs = c(1, 9), ws = c(0.9, 0.1))[1])
    expect_lt(.mix, .narrow)
    # ... and strictly MORE than the wide component alone
    expect_gt(.mix, 20 / (2 * 9))
    # monotone in quad, and finite far out in the tail where a naive
    # log-sum-exp would underflow to log(0) and return +Inf
    .far <- unname(.K("mixture", 1e6, 1, .p)[1])
    expect_true(is.finite(.far))
    expect_gt(.far, .mix)
  })
})
