nmTest({
  test_that("foceiControl() hessianMethod validation", {
    expect_equal(foceiControl()$hessianMethod, 3L)
    expect_equal(foceiControl(hessianMethod = "fd")$hessianMethod, 1L)
    expect_equal(foceiControl(hessianMethod = "bfgs")$hessianMethod, 2L)
    expect_equal(foceiControl(hessianMethod = "sr1")$hessianMethod, 3L)
    expect_equal(foceiControl(hessianMethod = "bofill")$hessianMethod, 4L)
    expect_equal(foceiControl(hessianMethod = 3L)$hessianMethod, 3L)
    expect_error(foceiControl(hessianMethod = "not-a-method"))
  })

  # A non-normal (Poisson) endpoint: the per-subject inner Hessian has no
  # Gaussian Gauss-Newton shortcut, so calcEtaHessian()'s needOptimHess
  # branch (finite-difference by default, or hessianMethod=)'s quasi-Newton
  # update) is the ONLY inner-Hessian path exercised here.
  .poisMod <- function() {
    ini({
      te0 <- log(5)
      eta.e0 ~ 0.3
    })
    model({
      e0 <- exp(te0 + eta.e0)
      effect <- e0
      effect ~ dpois(effect)
    })
  }
  .poisData <- data.frame(ID = rep(1:20, each = 3), TIME = rep(c(0, 1, 2), 20))
  .poisData$DV <- c(5L, 4L, 6L, 3L, 5L, 4L, 6L, 5L, 7L, 4L, 5L, 5L,
                     6L, 6L, 5L, 4L, 3L, 5L, 5L, 6L, 4L, 7L, 5L, 6L,
                     4L, 5L, 5L, 3L, 4L, 6L, 6L, 5L, 4L, 5L, 6L, 5L,
                     4L, 4L, 6L, 5L, 7L, 5L, 6L, 5L, 4L, 3L, 5L, 6L,
                     5L, 5L, 6L, 4L, 6L, 5L, 5L, 4L, 5L, 6L, 4L, 5L)

  test_that("hessianMethod converges close to fd on a non-normal endpoint", {
    skip_on_cran()
    .fFd <- .nlmixr(.poisMod, .poisData, est = "focei",
                    control = foceiControl(print = 0L, hessianMethod = "fd"))
    expect_true(is.finite(.fFd$objf))
    for (.hm in c("bfgs", "sr1", "bofill")) {
      .fT <- .nlmixr(.poisMod, .poisData, est = "focei",
                     control = foceiControl(print = 0L, hessianMethod = .hm))
      expect_true(is.finite(.fT$objf), info = .hm)
      expect_equal(.fT$objf, .fFd$objf, tolerance = 1e-2, info = .hm)
      expect_equal(unname(.fT$theta), unname(.fFd$theta), tolerance = 1e-2, info = .hm)
      # Positive evidence the quasi-Newton mechanism actually ran (mirrors
      # test-focei-trust-inner.R's own .nTrustInner() convention).
      expect_true(.nHessianQN() > 0L, info = .hm)
    }
  })

  .oneCmt <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka)
      cl <- exp(tcl)
      v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(centr) <- ka * depot - cl / v * centr
      cp <- centr / v
      cp ~ add(add.sd)
    })
  }

  test_that("hessianMethod has no effect on a normal-endpoint model (needOptimHess gates it off)", {
    skip_on_cran()
    .fN <- .nlmixr(.oneCmt, nlmixr2data::theo_sd, est = "focei",
                   control = foceiControl(print = 0L, hessianMethod = "sr1"))
    expect_true(is.finite(.fN$objf))
    # The Gauss-Newton inner Hessian is used unconditionally for a normal
    # endpoint -- the quasi-Newton mechanism should never engage here.
    expect_equal(.nHessianQN(), 0L)
  })

  test_that("hessianMethod is inert unless innerOpt='trust' (calcEtaHessian() is also reached from warmZm()/LikInner2() for other inner optimizers, where the per-attempt reset does not apply)", {
    skip_on_cran()
    .fN1qn1 <- .nlmixr(.poisMod, .poisData, est = "focei",
                       control = foceiControl(print = 0L, innerOpt = "n1qn1",
                                              hessianMethod = "sr1"))
    expect_true(is.finite(.fN1qn1$objf))
    # The quasi-Newton mechanism must never engage for a non-trust inner
    # optimizer, regardless of hessianMethod -- calcEtaHessian() is reached
    # from warmZm()'s one-time n1qn1 warm-start seed and from the final
    # objective recompute, neither of which is the repeated, reset-per-attempt
    # Newton-step loop the quasi-Newton state assumes.
    expect_equal(.nHessianQN(), 0L)
  })
})
