nmTest({
  skip_if_not(.nlmixr2estRxode2HasTte(), "rxode2 without tte() endpoints")

  # Gauss-Hermite nodes and weights for a standard normal
  .tteGh <- function(n = 40) {
    .J <- matrix(0, n, n)
    for (.i in seq_len(n - 1)) {
      .J[.i, .i + 1] <- .J[.i + 1, .i] <- sqrt(.i)
    }
    .e <- eigen(.J, symmetric = TRUE)
    list(x = .e$values, w = .e$vectors[1, ]^2)
  }

  # exact marginal -log-likelihood of a frailty model; `ll(i, eta)` is subject
  # i's log-likelihood at the eta nodes
  .tteMarginalNll <- function(nsub, ll) {
    .gh <- .tteGh()
    function(p) {
      -sum(vapply(
        seq_len(nsub),
        function(i) {
          .l <- ll(p, i, .gh$x * exp(p[length(p)] / 2))
          max(.l) + log(sum(.gh$w * exp(.l - max(.l))))
        },
        numeric(1)
      ))
    }
  }

  .tteFrailtyData <- function(mod) {
    .s <- withr::with_seed(7, rxode2::rxSolve(mod, rxode2::et(c(0, 20)) |> rxode2::et(id = 1:150)))
    .s <- as.data.frame(.s)
    data.frame(ID = .s$id, TIME = .s$time, DV = .s$sim, EVID = 0)
  }

  test_that("focei/laplace/agq keep the frailty variance of a repeated-event tte() model (#1465)", {
    skip_on_cran()
    .mod <- function() {
      ini({
        tlam <- log(0.3)
        eta.lam ~ 0.2
      })
      model({
        haz <- exp(tlam + eta.lam)
        ev ~ tte(haz)
      })
    }
    .d <- .tteFrailtyData(.mod())
    .n <- tapply(.d$DV, .d$ID, sum)
    .end <- tapply(.d$TIME, .d$ID, max)
    .ll <- function(p, i, eta) .n[i] * (p[1] + eta) - exp(p[1] + eta) * .end[i]
    .nll <- .tteMarginalNll(length(.n), .ll)
    .o <- stats::optim(c(log(0.5), log(0.1)), .nll, method = "BFGS")$par
    .start <- .mod |> ini(tlam = log(0.5), eta.lam ~ 0.1)
    for (.est in c("focei", "laplace", "agq")) {
      .fit <- suppressMessages(.nlmixr(.start, .d, est = .est, control = list(print = 0)))
      expect_equal(unname(.fit$theta[["tlam"]]), .o[1], tolerance = 0.02, label = .est)
      expect_equal(unname(.fit$omega[1, 1]), exp(.o[2]), tolerance = 0.1, label = .est)
    }
    # with the variance fixed at the truth, focei's rate is the exact one
    .oFix <- stats::optimize(function(t) .nll(c(t, log(0.2))), c(-3, 0))$minimum
    .fit <- suppressMessages(.nlmixr(
      .mod |> ini(tlam = log(0.5), eta.lam ~ fix(0.2)),
      .d,
      est = "focei",
      control = foceiControl(print = 0)
    ))
    expect_equal(unname(.fit$theta[["tlam"]]), .oFix, tolerance = 0.02)
    expect_gt(sd(.fit$eta$eta.lam), 0.2)
  })

  test_that("focei/agq keep the frailty variance of a Weibull repeated-event tte() model (#1465)", {
    skip_on_cran()
    .mod <- function() {
      ini({
        tlam <- log(0.3)
        tk <- log(1.4)
        eta.lam ~ 0.2
      })
      model({
        lam <- exp(tlam + eta.lam)
        k <- exp(tk)
        haz <- k * lam * (lam * time)^(k - 1)
        ev ~ tte(haz)
      })
    }
    .d <- .tteFrailtyData(.mod())
    .ev <- split(.d$TIME[.d$DV == 1], factor(.d$ID[.d$DV == 1], levels = unique(.d$ID)))
    .end <- tapply(.d$TIME, .d$ID, max)
    .ll <- function(p, i, eta) {
      .t <- .ev[[i]]
      .k <- exp(p[2])
      .lam <- p[1] + eta
      length(.t) * (p[2] + .k * .lam) + (.k - 1) * sum(log(.t)) - exp(.k * .lam) * .end[i]^.k
    }
    .o <- stats::optim(c(log(0.5), log(1.2), log(0.1)), .tteMarginalNll(length(.end), .ll), method = "BFGS")$par
    .start <- .mod |> ini(tlam = log(0.5), tk = log(1.2), eta.lam ~ 0.1)
    for (.est in c("focei", "agq")) {
      .fit <- suppressMessages(.nlmixr(.start, .d, est = .est, control = list(print = 0)))
      expect_equal(unname(.fit$theta[c("tlam", "tk")]), .o[1:2], tolerance = 0.02, label = .est)
      expect_equal(unname(.fit$omega[1, 1]), exp(.o[3]), tolerance = 0.1, label = .est)
    }
  })
})
