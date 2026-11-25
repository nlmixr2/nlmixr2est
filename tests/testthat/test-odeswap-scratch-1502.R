nmTest({
  # #1502: two per-thread buffers were thread_local std::vectors, whose
  # destructors MinGW-w64 runs at thread exit after the emulated thread-local
  # storage is gone: freeing them corrupted the heap, and R stopped when an
  # OpenMP worker exited at the next multi-threaded fit (Windows only).
  test_that("a fast ll() fit can be repeated on two cores in one process (#1502)", {
    skip_on_cran()
    m <- function() {
      ini({
        tka <- log(1.2); tcl <- log(2); tv <- log(20)
        eta.cl ~ 0.09; eta.v ~ 0.04
        lsd <- log(0.15)
      })
      model({
        ka <- exp(tka); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v); sd <- exp(lsd)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        ll(err) ~ -0.5 * log(2 * pi) - log(sd * cp) - 0.5 * ((DV - cp) / (sd * cp))^2
      })
    }
    set.seed(1498)
    d <- do.call(
      rbind,
      lapply(1:6, function(i) {
        x <- data.frame(
          ID = i,
          TIME = c(0, 0.5, 1, 2, 4, 8, 12, 24),
          EVID = c(1, rep(0, 7)),
          AMT = c(100, rep(0, 7)),
          CMT = c(1, rep(2, 7))
        )
        k <- 2 * exp(rnorm(1, 0, 0.3)) / (20 * exp(rnorm(1, 0, 0.2)))
        v <- 20 * exp(rnorm(1, 0, 0.2))
        cp <- 100 / v * 1.2 / (1.2 - k) * (exp(-k * x$TIME) - exp(-1.2 * x$TIME))
        x$DV <- ifelse(x$EVID == 0, cp * (1 + 0.15 * rnorm(8)), 0)
        x
      })
    )
    ctl <- foceiControl(fast = TRUE, print = 0L, covMethod = "", rxControl = rxode2::rxControl(cores = 2L))
    f1 <- suppressMessages(suppressWarnings(nlmixr2(m, d, "focei", ctl)))
    f2 <- suppressMessages(suppressWarnings(nlmixr2(m, d, "focei", ctl)))
    expect_true(is.finite(f1$objf))
    expect_equal(f2$objf, f1$objf)
  })
})
