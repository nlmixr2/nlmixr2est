nmTest({
  # #1499: a prediction clipped at zero with an eta on the threshold gives the
  # individual objective a kink the inner search can stop at.
  kink <- function() {
    ini({
      tth <- 1
      addSd <- 0.2
      eta.th ~ 0.5
    })
    model({
      th <- tth + eta.th
      if (TIME - th < 0) {
        f <- 0
      } else {
        f <- TIME - th
      }
      f ~ add(addSd)
    })
  }
  smooth <- function() {
    ini({
      tth <- 1
      addSd <- 0.2
      eta.th ~ 0.5
    })
    model({
      th <- tth + eta.th
      f <- TIME - th
      f ~ add(addSd)
    })
  }
  set.seed(2)
  n <- 20
  d <- data.frame(ID = rep(1:n, each = 6), TIME = rep(seq(0.5, 3, length.out = 6), n))
  th <- rep(1 + rnorm(n, 0, 0.7), each = 6)
  d$DV <- pmax(d$TIME - th, 0) + rnorm(nrow(d), 0, 0.2)

  .fit <- function(mod, ...) {
    suppressMessages(nlmixr2(
      mod,
      d,
      "focei",
      control = foceiControl(print = 0, covMethod = "", calcTables = FALSE, ...)
    ))
  }
  # -2 x the individual objective, up to a constant
  .objI <- function(fit, clip) {
    .th <- fixef(fit)
    .om <- fit$omega[1, 1]
    function(eta, id) {
      s <- d[d$ID == id, ]
      f <- s$TIME - (.th[["tth"]] + eta)
      if (clip) {
        f <- pmax(f, 0)
      }
      sum((s$DV - f)^2) / .th[["addSd"]]^2 + eta^2 / .om
    }
  }

  test_that("ebeCheck flags EBEs with a descent direction", {
    skip_on_cran()
    f <- .fit(kink, ebeCheck = TRUE)
    es <- f$ebeStationarity
    expect_s3_class(es, "data.frame")
    expect_equal(nrow(es), n)
    expect_named(es, c("ID", "dObjf", "slope", "stationary"))
    obj <- .objI(f, TRUE)
    h <- 0.1 * sqrt(f$omega[1, 1])
    eta <- f$eta$eta.th
    probe <- vapply(
      seq_len(n),
      function(i) {
        max(0, obj(eta[i], i) - min(obj(eta[i] + h, i), obj(eta[i] - h, i)))
      },
      numeric(1)
    )
    expect_equal(es$dObjf, probe, tolerance = 1e-4)
    # every flagged subject has a lower point on a brute-force grid
    grid <- seq(-1, 1, length.out = 2001) * sqrt(f$omega[1, 1])
    gridDrop <- vapply(
      seq_len(n),
      function(i) {
        obj(eta[i], i) - min(vapply(eta[i] + grid, obj, numeric(1), id = i))
      },
      numeric(1)
    )
    expect_true(any(!es$stationary))
    expect_true(all(gridDrop[!es$stationary] > 0.01))
    expect_true(any(grepl("non-stationary EBEs", f$runInfo)))
  })

  test_that("ebeRestart moves flagged EBEs to a lower individual objective", {
    skip_on_cran()
    f0 <- .fit(kink, ebeCheck = TRUE)
    f1 <- .fit(kink, ebeRestart = TRUE)
    es <- f1$ebeStationarity
    expect_named(
      es,
      c("ID", "dObjf", "slope", "stationary", "restarted", "dObjfRestart")
    )
    expect_true(any(es$restarted))
    expect_true(all(es$restarted <= !es$stationary))
    # thetas are not re-estimated; only the restarted subjects' EBEs move
    expect_equal(fixef(f1), fixef(f0))
    expect_equal(f1$eta$eta.th[!es$restarted], f0$eta$eta.th[!es$restarted])
    obj <- .objI(f1, TRUE)
    i <- which(es$restarted)
    drop <- vapply(
      i,
      function(j) {
        obj(f0$eta$eta.th[j], j) - obj(f1$eta$eta.th[j], j)
      },
      numeric(1)
    )
    expect_equal(es$dObjfRestart[i], drop, tolerance = 1e-4)
    expect_true(all(drop > 0))
    # the restarted EBE is at the brute-force grid minimum
    grid <- seq(-1, 1, length.out = 4001) * sqrt(f1$omega[1, 1])
    for (j in i) {
      v <- vapply(f0$eta$eta.th[j] + grid, obj, numeric(1), id = j)
      expect_lt(obj(f1$eta$eta.th[j], j) - min(v), 1e-3)
    }
    # the reported objective is taken at the moved EBEs
    expect_equal(
      f1$objf - f0$objf,
      sum(f1$etaObf$OBJI - f0$etaObf$OBJI),
      tolerance = 1e-4
    )
    expect_true(any(grepl("re-solved from their probe", f1$runInfo)))
    expect_true(f1$control$ebeCheck)
  })

  test_that("ebeCheck flags no subject of a smooth model", {
    skip_on_cran()
    f <- .fit(smooth, ebeCheck = TRUE)
    expect_true(all(f$ebeStationarity$stationary))
    expect_false(any(grepl("non-stationary EBEs", f$runInfo)))
  })

  test_that("ebeCheck is off by default", {
    skip_on_cran()
    f <- .fit(smooth)
    expect_null(f$ebeStationarity)
  })

  test_that("ebeCheck options are validated", {
    expect_false(foceiControl()$ebeCheck)
    expect_equal(foceiControl()$ebeCheckStep, 0.1)
    expect_equal(foceiControl()$ebeCheckTol, 0.01)
    expect_error(foceiControl(ebeCheck = NA))
    expect_error(foceiControl(ebeCheckStep = 0))
    expect_error(foceiControl(ebeCheckTol = -1))
    expect_false(foceiControl()$ebeRestart)
    expect_true(foceiControl(ebeRestart = TRUE)$ebeCheck)
    expect_error(foceiControl(ebeRestart = NA))
  })
})
