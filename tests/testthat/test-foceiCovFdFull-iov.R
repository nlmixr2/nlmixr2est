# The full FD covariance of a correlated IOV block estimates ONE occasion Omega,
# shared by every occasion's same() copy (#1316).  Weekly via .slowBatches.

.fdIovData <- function(seed = 42, nid = 40, nocc = 3) {
  set.seed(seed)
  .tt <- c(0.5, 1, 2, 4, 8, 12)
  .L <- t(chol(matrix(c(0.08, 0.02, 0.02, 0.06), 2)))
  .eta <- cbind(rnorm(nid, 0, sqrt(0.1)), rnorm(nid, 0, sqrt(0.05)))
  .d <- do.call(
    rbind,
    lapply(seq_len(nid), function(i) {
      do.call(
        rbind,
        lapply(seq_len(nocc), function(o) {
          .z <- .L %*% rnorm(2)
          .cl <- exp(1 + .eta[i, 1] + .z[1])
          .v <- exp(3.45 + .eta[i, 2] + .z[2])
          .ka <- exp(0.45)
          .k <- .cl / .v
          .cp <- 320 * .ka / (.v * (.ka - .k)) * (exp(-.k * .tt) - exp(-.ka * .tt))
          .t0 <- (o - 1) * 48
          rbind(
            data.frame(ID = i, TIME = .t0, EVID = 101, AMT = 320, DV = NA, occ = o),
            data.frame(ID = i, TIME = .t0 + .tt, EVID = 0, AMT = 0, DV = .cp + rnorm(length(.tt), 0, 0.3), occ = o)
          )
        })
      )
    })
  )
  .d
}

.fdIovMod <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    add.sd <- 0.3
    eta.cl ~ 0.1
    eta.v ~ 0.05
    iov.cl + iov.v ~ c(0.08, 0.02, 0.06) | occ
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl + iov.cl)
    v <- exp(tv + eta.v + iov.v)
    d / dt(depot) <- -ka * depot
    d / dt(center) <- ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

test_that("the full FD covariance ties the occasion copies of an IOV block (#1316)", {
  skip_on_cran()
  .d <- .fdIovData()
  .fit <- suppressWarnings(suppressMessages(
    nlmixr2(.fdIovMod, .d, "focei", control = foceiControl(print = 0, calcTables = FALSE))
  ))
  expect_identical(.fit$covMethod, "r,s (full)")
  .nm <- c("tka", "tcl", "tv", "add.sd", "om.eta.cl", "om.eta.v", "om.iov.cl", "cov.iov.v.iov.cl", "om.iov.v")
  expect_identical(rownames(.fit$cov), .nm)
  # the occasion level's SEs now come from $cov
  expect_identical(.ranParTable(.fit)$method, rep("cov", 5))

  # each tied direction moves every occasion's copy: its curvature is the
  # objective's second difference when the whole occasion block moves, at the
  # stage's own step (R = H / 2)
  .om <- .fit$omega$occ
  .x0 <- c(.om[1, 1], .om[1, 2], .om[2, 2])
  .ofv <- function(x) {
    .ui <- suppressMessages(eval(bquote(rxode2::ini(.fit$finalUi, iov.cl + iov.v ~ c(.(x[1]), .(x[2]), .(x[3])) | occ))))
    suppressWarnings(suppressMessages(nlmixr2(
      .ui,
      .d,
      "focei",
      control = foceiControl(
        print = 0,
        calcTables = FALSE,
        maxOuterIterations = 0,
        covMethod = ""
      )
    )))$objf
  }
  .f0 <- .ofv(.x0)
  .h <- .fit$env$.fdFullH[7:9]
  .r <- vapply(
    1:3,
    function(k) {
      .e <- replace(numeric(3), k, .h[k])
      0.5 * (.ofv(.x0 + .e) - 2 * .f0 + .ofv(.x0 - .e)) / .h[k]^2
    },
    numeric(1)
  )
  expect_equal(unname(diag(.fit$env$.fdFullR)[7:9]), .r, tolerance = 0.03)
})
