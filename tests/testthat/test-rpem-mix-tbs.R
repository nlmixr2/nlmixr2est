# Mixtures with a transform-both-sides residual (additive on the boxCox/yeoJohnson
# transformed scale + a dynamic lambda): the transformed-scale variance profiles out and
# lambda is golden-sectioned over the mixture's accepted (subject, sample, component)
# states with the log-Jacobian (same profile as the non-mixture TBS M-step; per-obs
# transform code/bounds come from the E-step).  The residual is shared across components.

test_that("RPEM classifies a mixture with a TBS (boxCox) residual", {
  mod <- function() {
    ini({ tka1 <- 0.2; tka2 <- 1.2; tcl <- fix(1.0); tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.2; lambda <- 1; eta.ka ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) + boxCox(lambda) })
  }
  cl <- .rpemClassify(rxode2::rxode2(mod))
  expect_equal(cl$errType, 3L)                    # TBS
  expect_false(is.null(cl$mix))
})

test_that("RPEM recovers a mixture with a TBS (boxCox) residual", {
  skip_on_cran()

  sim <- rxode2::rxode2({ ka <- exp(tka + eka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(42); nsub <- 150L; obsT <- seq(0.5, 24, by = 2); w1 <- 0.6; tkaTrue <- c(0.0, 1.4)
  trueLam <- 0.5; trueAdd <- 0.1
  bcI <- function(y, l) if (l == 0) exp(y) else (l * y + 1)^(1 / l)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    k <- if (stats::runif(1) < w1) 1L else 2L
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(sim, params = c(tka = tkaTrue[k], tcl = 1, tv = 3.45,
                                         eka = stats::rnorm(1, 0, sqrt(0.3))),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    tcp <- ((pmax(s$cp, 1e-3))^trueLam - 1) / trueLam
    d$DV[o] <- pmax(bcI(tcp + stats::rnorm(nrow(s), 0, trueAdd), trueLam), 1e-3)
    d
  }))
  mod <- function() {
    ini({ tka1 <- 0.2; tka2 <- 1.2; tcl <- fix(1.0); tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.2; lambda <- 1; eta.ka ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) + boxCox(lambda) })
  }
  rf <- .rpemFit(rxode2::rxode2(mod), dat,
                 rpemControl(nGauss = 400L, nMH = 80000L, mhBurn = 8000L,
                             niter = 30L, collect = 10L, seed = 7L, cores = 4L))
  # components separated
  expect_lt(rf$mix$muK[1], 0.6)
  expect_gt(rf$mix$muK[2], 1.0)
  # transformed-scale add.sd ~ 0.1 and lambda ~ 0.5 recovered (sane bands: cores>1 MH not
  # bit-reproducible)
  expect_true(rf$addSd > 0.06 && rf$addSd < 0.15)
  expect_true(rf$lambda > 0.3 && rf$lambda < 0.75)
  expect_true(is.finite(rf$omega) && rf$omega > 0)
})
