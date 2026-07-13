# Mixtures with a combined (add + prop) residual: the residual has no closed form, so
# rpemMstepMix accumulates per-(subject, sample, component) MH visit counts and runs the
# guarded (add.sd^2, prop.sd^2) Newton over the mixture's stored cp/dv (the same
# optimizer the non-mixture combined M-step uses).  The residual is shared across
# components (one add.sd + one prop.sd).

test_that("RPEM classifies a mixture with combined residual", {
  mod <- function() {
    ini({ tka1 <- 0.2; tka2 <- 1.2; tcl <- fix(1.0); tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.1; prop.sd <- 0.15; eta.ka ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) + prop(prop.sd) })
  }
  cl <- .rpemClassify(rxode2::rxode2(mod))
  expect_equal(cl$errType, 2L)                    # combined
  expect_false(is.null(cl$mix))
})

test_that("RPEM recovers a mixture with combined (add + prop) residual", {
  skip_on_cran()

  sim <- rxode2::rxode2({ ka <- exp(tka + eka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(42); nsub <- 150L; obsT <- seq(0.5, 24, by = 2); w1 <- 0.6; tkaTrue <- c(0.0, 1.4)
  dat <- do.call(rbind, lapply(seq_len(nsub), function(i) {
    k <- if (stats::runif(1) < w1) 1L else 2L
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(sim, params = c(tka = tkaTrue[k], tcl = 1, tv = 3.45,
                                         eka = stats::rnorm(1, 0, sqrt(0.3))),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + stats::rnorm(nrow(s), 0, sqrt(0.05^2 + (0.1 * s$cp)^2))
    d
  }))
  mod <- function() {
    ini({ tka1 <- 0.2; tka2 <- 1.2; tcl <- fix(1.0); tv <- fix(3.45)
          p1 <- 0.5; add.sd <- 0.1; prop.sd <- 0.15; eta.ka ~ 0.3 })
    model({ ka <- mix(exp(tka1 + eta.ka), p1, exp(tka2 + eta.ka))
            cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) + prop(prop.sd) })
  }
  rf <- .rpemFit(rxode2::rxode2(mod), dat,
                 rpemControl(nGauss = 400L, nMH = 80000L, mhBurn = 8000L,
                             niter = 30L, collect = 10L, seed = 7L, cores = 4L))
  # components separated
  expect_lt(rf$mix$muK[1], 0.6)
  expect_gt(rf$mix$muK[2], 1.0)
  # both residual parameters recovered (true add 0.05, prop 0.10); sane bands because
  # cores>1 MH is not bit-reproducible and small residual params have tight rel. tol.
  expect_true(rf$addSd > 0.02 && rf$addSd < 0.09)
  expect_true(rf$propSd > 0.06 && rf$propSd < 0.14)
  expect_true(is.finite(rf$omega) && rf$omega > 0)
})
