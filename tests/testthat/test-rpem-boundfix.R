# Bounded and fixed parameter support for est="rpem".
#
# Bounded: RPEM estimates on the unconstrained scale, so it sets the "unbounded" S3
# attribute (like saem/fsaem) and the shared .preProcessBoundedTransform hook rewrites a
# bounded theta to an internal unconstrained parameter (expit/exp) with a Jacobian-
# corrected back-transform.
#
# Fixed: a fix()ed typical value (mu-referenced theta) or residual parameter is held at
# its ini value across the M-step (the conjugate mu / residual updates would move it).

.bfData <- function(seed = 9L, nsub = 40L) {
  struct <- rxode2::rxode2({ ka <- exp(tka + eka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt() })
  set.seed(seed); obsT <- seq(0.5, 24, by = 2)
  do.call(rbind, lapply(seq_len(nsub), function(i) {
    eka <- rnorm(1, 0, sqrt(0.3))
    ev <- rxode2::et(amt = 100, cmt = "depot"); ev <- rxode2::et(ev, obsT)
    s <- rxode2::rxSolve(struct, params = c(tka = 0.45, tcl = 1, tv = 3.45, eka = eka),
                         events = ev, returnType = "data.frame", addDosing = FALSE)
    d <- as.data.frame(ev); d$id <- i; o <- d$evid == 0; d$DV <- 0
    d$DV[o] <- s$cp + rnorm(nrow(s), 0, 0.1); d
  }))
}

test_that("est=rpem declares the unbounded attribute (bounded-transform enabled)", {
  expect_true(isTRUE(attr(nlmixr2est:::nlmixr2Est.rpem, "unbounded")))
})

test_that("RPEM holds a fix()ed mu-referenced typical value and a fix()ed residual", {
  skip_on_cran()

  dat <- .bfData()
  # fixed typical value (with an eta): tka must stay at 0.30, only its omega is estimated
  modTka <- function() {
    ini({ tka <- fix(0.30); tcl <- fix(1.0); tv <- fix(3.45); add.sd <- 0.2; eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui <- rxode2::rxUiDecompress(rxode2::rxode2(modTka))
  cl <- .rpemClassify(ui)
  expect_true(cl$muFixFull[1])   # eta.ka's typical value is fixed
  rf <- .rpemFit(ui, dat, rpemControl(nGauss = 300L, nMH = 50000L, mhBurn = 5000L, niter = 25L,
                                      collect = 10L, seed = 1L, cores = 4L))
  expect_equal(unname(rf$mu["tka"]), 0.30, tolerance = 1e-6)   # held exactly
  expect_gt(rf$omega[1], 0)                                    # omega still estimated

  # fixed residual: add.sd must stay at 0.25
  modSd <- function() {
    ini({ tka <- 0.3; tcl <- fix(1.0); tv <- fix(3.45); add.sd <- fix(0.25); eta.ka ~ 0.6 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv); cp <- linCmt(); cp ~ add(add.sd) })
  }
  ui2 <- rxode2::rxUiDecompress(rxode2::rxode2(modSd))
  cl2 <- .rpemClassify(ui2)
  expect_true(cl2$addSdFix)
  rf2 <- .rpemFit(ui2, dat, rpemControl(nGauss = 300L, nMH = 50000L, mhBurn = 5000L, niter = 25L,
                                        collect = 10L, seed = 1L, cores = 4L))
  expect_equal(rf2$addSd, 0.25, tolerance = 1e-6)             # held exactly
  expect_equal(unname(rf2$mu["tka"]), 0.45, tolerance = 0.08) # tka still recovers
})

test_that("RPEM fits a bounded structural (likelihood) parameter within its bounds", {
  skip_on_cran()

  # Weibull TTE with a bounded shape via the full nlmixr2 path (bounded transform applies)
  mkWei <- function(seed = 1L, n = 200L, meanT = 40, k = 1.5) {
    set.seed(seed)
    do.call(rbind, lapply(seq_len(n), function(i) {
      lami <- meanT * exp(rnorm(1, 0, sqrt(0.15)))
      data.frame(ID = i, TIME = lami * (-log(runif(1)))^(1 / k), DV = 1, EVID = 0, CMT = 1)
    }))
  }
  weiB <- function() {
    ini({ tlam <- log(30); shape <- c(0.5, 1.0, 3.0); eta.lam ~ 0.2 })
    model({ lam <- exp(tlam + eta.lam)
            ll(dv) ~ log(shape) + (shape - 1) * log(time) - shape * log(lam) - (time / lam)^shape })
  }
  f <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(weiB, mkWei(1L), est = "rpem",
    control = rpemControl(nGauss = 400L, nMH = 50000L, mhBurn = 5000L, niter = 30L,
                          collect = 12L, seed = 1L, cores = 4L))))
  expect_s3_class(f, "nlmixr2FitData")
  .shape <- f$parFixedDf["shape", "Estimate"]
  expect_gte(.shape, 0.5); expect_lte(.shape, 3.0)             # respects the declared box
  # recovery band is wide: the mu M-step draws from the global threefry stream, so the
  # in-suite estimate drifts a little with test order; the band still confirms the shape
  # recovers in the right ballpark (true 1.5) within its declared bounds.
  expect_gt(.shape, 0.9); expect_lt(.shape, 2.3)
})
test_that("RPEM recovers a bounded mu-referenced typical value (centered-eta demotion)", {
  skip_on_cran()

  # The bounded transform rewrites tcl (bounded, mu-referenced by eta.cl) into a
  # structural theta rxBoundedTr.tcl + a model-computed tcl, so eta.cl becomes a
  # centered eta whose omega is still estimated by the conjugate M-step.
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45
      eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  # classification after the transform: eta.cl is centered, its theta is structural
  ui <- rxode2::rxUiDecompress(rxode2::rxode2(one.cmt))
  res <- .preProcessBoundedTransform(ui, "rpem", nlmixr2data::theo_sd, rpemControl())
  cl <- .rpemClassify(res$ui)
  expect_equal(cl$muRef, c(TRUE, FALSE, TRUE))
  expect_true("rxBoundedTr.tcl" %in% cl$thetaNames[cl$structIdx + 1L])

  f <- suppressMessages(suppressWarnings(nlmixr2est::nlmixr2(one.cmt, nlmixr2data::theo_sd,
    est = "rpem",
    control = rpemControl(nGauss = 300L, nMH = 50000L, mhBurn = 5000L, niter = 30L,
                          collect = 12L, seed = 7L, cores = 4L))))
  expect_s3_class(f, "nlmixr2FitData")
  .pf <- f$parFixedDf
  # tcl reported back-transformed inside its declared box; theo_sd reference ~ 1.01
  expect_lte(.pf["tcl", "Estimate"], log(100))
  expect_equal(unname(.pf["tcl", "Estimate"]), 1.01, tolerance = 0.15)
  expect_equal(unname(.pf["tka", "Estimate"]), 0.45, tolerance = 0.35)
  expect_equal(unname(.pf["tv", "Estimate"]), 3.45, tolerance = 0.1)
  .om <- diag(f$omega)
  expect_true(all(is.finite(.om) & .om > 0))       # centered eta.cl omega estimated
})
