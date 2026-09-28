nmTest({
  # #1132: a residual variance below sqrt(eps) used to be REPLACED by 1, putting a
  # ~+16 cliff per observation in the FOCEi objective where a proportional-error
  # prediction crossed ~1.2e-3.  It is now floored, so the objective is continuous.
  test_that("focei objective is continuous where the variance crosses sqrt(eps)", {
    m <- function() {
      ini({
        lf <- log(0.0012)
        eta.f ~ 0.1
        prop.sd <- 0.1
      })
      model({
        ipred <- exp(lf + eta.f)
        ipred ~ prop(prop.sd)
      })
    }
    d <- data.frame(ID = 1, TIME = c(1, 2), DV = c(0.0012, 0.0011), EVID = 0, AMT = 0)
    sd <- 0.1
    omega <- 0.1
    # eta at which the prediction's variance (sd * f)^2 equals sqrt(eps)
    fCross <- sqrt(sqrt(.Machine$double.eps)) / sd
    etaCross <- log(fCross / 0.0012)

    objAt <- function(eta) {
      ctl <- foceiControl(print = 0, covMethod = "", calcTables = FALSE,
                          maxOuterIterations = 0, maxInnerIterations = 0,
                          etaMat = matrix(eta, 1))
      suppressWarnings(nlmixr2(m, d, "focei", ctl))$objf
    }
    # R-side FOCEi objective: df/deta = f, dR/deta = 2R
    objR <- function(eta) {
      f <- 0.0012 * exp(eta)
      r <- (sd * f)^2
      h <- 1 / omega + nrow(d) * (f^2 / r + 0.5 * (2 * r)^2 / r^2)
      sum(log(r) + (d$DV - f)^2 / r) + eta^2 / omega + log(omega) + log(h)
    }

    below <- objAt(etaCross - 1e-4)
    above <- objAt(etaCross + 1e-4)
    # differential pair straddling the threshold: the old code jumped ~32 here
    expect_lt(abs(above - below), 0.05)
    # and the unfloored side matches the closed form
    expect_equal(above, objR(etaCross + 1e-4), tolerance = 1e-4)
  })
})
