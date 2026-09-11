test_that("the subject information kernel differentiates the FOCEI objective", {
  y <- c(1.5, 1.1, 0.7)
  timeFactor <- c(1, 0.8, 0.5)
  mode <- function(p) {
    uniroot(function(eta) {
      z <- y/(exp(p[1]+eta)*timeFactor)
      sum(1-z*(z-1)/p[2]^2)+eta/p[3]
    }, c(-10, 10), tol = 1e-12)$root
  }
  objective <- function(p) {
    eta <- mode(p)
    f <- exp(p[1]+eta)*timeFactor
    variance <- (p[2]*f)^2
    sum((y-f)^2/variance+log(variance))+eta^2/p[3]+log(p[3])+
      log(length(y)*(1/p[2]^2+2)+1/p[3])
  }
  for (p in list(c(0.2, 0.2, 0.2), c(1, 0.5, 0.3))) {
    eta <- mode(p)
    f <- exp(p[1]+eta)*timeFactor
    variance <- (p[2]*f)^2
    a <- cbind(f, 0)
    A <- array(0, c(3, 2, 2)); A[, 1, 1] <- f
    third <- array(0, c(3, 1, 4)); third[, 1, 1] <- f
    aR <- cbind(2*variance, 2*p[2]*f^2)
    AR <- array(0, c(3, 2, 2))
    AR[, 1, 1] <- 4*variance
    AR[, 1, 2] <- AR[, 2, 1] <- 4*p[2]*f^2
    AR[, 2, 2] <- 2*f^2
    thirdR <- array(0, c(3, 1, 4))
    thirdR[, 1, 1] <- 8*variance
    thirdR[, 1, 2] <- thirdR[, 1, 3] <- 8*p[2]*f^2
    thirdR[, 1, 4] <- 4*f^2
    information <- foceiSubjectRFR_(a, A, third, aR, AR, thirdR,
      matrix(0, 3, 0), matrix(0, 3, 0), integer(), numeric(), f, y, variance,
      eta, matrix(1/p[3]), array(-1/p[3]^2, c(1, 1, 1)),
      array(2/p[3]^3, c(1, 1, 1)), matrix(-1/p[3]^2),
      1L, 2L, 2L, 1L, c(1L, 2L))
    reference <- numDeriv::hessian(objective, p)
    expect_equal(2*information, reference, tolerance = 1e-5)
    expect_equal(information, t(information), tolerance = 1e-12)
  }
})
