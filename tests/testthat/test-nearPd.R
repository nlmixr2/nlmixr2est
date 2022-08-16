test_that("test nearPD with same functions as Matrix", {
  # https://github.com/cran/Matrix/blob/65c37738919156f6bbb682a1d8198d715a82a19f/tests/dpo-test.R#L69-L136
  # Testing nearPD() --- this is partly in  ../man/nearPD.Rd :
  
  pr <- matrix(c(1,     0.477, 0.644, 0.478, 0.651, 0.826,
                 0.477, 1,     0.516, 0.233, 0.682, 0.75,
                 0.644, 0.516, 1,     0.599, 0.581, 0.742,
                 0.478, 0.233, 0.599, 1,     0.741, 0.8,
                 0.651, 0.682, 0.581, 0.741, 1,     0.798,
                 0.826, 0.75,  0.742, 0.8,   0.798, 1),
               nrow = 6, ncol = 6)

  pr <- 0.5*(pr + t(pr))

  .npd <- function(x, conv.tol = 1e-7,...) {
    .ret <- as.matrix(Matrix::nearPD(x, conv.tol=conv.tol, ...)$mat)
    dimnames(.ret) <- NULL
    .ret
  }

  expect_equal(nmNearPD(pr, conv.tol = 1e-7), .npd(pr))
  expect_equal(nmNearPD(pr, conv.tol = 1e-7, doDykstra=FALSE), .npd(pr, doDykstra=FALSE))
})
