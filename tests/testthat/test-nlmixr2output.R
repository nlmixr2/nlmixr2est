test_that(".updateParFixedGetEtaRow returns correct values", {
  envPrep <- new.env()
  expect_equal(
    .updateParFixedGetEtaRow(
      .eta = "iivemax",
      .env = envPrep,
      .ome = matrix(25, nrow = 1, dimnames = list("iivemax", "iivemax")),
      .omegaFix = c(iivemax = FALSE),
      .muRefCurEval = data.frame(parameter = "iivemax", curEval = "", low = NA_real_, hi = NA_real_),
      .sigdig = 3L
    ),
    data.frame(ch = "5.00", v = 25)
  )
  expect_false(envPrep$.cvOnly)

  envPrep <- new.env()
  expect_equal(
    .updateParFixedGetEtaRow(
      .eta = "iivemax",
      .env = envPrep,
      .ome = matrix(0.4, nrow = 1, dimnames = list("iivemax", "iivemax")),
      .omegaFix = c(iivemax = FALSE),
      .muRefCurEval = data.frame(parameter = "iivemax", curEval = "exp", low = NA_real_, hi = NA_real_),
      .sigdig = 3L
    ),
    data.frame(ch = "70.1", v = sqrt(exp(0.4) - 1) * 100)
  )
  expect_false(envPrep$.sdOnly)
})
