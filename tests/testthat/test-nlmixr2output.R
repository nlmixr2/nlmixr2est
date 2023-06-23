test_that(".updateParFixedGetEtaRow returns correct values", {
  expect_equal(
    .updateParFixedGetEtaRow(
      .eta = "iivemax",
      .env = new.env(),
      .ome = matrix(25, nrow = 1, dimnames = list("iivemax", "iivemax")),
      .omegaFix = c(iivemax = FALSE),
      .muRefCurEval = data.frame(parameter = "iivemax", curEval = "", low = NA_real_, hi = NA_real_),
      .sigdig = 3L
    ),
    data.frame(.ch = "5.00", .v = 25)
  )
})
