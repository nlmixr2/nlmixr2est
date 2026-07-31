test_that("vpcSimExpand ignores unknown extra columns (#830)", {
  sim <- data.frame(id = c(1, 1, 2), time = c(0, 1, 0), sim = c(1, 2, 3))
  obs <- data.frame(ID = c(1, 2), TIME = c(0, 0), DV = c(1, 3), WT = c(70, 80))
  # a column in neither frame warns and returns the simulation unchanged
  expect_warning(out <- vpcSimExpand(NULL, sim, "NOSUCHCOL", obs), "NOSUCHCOL")
  expect_identical(out, sim)
  # a column already in the simulation is left alone (no merge, no warning)
  expect_identical(vpcSimExpand(NULL, sim, "time", obs), sim)
  expect_identical(vpcSimExpand(NULL, sim, NULL, obs), sim)
})
