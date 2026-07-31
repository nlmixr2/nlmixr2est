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

test_that("vpcSimExpand merges only the requested columns (#830)", {
  # vpcNameDataCmts needs a real fit; the merge logic does not
  local_mocked_bindings(vpcNameDataCmts = function(object, data) data)
  sim <- data.frame(id = c(1, 1, 2), nlmixrRowNums = c(1, 2, 3),
                    time = c(0, 1, 0), sim = c(1, 2, 3))
  obs <- data.frame(ID = c(1, 1, 2), time = c(0, 1, 0), DV = c(1, 2, 3),
                    WT = c(70, 70, 80))
  out <- vpcSimExpand(NULL, sim, "WT", obs)
  expect_true("WT" %in% names(out))
  expect_false("DV" %in% names(out))
  expect_false(any(grepl("\\.(x|y)$", names(out))))
  expect_equal(nrow(out), 3L)

  # requesting the id column itself must not crash or duplicate it, even
  # when its case differs between the simulation and the data
  simU <- data.frame(ID = c(1, 2), nlmixrRowNums = c(1, 2), sim = c(1, 2))
  obsL <- data.frame(id = c(1, 2), WT = c(70, 80))
  out2 <- vpcSimExpand(NULL, simU, c("id", "WT"), obsL)
  expect_equal(sum(tolower(names(out2)) == "id"), 1L)
  expect_true("WT" %in% names(out2))

  simL <- data.frame(id = c(1, 2), nlmixrRowNums = c(1, 2), sim = c(1, 2))
  obsU <- data.frame(ID = c(1, 2), WT = c(70, 80))
  out3 <- vpcSimExpand(NULL, simL, c("ID", "WT"), obsU)
  expect_equal(sum(tolower(names(out3)) == "id"), 1L)
  expect_true("WT" %in% names(out3))
})
