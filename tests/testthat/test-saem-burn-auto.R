test_that("saemControl(nBurn = \"auto\") validates, round-trips and deparses", {
  .c <- saemControl(nBurn = "auto", nBurnMin = 120, nEm = 50)
  expect_true(.c$mcmc$nBurnAuto)
  expect_equal(.c$mcmc$nBurnMin, 120L)
  expect_equal(.c$mcmc$niter, c(500, 50))
  expect_null(saemControl()$mcmc$nBurnAuto)
  expect_identical(do.call(saemControl, unclass(.c))$mcmc, .c$mcmc)
  for (.x in list(.c, saemControl(nBurn = "auto"), saemControl(nBurn = "auto", nBurnMax = 300L))) {
    .e <- rxode2::rxUiDeparse(.x, "ctl")
    expect_identical(eval(.e[[3]])$mcmc, .x$mcmc, info = deparse1(.e))
  }
  expect_error(saemControl(nBurn = "auto", nBurnMax = 50))
  expect_error(saemControl(nBurn = "auto", nBurnMin = 400, nBurnMax = 300))
  expect_error(saemControl(nBurn = "fast"))
})

test_that("nBurn = \"auto\" ends the burn-in between nBurnMin and nBurnMax", {
  skip_on_cran()
  one.cmt <- function() {
    ini({
      tka <- 0.45
      tcl <- log(c(0, 2.7, 100))
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .fit <- function(...) {
    suppressMessages(nlmixr2(
      one.cmt,
      nlmixr2data::theo_sd,
      "saem",
      control = saemControl(nEm = 30, print = 0, calcTables = FALSE, covMethod = "", ...)
    ))
  }
  .f <- .fit(nBurn = "auto", nBurnMin = 100, nBurnMax = 300)
  .note <- grep("^auto burn-in used [0-9]+ iterations$", .f$runInfo, value = TRUE)
  expect_length(.note, 1L)
  .used <- as.integer(gsub("[^0-9]", "", .note))
  expect_gte(.used, 100L)
  expect_lte(.used, 300L)
  # skipped schedule rows are dropped from every per-iteration record
  expect_equal(nrow(.f$parHist), .used + 30L)
  expect_equal(attr(class(.f$parHistData), "niter"), .used)
  expect_equal(dim(.f$env$phiM)[3], .used + 30L)
  expect_equal(nrow(.f$mcmcAccept), .used + 30L)
  expect_true(all(is.finite(.f$theta)))

  # nBurnMin == nBurnMax never skips, so it is the fixed burn-in exactly
  .a <- .fit(nBurn = "auto", nBurnMin = 100, nBurnMax = 100)
  .i <- .fit(nBurn = 100)
  expect_identical(.a$parHist, .i$parHist)
  expect_identical(.a$theta, .i$theta)
  expect_identical(.a$omega, .i$omega)
  expect_false(any(grepl("burn-in", .i$runInfo)))

  # omegas started far too wide still shrink at nBurnMax and are named
  one.cmt <- one.cmt |> rxode2::ini(eta.cl ~ 10, eta.v ~ 5)
  .m <- suppressMessages(.fit(nBurn = "auto", nBurnMin = 100, nBurnMax = 100))
  .w <- grep("^burn-in hit nBurnMax, still trending: ", .m$runInfo, value = TRUE)
  expect_length(.w, 1L)
  expect_match(.w, "V(eta.cl)", fixed = TRUE)
  expect_lt(nchar(.w), 75L)
})
