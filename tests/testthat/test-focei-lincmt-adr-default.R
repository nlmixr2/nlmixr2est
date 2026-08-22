# rxode2's default linCmt() sensitivity method resolves to reverse-mode AD
# (linCmtSensType="auto" -> 31, "ADr").  nlmixr2est does not pin
# linCmtSensType anywhere, so a FOCEi fit must really take that path; the
# observable is which sensType codes linCmtB() computed a Jacobian with.
test_that("a linCmt() FOCEi fit uses rxode2's reverse-mode AD default", {
  skip_on_cran()
  skip_if_not(exists("linCmtBSensTypesSeen", envir = asNamespace("rxode2")),
              "rxode2 without linCmtBSensTypesSeen()")
  seen <- function() rxode2:::linCmtBSensTypesSeen(TRUE)
  one.cmt <- function() {
    ini({
      tka <- 0.45; tcl <- log(c(0, 2.7, 100)); tv <- 3.45
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
  ctl <- function(st) {
    foceiControl(print = 0, maxOuterIterations = 0L, covMethod = "",
                 calcTables = FALSE, rxControl = rxode2::rxControl(linCmtSensType = st))
  }
  fit <- function(st) {
    invisible(seen())
    f <- suppressWarnings(suppressMessages(
      nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "focei", control = ctl(st))))
    list(objf = f$objf, seen = seen())
  }
  auto <- fit("auto")
  expect_true(31L %in% auto$seen)
  expect_false(3L %in% auto$seen)
  fwd <- fit("AD")
  expect_true(3L %in% fwd$seen)
  expect_false(31L %in% fwd$seen)
  expect_equal(auto$objf, fwd$objf, tolerance = 1e-6)
})
