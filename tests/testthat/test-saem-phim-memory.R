nmTest({
  test_that("saem keeps phiM in memory, not in a temp file (#1276)", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tcl <- 1
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
    .phiFiles <- function() {
      list.files(rxode2::rxTempDir(), pattern = "^phi-.*[.]phi$", recursive = TRUE)
    }
    .before <- .phiFiles()
    fit <- .nlmixr(
      one.cmt,
      theo_sd,
      "saem",
      saemControl(nBurn = 5, nEm = 5, print = 0, logLik = TRUE)
    )
    expect_equal(.phiFiles(), .before)
    expect_null(attr(fit$saem, "saem.cfg")$phiMFile)
    expect_null(fit$saem$phiM)
    .phiM <- fit$phiM
    expect_equal(dim(.phiM), c(12L, 3L, 10L, 3L))
    # full precision, not arma's 4-significant-digit text output
    expect_false(all(.phiM == signif(.phiM, 4)))
    expect_true(is.finite(fit$objf))
  })
})
