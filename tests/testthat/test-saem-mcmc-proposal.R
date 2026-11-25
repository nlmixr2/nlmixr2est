nmTest({
  test_that("saemControl(mcmcProposal=) is validated and survives a rebuild", {
    expect_equal(saemControl()$mcmc$proposal, "diag")
    .ctl <- saemControl(mcmcProposal = "omega")
    expect_equal(.ctl$mcmc$proposal, "omega")
    expect_equal(do.call(saemControl, .ctl)$mcmc$proposal, "omega")
    expect_error(saemControl(mcmcProposal = "full"))
    .d <- rxUiDeparse(.ctl, "ctl")
    expect_equal(eval(.d[[3]])$mcmc$proposal, "omega")
    expect_false(any(grepl("mcmcProposal", deparse(rxUiDeparse(saemControl(), "ctl")))))
  })

  test_that("full-omega random walk only differs when omega is correlated", {
    skip_on_cran()
    .mod <- function(cor = TRUE) {
      .f <- function() {
        ini({
          tka <- 0.45
          tcl <- 1
          tv <- 3.45
          eta.ka ~ 0.6
          eta.cl + eta.v ~ c(0.3, 0.1, 0.1)
          add.sd <- 0.7
        })
        model({
          ka <- exp(tka + eta.ka)
          cl <- exp(tcl + eta.cl)
          v <- exp(tv + eta.v)
          linCmt() ~ add(add.sd)
        })
      }
      if (cor) {
        return(.f)
      }
      .f |> ini(eta.cl ~ 0.3, eta.v ~ 0.1)
    }
    .fit <- function(m, p) {
      suppressMessages(suppressWarnings(nlmixr2(
        m,
        theo_sd,
        est = "saem",
        control = saemControl(
          print = 0, nBurn = 30, nEm = 30, seed = 42L, covMethod = "",
          calcTables = FALSE, mcmcProposal = p
        )
      )))
    }
    # a diagonal omega has a diagonal Cholesky factor: identical walk
    .m0 <- .mod(FALSE)
    expect_equal(.fit(.m0, "diag")$objective, .fit(.m0, "omega")$objective)
    .m1 <- .mod(TRUE)
    .d <- .fit(.m1, "diag")
    .o <- .fit(.m1, "omega")
    expect_false(isTRUE(all.equal(.d$omega, .o$omega)))
    expect_true(all(is.finite(.o$mcmcAccept[, "rw"])))
  })
})
