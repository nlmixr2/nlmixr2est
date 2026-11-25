nmTest({
  test_that("saemControl(rwProposal=) is validated and survives a rebuild", {
    expect_equal(saemControl()$mcmc$proposal, "diag")
    .ctl <- saemControl(rwProposal = "omega")
    expect_equal(.ctl$mcmc$proposal, "omega")
    expect_equal(do.call(saemControl, .ctl)$mcmc$proposal, "omega")
    expect_error(saemControl(rwProposal = "full"))
    .d <- rxUiDeparse(.ctl, "ctl")
    expect_equal(eval(.d[[3]])$mcmc$proposal, "omega")
    expect_false(any(grepl("rwProposal", deparse(rxUiDeparse(saemControl(), "ctl")))))
  })

  test_that("full-omega random walk only differs when omega is correlated", {
    skip_on_cran()
    .diag <- function() {
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
    .block <- function() {
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
    .fit <- function(m, p) {
      suppressMessages(suppressWarnings(nlmixr2(
        m,
        theo_sd,
        est = "saem",
        control = saemControl(
          print = 0, nBurn = 30, nEm = 30, seed = 42L, covMethod = "",
          calcTables = FALSE, rwProposal = p
        )
      )))
    }
    # a diagonal omega has a diagonal Cholesky factor: the same walk
    expect_equal(.fit(.diag, "diag")$omega, .fit(.diag, "omega")$omega)
    .d <- .fit(.block, "diag")
    .o <- .fit(.block, "omega")
    expect_false(isTRUE(all.equal(.d$omega, .o$omega)))
    expect_true(all(is.finite(.o$mcmcAccept[, "rw"])))
  })
})
