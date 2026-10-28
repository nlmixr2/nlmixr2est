# a dopred() that stops, as a solve failing inside the calculation does
.saemAuxStopPred <- function(...) {
  stop("the prediction stopped", call. = FALSE)
}

.saemAuxModel <- function() {
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

nmTest({
  # calc.2LL() and calc.COV() show a progress bar; the error that stops them
  # reaches the caller, not the bar's generic abort message
  test_that("an error inside saem's quadrature and linearization is the one raised", {
    .ctl <- saemControl(nBurn = 10, nEm = 10, print = 0, seed = 1L, covMethod = "", calcTables = FALSE)
    .fit <- .nlmixr(.saemAuxModel, nlmixr2data::theo_sd, est = "saem", control = .ctl)
    .s <- .fit$env$saem
    attr(.s, "env") <- .fit$env
    attr(.s, "dopred") <- .saemAuxStopPred
    expect_error(
      suppressMessages(calc.2LL(.s, nnodes.gq = 3, nsd.gq = 1.6, .fit$phiM)),
      "the prediction stopped"
    )
    expect_error(suppressMessages(calc.COV(.s)), "the prediction stopped")
  })
})
