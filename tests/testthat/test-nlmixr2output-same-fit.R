# The model with eta.v ~ same() of eta.cl's block, with or without IOV
.sameFitMod <- function(iov) {
  eval(parse(
    text = sprintf(
      "function() {
        ini({
          tka <- 0.45; tcl <- 1; tv <- 3.45; add.sd <- c(0, 0.7)
          eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ same() %s
        })
        model({
          ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl %s); v <- exp(tv + eta.v)
          linCmt() ~ add(add.sd)
        })
      }",
      if (iov) "; iov.cl ~ 0.1 | occ" else "",
      if (iov) "+ iov.cl" else ""
    )
  ))
}

# The error message of `expr`, or NULL
.sameFitErr <- function(expr) {
  tryCatch(
    {
      expr
      NULL
    },
    error = conditionMessage
  )
}

test_that("a fit's model keeps a repeated (same()) block, so a refit keeps the tie", {
  skip_on_cran()
  # eta.v repeats eta.cl's block ("id:same:eta.cl"); the output rebuilds the
  # omega rows from the estimated matrix, which carries no such condition
  .d <- nlmixr2data::theo_sd
  .d$occ <- 1 + (.d$TIME >= 5)
  .ctl <- foceiControl(print = 0, covMethod = "")
  .refusal <- paste0(
    "'saem' does not honour a repeated (same()) omega block\n",
    "use fo, foi, foce, focei, focep, laplace, agq, emvi, fbvi or vae"
  )
  for (.iov in c(FALSE, TRUE)) {
    .fit <- suppressWarnings(suppressMessages(nlmixr2(.sameFitMod(.iov)(), .d, "focei", .ctl)))
    .idf <- .fit$ui$iniDf
    expect_identical(.idf$condition[.idf$name == "eta.v"], "id:same:eta.cl", label = .iov)
    expect_identical(.fit$ui$omegaSameMap, if (.iov) c(0L, 0L, 2L, 0L) else c(0L, 0L, 2L), label = .iov)
    .refit <- suppressWarnings(suppressMessages(nlmixr2(.fit$ui, .d, "focei", .ctl)))
    .om <- if (is.list(.refit$omega)) .refit$omega$id else .refit$omega
    expect_equal(.om["eta.v", "eta.v"], .om["eta.cl", "eta.cl"], tolerance = 1e-12, label = .iov)
    # a method that would drop the tie now sees it in a fit too
    expect_identical(.sameFitErr(suppressMessages(nlmixr2(.fit, est = "saem"))), .refusal, label = .iov)
    expect_identical(
      .sameFitErr(suppressMessages(nlmixr2(rxode2::ini(.fit, tka = 0.5), .d, "saem"))),
      .refusal,
      label = .iov
    )
  }
})
