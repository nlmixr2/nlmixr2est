nmTest({
  # #1287: a fix()ed theta without an eta (a phi0 column) was sampled per
  # subject like a pseudo random effect, so the same model fit differently
  # depending on how many fixed terms it carried.
  .d <- nlmixr2data::theo_sd
  .d$lwt <- log(.d$WT / 70)
  .d$SW <- as.numeric(.d$TIME >= 6) # a time-varying flag

  .base <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      wt.cl <- fix(0.75)
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + wt.cl * lwt + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  # the same model: the flag effect is fixed at 0
  .flag0 <- .base |>
    model(cl <- exp(tcl + wt.cl * lwt + b.cl * SW + eta.cl)) |>
    ini(b.cl = fix(0), append = "wt.cl")

  .ctl <- function(literalFix = FALSE) {
    saemControl(
      nBurn = 10,
      nEm = 10,
      seed = 42,
      print = 0,
      covMethod = "",
      calcTables = FALSE,
      literalFix = literalFix
    )
  }
  .est <- function(fit) c(fit$theta[c("tka", "tcl", "tv", "add.sd")], diag(fit$omega))

  .fitBase <- suppressMessages(nlmixr2(.base, .d, est = "saem", control = .ctl()))
  .fitFlag0 <- suppressMessages(nlmixr2(.flag0, .d, est = "saem", control = .ctl()))
  .fitLit <- suppressMessages(nlmixr2(.base, .d, est = "saem", control = .ctl(TRUE)))

  test_that("fixed phi0 thetas stay at their value for every subject", {
    .phiM <- .fitFlag0$phiM # subject x chain x iteration x phi
    .nphi <- dim(.phiM)[4]
    # wt.cl and b.cl are the last two phi columns
    for (.j in c(.nphi - 1L, .nphi)) {
      expect_equal(unique(as.vector(.phiM[, , , .j])), c(0.75, 0)[.j - .nphi + 2L])
    }
    expect_equal(.fitFlag0$theta[["b.cl"]], 0)
    expect_equal(.fitFlag0$theta[["wt.cl"]], 0.75)
  })

  test_that("a fixed-at-0 covariate term gives the model without it", {
    expect_equal(.est(.fitFlag0), .est(.fitBase))
  })

  test_that("held fixed thetas fit the same as literally substituted ones", {
    expect_equal(.est(.fitLit), .est(.fitBase))
  })
})
