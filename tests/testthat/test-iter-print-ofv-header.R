# Regression tests for the unified iteration-time printer (src/scale.h):
# focei's Function Val. column and the periodic header's Key-legend suppression.
#
# The iteration trace is written from C and surfaces on the message stream,
# so it is captured with a dual output+message sink; do not wrap in
# suppressMessages() -- that muffles the trace itself.
nmTest({
  .captureIterTrace <- function(expr) {
    tf <- tempfile(fileext = ".txt")
    con <- file(tf, "w")
    .unwind <- function() {
      if (sink.number(type = "message") > 0L) sink(type = "message")
      if (sink.number() > 0L) sink()
    }
    sink(con)
    sink(con, type = "message")
    on.exit({
      .unwind()
      if (isOpen(con)) close(con)
    }, add = TRUE)
    force(expr)
    .unwind()
    close(con)
    on.exit()
    readLines(tf, warn = FALSE)
  }

  .wang2007PropFit <- function(estName, control) {
    f <- function() {
      ini({
        tke <- 0.5
        eta.ke ~ 0.04
        prop.sd <- sqrt(0.1)
      })
      model({
        ke <- tke * exp(eta.ke)
        ipre <- 10 * exp(-ke * t)
        ipre ~ prop(prop.sd)
      })
    }
    .d <- Wang2007
    .d$DV <- .d$Y
    # suppressWarnings (but never suppressMessages) so the trace is preserved.
    .captureIterTrace(
      suppressWarnings(nlmixr2(.nlmixr(f), .d, est = estName, control = control)))
  }

  test_that("focei shows the Function Val. objective column for every outer optimizer", {
    skip_on_cran()
    for (.opt in c("bobyqa", "lbfgsb3c")) {
      .out <- .wang2007PropFit(
        "focei",
        foceiControl(outerOpt = .opt,
                     maxOuterIterations = 15L, maxInnerIterations = 15L,
                     covMethod = "", calcTables = FALSE,
                     print = iterPrintControl(every = 1L, headerEvery = 3L)))
      .info <- paste0("outerOpt=", .opt)
      # the objective column header is present ...
      expect_true(any(grepl("|    #| Function Val. |", .out, fixed = TRUE)),
                  info = paste0(.info, ": missing 'Function Val.' header"))
      # ... and the # rows actually carry an objective value in that slot.
      expect_true(any(grepl("^\\|\\s*[0-9]+\\|\\s*[0-9.eE+-]+ \\|", .out)),
                  info = paste0(.info, ": # rows carry no objective value"))
    }
  })

  test_that("focei periodic header repeats column labels but not the Key legend", {
    skip_on_cran()
    .out <- .wang2007PropFit(
      "focei",
      foceiControl(outerOpt = "bobyqa",
                   maxOuterIterations = 15L, maxInnerIterations = 15L,
                   covMethod = "", calcTables = FALSE,
                   print = iterPrintControl(every = 1L, headerEvery = 3L)))
    .nKey    <- sum(grepl("^Key:", .out))
    .nHeader <- sum(grepl("^\\|\\s*#\\|", .out))
    # the legend is printed exactly once, at fit start ...
    expect_equal(.nKey, 1L)
    # ... while the column-label header re-emits several more times without it.
    expect_gt(.nHeader, 1L)
  })

  test_that("saem keeps no Function Val. column and re-emits its header without the Key", {
    skip_on_cran()
    one.cmt <- function() {
      ini({
        tka <- 0.45; tcl <- log(c(0, 2.7, 100)); tv <- 3.45
        add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }
    .out <- .captureIterTrace(
      suppressWarnings(nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "saem",
        control = saemControl(nBurn = 30L, nEm = 30L, nmc = 2L,
                              print = iterPrintControl(every = 2L, headerEvery = 3L)))))
    # saem has no per-iteration objective, so the objective column stays off.
    expect_false(any(grepl("Function Val.", .out, fixed = TRUE)))
    # the legend appears at most once (startup) ...
    expect_lte(sum(grepl("^Key:", .out)), 1L)
    # ... and the header is still re-emitted for readability.
    expect_gt(sum(grepl("^\\|\\s*#\\|", .out)), 1L)
  })
})
