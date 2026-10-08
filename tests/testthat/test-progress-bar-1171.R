# Issue #1171: the covariance-step progress bar stopped short of 100% when an
# earlier bar had already reached 100% (rxode2's par_progress() latch), and an
# in-place bar left the next message on its own line.
nmTest({
  test_that("covariance progress bar reaches 100% and ends its line (#1171)", {
    .setProg <- getFromNamespace("setProgSupported", "rxode2")
    .getProg <- getFromNamespace("getProgSupported", "rxode2")
    .oldProg <- .getProg()
    on.exit(.setProg(.oldProg), add = TRUE)

    one.compartment <- function() {
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
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    tf <- tempfile(fileext = ".txt")
    con <- file(tf, "w")
    .unwind <- function() {
      if (sink.number(type = "message") > 0L) {
        sink(type = "message")
      }
      if (sink.number() > 0L) sink()
    }
    sink(con)
    sink(con, type = "message")
    on.exit(
      {
        .unwind()
        try(close(con), silent = TRUE)
      },
      add = TRUE
    )
    # attaching rxode2 resets the bar style when not interactive, so attach first;
    # then draw the in-place bar, and leave the 100% latch set by a finished bar
    suppressPackageStartupMessages(library(rxode2))
    .setProg(1L)
    rxode2::rxProgress(1L)
    rxode2::rxTick()
    rxode2::rxProgressStop()
    fit <- suppressWarnings(nlmixr2(
      one.compartment,
      nlmixr2data::theo_sd,
      est = "focei",
      control = foceiControl(print = 0, maxOuterIterations = 0)
    ))
    # testthat captures the fit's own messages, so mark where the next output lands
    cat("nextOutput\n")
    .unwind()
    close(con)
    # readLines() would also split at the bar's carriage returns
    lines <- strsplit(rawToChar(readBin(tf, "raw", file.size(tf))), "\n", fixed = TRUE)[[1]]

    expect_true(is.matrix(fit$cov))
    .i <- grep("calculating covariance matrix", lines, fixed = TRUE)
    expect_length(.i, 1L)
    .bar <- lines[.i + 1L]
    .last <- utils::tail(strsplit(.bar, "\r", fixed = TRUE)[[1]], 1L)
    expect_match(.last, "100%", fixed = TRUE)
    expect_false(grepl("nextOutput", .bar, fixed = TRUE))
  })
})
