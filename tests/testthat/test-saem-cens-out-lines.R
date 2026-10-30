nmTest({
  # .saemCensOutLines() reads the censoring flag a general-likelihood fit stages
  # on its ui environment.  rxUiGet.saemModel() is also called with a compressed
  # ui, which is not an environment and carries no flag.
  .mkUi <- function() {
    rxode2::rxode2(function() {
      ini({
        tka <- 0.45
        tcl <- 1
        tv <- 3.45
        eta.ka ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        cp <- central / v
        cp ~ add(add.sd)
      })
    })
  }

  test_that(".saemCensOutLines() is empty for a compressed ui", {
    .ui <- .mkUi()
    expect_false(is.environment(.ui))
    expect_identical(.saemCensOutLines(.ui, new.env()), character(0))
  })

  test_that(".saemCensOutLines() is empty for a ui environment with no flag staged", {
    .ui <- rxode2::rxUiDecompress(.mkUi())
    expect_true(is.environment(.ui))
    expect_false(exists(".saemCensOut", envir = .ui, inherits = FALSE))
    expect_identical(.saemCensOutLines(.ui, new.env()), character(0))
  })

  test_that(".saemCensOutLines() is empty when the model has no censoring outputs", {
    .ui <- rxode2::rxUiDecompress(.mkUi())
    for (.flag in 1L:2L) {
      assign(".saemCensOut", .flag, envir = .ui)
      expect_identical(.saemCensOutLines(.ui, new.env()), character(0), info = .flag)
      rm(list = ".saemCensOut", envir = .ui)
    }
  })

  test_that("rxUiGet.saemModel() builds the model text from a compressed ui", {
    .txt <- rxUiGet.saemModel(list(.mkUi()))
    expect_type(.txt, "character")
    expect_length(.txt, 1L)
    expect_match(.txt, "rx_pred_=", fixed = TRUE)
  })
})
