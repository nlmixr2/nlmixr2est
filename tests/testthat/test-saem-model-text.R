nmTest({
  test_that("saem model text builds from a compressed ui (#1236)", {
    one.cmt <- function() {
      ini({
        tka <- 0.45
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot
        cp <- central / exp(tv)
        cp ~ add(add.sd)
      })
    }
    .ui <- one.cmt()
    expect_false(is.environment(.ui))
    .txt <- rxUiGet.saemModel(list(.ui))
    expect_type(.txt, "character")
    expect_false(any(grepl("rx_cx_", .txt, fixed = TRUE)))
    expect_true(is.character(rxUiGet.saemModel(list(rxode2::rxUiDecompress(.ui)))))
  })
})
