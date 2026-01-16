nmTest({
  test_that("model piping", {

    KA1Lode <- function() {
      ini({
        ltlag <- log(0.2)
        lka  <- log(1.15)
        lcl  <- log(0.135)
        lv   <- log(8)
        prop.err <- 0.15
        add.err  <- 0.6
        eta.tlag ~ 0.5
        eta.ka ~ 0.5
        eta.cl ~ 0.1
        eta.v  ~ 0.1
      })
      model({
        tlag <- exp(ltlag + eta.tlag)
        ka <- exp(lka + eta.ka)
        cl <- exp(lcl + eta.cl)
        v  <- exp(lv + eta.v)
        d/dt(gut) <- -ka*gut
        d/dt(central) <- ka*gut - (cl/v)*central
        lag(gut) <- tlag
        cp <- central/v
        cp ~ prop(prop.err) + add(add.err)
      })
    }

    d <- nlmixr2data::warfarin |>
      dplyr::filter(dvid=="cp")

    f <- .nlmixr(KA1Lode, data = d, est = "saem", control = saemControlFast)

    # General piping model updates work
    suppressMessages(expect_error(
      fUpV <-
        f |>
        model(v <- exp(lv)),
      NA
    ))
    expect_equal(
      methods::functionBody(as.function(fUpV))[[3]][[2]][[5]],
      str2lang("v <- exp(lv)")
    )

    # piping model updates work with append
    suppressMessages(expect_error(
      fUpFoo <-
        f |>
        model(foo <- exp(lv), append = TRUE),
      NA
    ))
    expect_equal(
      methods::functionBody(as.function(fUpFoo))[[3]][[2]][[11]],
      str2lang("foo <- exp(lv)")
    )
  })
})
