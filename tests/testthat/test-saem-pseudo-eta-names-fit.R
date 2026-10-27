# An occasion eta whose name starts with "eta." (eta.ka.occ) is rx.eta.ka.occ.<k>
# in the IOV rewrite, which saem's temporary-eta clean-up took for one of its
# own and dropped, so .uiFinalizeIov() stopped ("incorrect number of
# dimensions").  The fit equals that of the same model with another spelling.

# theo_md, a new occasion at each dose; the occasion eta on ka is `occEta`
.pseudoNameMod <- function(occEta, extraIni = "", errLine = "cp ~ add(add.sd)") {
  eval(parse(
    text = sprintf(
      "function() {
        ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(30); eta.cl ~ 0.1; eta.v ~ 0.1; %s ~ 0.1 | OCC; add.sd <- 0.7; %s })
        model({
          ka <- exp(tka + %s)
          cl <- exp(tcl + eta.cl)
          v <- exp(tv + eta.v)
          d/dt(depot) <- -ka * depot
          d/dt(center) <- ka * depot - cl / v * center
          cp <- center / v
          %s
        })
      }",
      occEta,
      extraIni,
      occEta,
      errLine
    )
  ))
}

.pseudoNameData <- function() {
  .d <- nlmixr2data::theo_md
  .d$OCC <- ave(.d$EVID != 0, .d$ID, FUN = cumsum)
  .d
}

.pseudoNameFit <- function(occEta, est) {
  .ctl <- if (est == "focei") {
    foceiControl(print = 0, maxOuterIterations = 0L, covMethod = "")
  } else {
    saemControl(print = 0, nBurn = 20, nEm = 20, covMethod = "", seed = 42L)
  }
  suppressWarnings(suppressMessages(nlmixr2(.pseudoNameMod(occEta), .pseudoNameData(), est, control = .ctl)))
}

nmTest({
  test_that("an occasion eta named eta.<x> fits as any other spelling does", {
    for (.est in c("focei", "saem")) {
      .ref <- .pseudoNameFit("iovKaOcc", .est)
      for (.nm in c("eta.ka.occ", "eta.kaocc")) {
        .lab <- paste(.est, .nm)
        .f <- .pseudoNameFit(.nm, .est)
        expect_equal(.f$objf, .ref$objf, tolerance = 1e-10, label = .lab)
        expect_identical(names(.f$ranef), names(.ref$ranef), label = .lab)
        expect_equal(.f$ranef, .ref$ranef, tolerance = 1e-10, label = .lab)
        # the table, with the occasion eta's column renamed
        .t <- as.data.frame(.f)
        names(.t)[names(.t) == .nm] <- "iovKaOcc"
        expect_identical(names(.t), names(as.data.frame(.ref)), label = .lab)
        expect_equal(.t, as.data.frame(.ref), tolerance = 1e-10, label = .lab)
        # the occasion estimates
        .iov <- .f$iov$OCC
        names(.iov)[names(.iov) == .nm] <- "iovKaOcc"
        expect_equal(.iov, .ref$iov$OCC, tolerance = 1e-10, label = .lab)
      }
    }
  })

  test_that("saem drops its temporary etas and keeps an occasion eta named eta.<x>", {
    # boxCox(lam) under dnorm() is an eta-less likelihood theta: saem fits it
    # with a temporary eta, rx.eta.lam
    .m <- .pseudoNameMod("eta.ka.occ", "lam <- 0.5", "cp ~ add(add.sd) + boxCox(lam) + dnorm()")
    .d <- .pseudoNameData()
    .d <- .d[.d$EVID != 0 | .d$DV > 0, ]
    .f <- suppressWarnings(suppressMessages(nlmixr2(
      .m,
      .d,
      "saem",
      control = saemControl(print = 0, nBurn = 20, nEm = 20, covMethod = "", seed = 42L)
    )))
    expect_true(any(grepl("temporary eta for eta-less likelihood theta", .f$runInfo, fixed = TRUE)))
    expect_false(any(grepl(
      "^rx[.](eta|l)[.]|^rxBoundedTr",
      c(names(.f), names(.f$eta), rownames(.f$omega$id), names(fixef(.f)))
    )))
    expect_true("eta.ka.occ" %in% names(.f))
    expect_identical(names(.f$iov$OCC), c("ID", "OCC", "eta.ka.occ"))
    expect_identical(names(.f$ranef), c("ID", "eta.cl", "eta.v"))
  })
})
