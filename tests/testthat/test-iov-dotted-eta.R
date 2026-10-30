.iovDottedModel <- function(occ, err = "cp ~ add(add.sd)", ini = "") {
  eval(str2lang(sprintf(
    "function() {
    ini({
      tka <- log(1.5); tcl <- log(2.7); tv <- log(30)
      eta.cl ~ 0.1; eta.v ~ 0.1
      %s ~ 0.1 | OCC
      add.sd <- 0.7
      %s
    })
    model({
      ka <- exp(tka + %s); cl <- exp(tcl + eta.cl); v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      %s
    })
  }",
    occ,
    ini,
    occ,
    err
  )))
}

.iovDottedData <- function() {
  d <- nlmixr2data::theo_md
  d$OCC <- ave(d$EVID != 0, d$ID, FUN = cumsum)
  d
}

.iovDottedFit <- function(model, est) {
  .ctl <- switch(
    est,
    focei = foceiControl(maxOuterIterations = 2, covMethod = "", print = 0),
    saem = saemControl(print = 0, nBurn = 5, nEm = 5)
  )
  suppressMessages(suppressWarnings(nlmixr2(model, .iovDottedData(), est = est, control = .ctl)))
}

test_that("occasion etas spelled like saem's temporary etas fit like any other name (#1233)", {
  skip_on_cran()
  # the occasion columns are named after the eta, so compare values only
  .vals <- function(fit) {
    list(
      theta = unname(fit$theta),
      omega = lapply(fit$omega, unname),
      ranef = fit$ranef,
      shrink = lapply(fit$shrink, unname),
      iov = lapply(fit$iov, function(d) unname(as.data.frame(d)))
    )
  }
  for (.est in c("focei", "saem")) {
    .ref <- .vals(.iovDottedFit(.iovDottedModel("iovKaOcc"), .est))
    for (.occ in c("eta.ka.occ", "eta.kaocc")) {
      .fit <- .iovDottedFit(.iovDottedModel(.occ), .est)
      expect_true("OCC" %in% names(.fit$omega), label = paste(.est, .occ))
      expect_equal(.vals(.fit), .ref, tolerance = 1e-6, label = paste(.est, .occ))
    }
  }
})

test_that("saem drops its temporary etas and keeps a dotted occasion eta (#1233)", {
  skip_on_cran()
  .fit <- .iovDottedFit(
    .iovDottedModel("eta.ka.occ", "cp ~ add(add.sd) + boxCox(lam) + dnorm()", "lam <- 0.5"),
    "saem"
  )
  .tmp <- c("rx.eta.add.sd", "rx.eta.lam")
  expect_false(any(.tmp %in% names(.fit$ranef)))
  expect_false(any(.tmp %in% names(.fit$shrink$id)))
  expect_false(any(.tmp %in% names(.fit)))
  expect_true("OCC" %in% names(.fit$omega))
  expect_true("eta.ka.occ" %in% names(.fit$iov$OCC))
})

test_that("no pseudoEta transforms records no temporary etas (#1233)", {
  expect_identical(.saemPseudoEtaRecord(NULL), character(0))
  expect_identical(.saemPseudoEtaRecord(list(list(name = "ka"))), character(0))
  expect_identical(
    .saemPseudoEtaRecord(list(list(name = "add.sd", pseudoEta = TRUE), list(name = "ka"))),
    "rx.eta.add.sd"
  )
})
