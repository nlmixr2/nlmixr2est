nmTest({
  test_that("nls supports interp", {
    one.cmt <- function() {
      ini({
        tka <- fix(0.45)
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl) + wt
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    f <- one.cmt()

    expect_false(grepl("linear\\(wt\\)", rxode2::rxNorm(f$nlsRxModel$predOnly)))

    one.cmt <- function() {
      ini({
        tka <- fix(0.45)
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        linear(wt)
        ka <- exp(tka)
        cl <- exp(tcl) + wt
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    f <- one.cmt()

    expect_true(grepl("linear\\(wt\\)", rxode2::rxNorm(f$nlsRxModel$predOnly)))
  })

  test_that("nls all 1 issue", {
    pheno <- function() {
      ini({
        tcl <- log(1) # typical value of clearance
        tv <-  log(1)   # typical value of volume
        add.err <- 0.1    # residual variability
      })
      model({
        cl <- exp(tcl ) # individual value of clearance
        v <- exp(tv)    # individual value of volume
        ke <- cl / v            # elimination rate constant
        d/dt(A1) = - ke * A1    # model differential equation
        cp = A1 / v             # concentration in plasma
        cp ~ add(add.err)       # define error model
      })
    }

    expect_error(.nlmixr(pheno, nlmixr2data::pheno_sd, est = "nls", nlsControl(algorithm = "LM", print = 0L)), NA)
  })

  test_that("nls makes sense", {
    d <- nlmixr2data::theo_sd

    d <- d[d$AMT != 0 | d$DV != 0, ]

    one.cmt <- function() {
      ini({
        tka <- fix(0.45)
        tcl <- log(c(0, 2.7, 100))
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        linCmt() ~ add(add.sd)
      })
    }

    fit1 <- .nlmixr(one.cmt, d, est = "nls", list(print = 0L))

    expect_true(inherits(fit1, "nlmixr2.nls"))

    fit1 <- .nlmixr(one.cmt, d, est = "nls", nlsControl(solveType = "fun", print = 0L))

    Treated <- Puromycin[Puromycin$state == "treated", ]
    names(Treated) <- gsub("rate", "DV", gsub("conc", "time", names(Treated)))
    Treated$ID <- 1

    f <- function() {
      ini({
        Vm <- 200
        K <- 0.1
        prop.sd <- 0.1
      })
      model({
        pred <- (Vm * time)/(K + time)
        pred ~ prop(prop.sd)
      })
    }

    fit1 <- .nlmixr(f, Treated, est = "nls", control = nlsControl(algorithm = "default", print = 0L))

    expect_true(inherits(fit1, "nlmixr2.nls"))
  })

  test_that("nls fits a delay() model with its past() pre-history", {
    # y' = -k*delay(y, 1), with the history y = a before time 0
    dde <- function() {
      ini({
        tk <- log(0.3)
        ta <- log(2)
        add.sd <- 0.05
      })
      model({
        k <- exp(tk)
        a <- exp(ta)
        y(0) <- a
        d/dt(y) <- -k * delay(y, 1)
        past(y, 1) <- a
        y ~ add(add.sd)
      })
    }
    ui <- rxode2::rxode2(dde)
    .past <- "^(past\\([^)]*\\))=.*$"
    .lines <- strsplit(rxode2::rxNorm(suppressMessages(ui$nlsRxModel)$predOnly), "\n")[[1]]
    expect_equal(sub(.past, "\\1", grep(.past, .lines, value = TRUE)), "past(y,1)")
    .sens <- suppressMessages(ui$nlsSensModel)
    .lines <- strsplit(rxode2::rxNorm(.sens$predOnly), "\n")[[1]]
    expect_equal(sub(.past, "\\1", grep(.past, .lines, value = TRUE)), "past(y,1)")
    # the history depends on ta (THETA[2]) only, so only its sensitivity has one
    .lines <- strsplit(rxode2::rxNorm(.sens$thetaGrad), "\n")[[1]]
    expect_equal(
      sub(.past, "\\1", grep(.past, .lines, value = TRUE)),
      c("past(y,1)", "past(rx__sens_y_BY_THETA_2___,1)")
    )

    .dat <- rxode2::rxWithSeed(42, {
      .s <- rxode2::rxSolve(
        rxode2::rxode2("k=0.3\na=2\ny(0)<-a\nd/dt(y)<- -k*delay(y,1)\npast(y,1)<-a\n"),
        rxode2::et(seq(0.5, 5, by = 0.5)),
        atol = 1e-9,
        rtol = 1e-9
      )
      data.frame(
        ID = rep(1:4, each = nrow(.s)),
        TIME = rep(.s$time, 4),
        DV = rep(.s$y, 4) + stats::rnorm(4 * nrow(.s), 0, 0.05)
      )
    })
    .nlm <- .nlmixr(dde, .dat, est = "nlm", control = nlmControl(print = 0L, calcTables = FALSE))
    for (.st in c("grad", "fun")) {
      .fit <- .nlmixr(dde, .dat, est = "nls", control = nlsControl(print = 0L, solveType = .st))
      # the residuals nls minimized are those of the fitted model (its table);
      # without the past() lines they were those of another history
      expect_equal(unname(.fit$nls$fvec), .fit$IRES, tolerance = 1e-4)
      # least squares and the nlm maximum likelihood share the optimum
      expect_equal(.fit$theta[c("tk", "ta")], .nlm$theta[c("tk", "ta")], tolerance = 1e-3)
    }
  })

  test_that("nls solveType='fun' fits a model that uses lag()", {
    lagMod <- function() {
      ini({
        tka <- 0.45
        tcl <- 1.0
        tv <- 3.45
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka)
        cl <- exp(tcl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(central) <- ka * depot - cl / v * central
        c0 <- central / v
        cp <- 0.5 * c0 + 0.5 * lag(c0)
        cp ~ add(add.sd)
      })
    }
    # the objective-only model defines the lag()-referenced c0
    .lines <- strsplit(rxode2::rxNorm(suppressMessages(rxode2::rxode2(lagMod)$nlsRxModel)$predOnly), "\n")[[1]]
    expect_equal(grep("^c0=", .lines, value = TRUE), "c0=exp(-THETA[3])*central;")

    .m <- rxode2::rxode2(
      "ka=exp(tka)\ncl=exp(tcl)\nv=exp(tv)\nd/dt(depot)=-ka*depot\nd/dt(central)=ka*depot-cl/v*central\nc0=central/v\ncp=0.5*c0+0.5*lag(c0)\n"
    )
    .ev <- rxode2::et(amt = 320, cmt = "depot") |> rxode2::et(seq(0.5, 24, by = 1.5))
    .s <- rxode2::rxSolve(.m, .ev, params = c(tka = 0.6, tcl = 1.1, tv = 3.6), returnType = "data.frame")
    .dat <- rxode2::rxWithSeed(1234, {
      .d <- rbind(
        data.frame(ID = 1:4, TIME = 0, DV = NA, AMT = 320, EVID = 1),
        data.frame(
          ID = rep(1:4, each = nrow(.s)),
          TIME = rep(.s$time, 4),
          DV = rep(.s$cp, 4) + stats::rnorm(4 * nrow(.s), 0, 0.3),
          AMT = 0,
          EVID = 0
        )
      )
      .d[order(.d$ID, .d$TIME, -.d$EVID), ]
    })
    .fit <- .nlmixr(lagMod, .dat, est = "nls", control = nlsControl(print = 0L, solveType = "fun"))
    # the residuals nls minimized are those of the model at its estimates
    .ref <- rxode2::rxSolve(.m, .dat, params = .fit$theta[c("tka", "tcl", "tv")], returnType = "data.frame")
    expect_equal(unname(.fit$nls$fvec), .dat$DV[.dat$EVID == 0] - .ref$cp, tolerance = 1e-5)
  })
})
