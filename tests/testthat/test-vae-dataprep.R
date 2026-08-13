## Regression for the neonate-style data path: .vaeDataPrep must not crash when
## the data has neither an AMT nor an EVID column (dose-free observations). The
## former `ifelse(is.na(d$AMT) | d$AMT == 0, ...)` yielded a length-0 vector when
## d$AMT was NULL, giving "replacement has 0 rows, data has N". Fast: no fit.

nmTest({
  test_that(".vaeDataPrep handles data with no AMT and no EVID columns", {
    theo <- function() {
      ini({
        lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
        add.err <- 0.7
      })
      model({
        ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V
        cp ~ add(add.err)
      })
    }
    ui <- rxode2::assertRxUi(theo)

    d <- nlmixr2data::theo_sd
    d <- d[d$EVID == 0, , drop = FALSE]   # observations only
    d$AMT <- NULL
    d$EVID <- NULL
    expect_false("AMT" %in% names(d))
    expect_false("EVID" %in% names(d))

    prep <- .vaeDataPrep(ui, d)
    expect_equal(prep$N, length(unique(d$ID)))
    expect_equal(prep$Nobs, nrow(d))       # every row treated as an observation
    expect_true(all(vapply(prep$subj, function(s) all(s$ev$EVID == 0L), logical(1))))
  })

  ## residOptimize="twoStage" stage-2 eligibility.  Stage 2 pins the ODE states
  ## from stage 1, so a parameter belongs there exactly when the state trajectory
  ## cannot depend on it.  Fast: classification only, no fit.
  .llMod <- function() {
    ini({ lka <- 0.4; lcl <- -0.9; lb <- log(3); lsd <- log(1.2); eta.b ~ 0.1 })
    model({
      ka <- exp(lka); cl <- exp(lcl)
      d / dt(depot) <- -ka * depot
      d / dt(central) <- ka * depot - cl * central
      mu <- exp(lb + eta.b) * central
      sd <- exp(lsd)
      ll(cp) ~ -0.5 * log(2 * pi) - log(sd) - 0.5 * ((DV - mu) / sd)^2
    })
  }
  .gMod <- function() {
    ini({ tka <- 0.45; tcl <- 1; tv <- c(2, 3.45, 5); add.sd <- c(0, 0.7, 5)
      eta.ka ~ 0.6; eta.cl ~ 0.3 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
      d / dt(depot) <- -ka * depot
      d / dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd) })
  }

  test_that("the ODE-invariance scan separates structural from log-density thetas", {
    ui <- rxode2::assertRxUi(.llMod())
    ## lka/lcl reach a d/dt right-hand side; lb/lsd are read only by the density
    expect_equal(.vaeOdeFreeThetas(ui, c("lka", "lcl", "lb", "lsd")),
                 c(FALSE, FALSE, TRUE, TRUE))
  })

  test_that("a variable feeding d/dt keeps its dependency through the endpoint name", {
    ## The dangerous direction: `conc` is BOTH the endpoint variable and an input
    ## to d/dt.  Its `<-` definition must stay in the dependency map, or lv would
    ## be called ODE-free and stage 2 would optimize it against a frozen solve.
    .trapMod <- function() {
      ini({ lk <- -1; lv <- 1; lsd <- log(0.7); eta.b ~ 0.1 })
      model({
        k <- exp(lk); v <- exp(lv)
        conc <- central / v
        d / dt(central) <- -k * conc + eta.b * 0
        sd <- exp(lsd)
        ll(conc) ~ -log(sd) - 0.5 * ((DV - conc) / sd)^2
      })
    }
    ui <- suppressWarnings(rxode2::assertRxUi(.trapMod()))
    expect_equal(.vaeOdeFreeThetas(ui, c("lk", "lv", "lsd")),
                 c(FALSE, FALSE, TRUE))
  })

  ## Three ways a theta can reach the solve that a top-level, assignment-only
  ## scan misses.  Each misclassified the theta as ODE-free -- the UNSAFE
  ## direction, since stage 2 would then optimize it against a frozen solve.
  test_that("a state_0 initial condition is a solve dependency", {
    ui <- rxode2::assertRxUi(function() {
      ini({ lcl <- -1; linit <- 1; lsd <- 0; eta.b ~ 0.1 })
      model({
        cl <- exp(lcl)
        central_0 <- exp(linit)          # rxode2's other spelling of central(0)
        d / dt(central) <- -cl * central
        sd <- exp(lsd)
        ll(cp) ~ -log(sd) - 0.5 * ((DV - central * exp(eta.b)) / sd)^2
      })
    })
    expect_equal(.vaeOdeFreeThetas(ui, c("lcl", "linit", "lsd")),
                 c(FALSE, FALSE, TRUE))
  })

  test_that("an assignment inside if/else is a solve dependency", {
    ui <- rxode2::assertRxUi(function() {
      ini({ lcl <- -1; lsd <- 0; eta.b ~ 0.1 })
      model({
        if (TIME > 0) { cl <- exp(lcl) } else { cl <- exp(lcl) }
        d / dt(central) <- -cl * central
        sd <- exp(lsd)
        ll(cp) ~ -log(sd) - 0.5 * ((DV - central * exp(eta.b)) / sd)^2
      })
    })
    expect_equal(.vaeOdeFreeThetas(ui, c("lcl", "lsd")), c(FALSE, TRUE))
  })

  test_that("a `~` defined variable feeding d/dt is a solve dependency", {
    ## `conc` is both the endpoint variable and a d/dt input; skipping the `~`
    ## line because the name is an endpoint would drop the edge to `lv`.
    ui <- suppressWarnings(rxode2::assertRxUi(function() {
      ini({ lk <- -1; lv <- 1; lsd <- 0; eta.b ~ 0.1 })
      model({
        k <- exp(lk); v <- exp(lv)
        conc ~ central / v
        d / dt(central) <- -k * central
        d / dt(effect) <- 1 - conc * effect
        sd <- exp(lsd)
        ll(conc) ~ -log(sd) - 0.5 * ((DV - conc * exp(eta.b)) / sd)^2
      })
    }))
    expect_equal(.vaeOdeFreeThetas(ui, c("lk", "lv", "lsd")),
                 c(FALSE, FALSE, TRUE))
  })

  test_that("a linCmt() parameter is never called ODE-free", {
    ## linCmt() solves PK compartments from parameters read by NAME, which the
    ## assignment-graph scan cannot trace.  In a mixed linCmt() + ODE model the
    ## PK thetas (tcl, tv) must NOT be reported ODE-free, or stage 2 would
    ## optimize them against a frozen PK solve.
    ui <- suppressWarnings(rxode2::assertRxUi(function() {
      ini({ tcl <- 1.1; tv <- 3.9; te0 <- 2.3; tkout <- -2.3; prop.err <- 0.1
        eta.e0 ~ 0.1 })
      model({
        cl <- exp(tcl); v <- exp(tv)
        cp <- linCmt()
        e0 <- exp(te0 + eta.e0); kout <- exp(tkout); kin <- e0 * kout
        d / dt(pd) <- kin - kout * pd
        pd ~ prop(prop.err)
      })
    }))
    ## the whole model bails to "nothing ODE-free" -- conservative and safe
    expect_equal(.vaeOdeFreeThetas(ui, c("tcl", "tv", "tkout")),
                 c(FALSE, FALSE, FALSE))
  })

  test_that("stage-2 eligibility keeps every err parameter", {
    ## the historic half of the rule: an err parameter is always stage 2, and a
    ## structural theta that reaches d/dt is not
    gui <- rxode2::assertRxUi(.gMod())
    expect_equal(.vaeRegressStage2(gui, c("tv", "add.sd"), c(-1L, 0L)),
                 c(0L, 1L))
    ## an ll() model has no err rows, so only the structural proxy contributes
    lui <- rxode2::assertRxUi(.llMod())
    expect_equal(.vaeRegressStage2(lui, c("lka", "lcl", "lsd"), c(-1L, -1L, -1L)),
                 c(0L, 0L, 1L))
    ## no regressed parameters at all -> empty, not an error
    expect_equal(.vaeRegressStage2(lui, character(0), integer(0)), integer(0))
  })

  test_that("a mismatched err-index length is rejected, not recycled", {
    ## Recycling would silently mis-mask: a structural theta labelled stage 2 is
    ## then optimized against a frozen ODE, which is wrong rather than slow.
    lui <- rxode2::assertRxUi(.llMod())
    expect_error(.vaeRegressStage2(lui, c("lka", "lcl", "lsd"), c(-1L, -1L)),
                 "must match")
  })

  test_that("a model with BOTH err rows and an ll() endpoint gets both in stage 2", {
    ## The case a model-level "does this model have err parameters" short-circuit
    ## would get wrong: `add.sd` is flagged err, `lsd` is the ll() endpoint's
    ## log-density-only scale.  Both belong in stage 2; the structural thetas do
    ## not.  Under an err-only rule `lsd` would silently stay in stage 1.
    .mixMod <- function() {
      ini({ lka <- 0.4; lcl <- -0.9; lb <- log(3); add.sd <- 1.2; lsd <- log(0.5)
        eta.b ~ 0.1 })
      model({
        ka <- exp(lka); cl <- exp(lcl)
        d / dt(depot) <- -ka * depot
        d / dt(central) <- ka * depot - cl * central
        cp <- central
        cp ~ add(add.sd)
        mu <- exp(lb + eta.b) * central
        sd <- exp(lsd)
        ll(ef) ~ -log(sd) - 0.5 * ((DV - mu) / sd)^2
      })
    }
    ui <- rxode2::assertRxUi(.mixMod())
    expect_true(sum(!is.na(ui$iniDf$err) & !is.na(ui$iniDf$ntheta)) > 0L)  # has err rows
    expect_equal(.vaeRegressStage2(ui, c("lka", "lcl", "add.sd", "lsd"),
                                   c(-1L, -1L, 0L, -1L)),
                 c(0L, 0L, 1L, 1L))
  })

  test_that(".vaeDataPrep surfaces the stage-2 mask for an ll() model", {
    ui <- rxode2::assertRxUi(.llMod())
    d <- data.frame(ID = rep(1:3, each = 4), TIME = rep(c(1, 2, 4, 8), 3),
                    DV = stats::rnorm(12, 5, 1), EVID = 0L, AMT = 0)
    prep <- .vaeDataPrep(ui, d, vaeControl(nonMuTheta = "regress"))
    expect_length(prep$a, 0L)                                # no err parameter
    expect_equal(length(prep$regressStage2), length(prep$regressNames))
    ## lsd is the only stage-2 parameter; every structural one stays in stage 1
    expect_equal(prep$regressNames[prep$regressStage2 > 0L], "lsd")
  })
})
