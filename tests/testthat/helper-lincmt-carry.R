# Shared fixtures for the linCmt() sensitivity-carry tests
# (test-focei-lincmt-carry-*.R): the canonical covariate models, a control
# setter, and the non-uniform event tables.

# The model functions are rxode2 ini()/model() DSL blocks, which lintr's
# object_usage cannot follow (every assignment looks unused, every covariate
# undefined).
# nolint start: object_usage_linter.
.carryModCov <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    eta.cl ~ 0.1
    add.sd <- 0.5
  })
  model({
    cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

.carryModTwoPair <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    eta.cl ~ 0.1
    eta.v ~ 0.1
    add.sd <- 0.5
  })
  model({
    cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
    v <- exp(tv) * (wt / 70) * exp(eta.v)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

.carryModNoCov <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    eta.cl ~ 0.1
    add.sd <- 0.5
  })
  model({
    cl <- exp(tcl) * exp(eta.cl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ add(add.sd)
  })
}

.carryModProp <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    eta.cl ~ 0.1
    prop.sd <- 0.1
  })
  model({
    cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
    v <- exp(tv)
    cp <- linCmt()
    cp ~ prop(prop.sd)
  })
}
# nolint end

.carryUiCov <- function() nlmixr2est::nlmixr2(.carryModCov)
.carryUiTwoPair <- function() nlmixr2est::nlmixr2(.carryModTwoPair)
.carryUiNoCov <- function() nlmixr2est::nlmixr2(.carryModNoCov)

.carrySetControl <- function(ui, value) {
  .ui <- rxode2::.copyUi(ui)
  .ctl <- nlmixr2est::foceiControl()
  .ctl$linCmtSensCarry <- value
  assign("control", .ctl, envir = .ui)
  .ui
}

# Non-uniform dose/observation spacing with a within-subject wt change
# (uniform spacing has hidden real pairing bugs in this effort before).
.carryEv <- function() {
  .ev <- data.frame(id = 1,
                    time = c(0, 3, 7, 15, 24, 30, 41, 50),
                    amt = c(100, 0, 100, 0, 100, 0, 100, 0),
                    evid = c(1, 0, 1, 0, 1, 0, 1, 0),
                    cmt = 1)
  .ev$wt <- ifelse(.ev$time < 20, 70, ifelse(.ev$time < 40, 85, 100))
  .ev
}

# Multi-subject version: per-subject time offsets and wt trajectories
.carryFitDat <- function(nid = 6L) {
  do.call(rbind, lapply(seq_len(nid), function(i) {
    tim <- c(0, 3, 7, 15, 24, 30, 41, 50) + (i - 1) * 0.5
    d <- data.frame(id = i, time = tim,
                    amt = c(100, 0, 100, 0, 100, 0, 100, 0),
                    evid = c(1, 0, 1, 0, 1, 0, 1, 0), cmt = 1)
    w0 <- 60 + 5 * i
    d$wt <- ifelse(d$time < 20, w0, ifelse(d$time < 40, w0 + 15, w0 + 30))
    d
  }))
}

.carryFitCtl <- function(carry, maxOut = 0L) {
  nlmixr2est::foceiControl(print = 0, maxOuterIterations = maxOut,
                           covMethod = "", calcTables = FALSE,
                           sigdig = 8, etaNudge = 0, etaNudge2 = 0,
                           rxControl = rxode2::rxControl(covsInterpolation = "nocb"),
                           linCmtSensCarry = carry)
}
