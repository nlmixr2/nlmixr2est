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

# Jump channels: eta.cl on a covariate-driven slot, eta.f on a bioavailability
# whose dlnF depends on wt (logit-additive), eta.lag on a modeled lag over a
# time-varying kernel
.carryModJump <- function() {
  ini({
    tcl <- log(2)
    tv <- log(20)
    tka <- log(1.2)
    tf <- -0.5
    tlag <- log(0.5)
    eta.cl ~ 0.1
    eta.f ~ 0.1
    eta.lag ~ 0.1
    add.sd <- 0.5
  })
  model({
    cl <- exp(tcl) * (wt / 70)^0.75 * exp(eta.cl)
    v <- exp(tv)
    ka <- exp(tka)
    f(depot) <- expit(tf + eta.f + (wt - 70) / 70)
    alag(depot) <- exp(tlag + eta.lag)
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
.carryUiJump <- function() suppressMessages(nlmixr2est::nlmixr2(.carryModJump))

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

# Jump-channel fixtures shared by test-focei-lincmt-carry-jump.R and
# test-focei-lincmt-carry-trans.R
.carryJumpPars <- c(`THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = log(1.2),
                    `THETA[4]` = -0.5, `THETA[5]` = log(0.5), `THETA[6]` = 0.5,
                    `ETA[1]` = 0.3, `ETA[2]` = -0.2, `ETA[3]` = 0.1)

# max relative error of every eta's substituted gradient vs central FD
.carryJumpFd <- function(mod, pars, ev, carry = "auto", interp = "nocb") { # nolint: object_usage_linter.
  ui <- suppressMessages(nlmixr2est::nlmixr2(mod))
  u <- rxode2::.copyUi(ui)
  assign("control", nlmixr2est::foceiControl(linCmtSensCarry = carry), envir = u)
  txt <- suppressMessages(u$foceiEnv)$..inner
  m <- suppressWarnings(rxode2::rxode2(txt))
  slv <- function(q) {
    rxode2::rxSolve(m, params = q, events = ev, returnType = "data.frame",
                    covsInterpolation = interp)
  }
  r0 <- slv(pars)
  h <- 1e-5
  neta <- sum(grepl("^ETA", names(pars)))
  err <- vapply(seq_len(neta), function(k) {
    en <- paste0("ETA[", k, "]")
    a <- pars; a[en] <- a[en] + h
    b <- pars; b[en] <- b[en] - h
    fd <- (slv(a)$rx_pred_ - slv(b)$rx_pred_) / (2 * h)
    got <- r0[[paste0("rx__sens_rx_pred__BY_ETA_", k, "___")]]
    max(abs(got - fd) / (abs(fd) + 1e-8))
  }, numeric(1))
  list(err = err, txt = txt)
}
