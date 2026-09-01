# linCmt() sensitivity carry on the other linCmt() parameterizations
# (trans 2/3/4/5/10/11): eta+covariate on a kernel slot or on a slot that
# enters the observation scaling (Vc), gradient vs finite differences
# through a real rxSolve() of the generated inner model.  The f()/alag()
# channels are in test-focei-lincmt-carry-jump.R (); the FD helper is in
# helper-lincmt-carry.R; shared fixtures in helper-lincmt-carry.R.

test_that("the other parameterizations carry exactly (eta+covariate on each slot kind)", {
  skip_if_not(.rxFoceiLinCmtCarryCapable())
  ev <- .carryEv()
  cases <- list(
    list(name = "1-cmt k/v", mod = function() { # nolint: indentation_linter.
      ini({
        tk <- log(0.1)
        tv <- log(20)
        eta.k ~ 0.1
        add.sd <- 0.5
      })
      model({
        k <- exp(tk) * (wt / 70)^0.75 * exp(eta.k)
        v <- exp(tv)
        cp <- linCmt()
        cp ~ add(add.sd)
      })
    }, pars = c(`THETA[1]` = log(0.1), `THETA[2]` = log(20), `THETA[3]` = 0.5, `ETA[1]` = 0.3)),
    list(name = "2-cmt k/k12/k21", mod = function() {
      ini({
        tk <- log(0.1)
        tv <- log(20)
        tk12 <- log(0.05)
        tk21 <- log(0.03)
        eta.k12 ~ 0.1
        add.sd <- 0.5
      })
      model({
        k <- exp(tk)
        v <- exp(tv)
        k12 <- exp(tk12) * (wt / 70) * exp(eta.k12)
        k21 <- exp(tk21)
        cp <- linCmt()
        cp ~ add(add.sd)
      })
    }, pars = c(
      `THETA[1]` = log(0.1), `THETA[2]` = log(20), `THETA[3]` = log(0.05),
      `THETA[4]` = log(0.03), `THETA[5]` = 0.5, `ETA[1]` = 0.3
    )),
    list(name = "2-cmt cl/v/q/vss", mod = function() { # nolint: indentation_linter.
      ini({
        tcl <- log(2)
        tv <- log(20)
        tq <- log(1)
        tvss <- log(60)
        eta.vss ~ 0.1
        add.sd <- 0.5
      })
      model({
        cl <- exp(tcl)
        v <- exp(tv)
        q <- exp(tq)
        vss <- exp(tvss) * (wt / 70) * exp(eta.vss)
        cp <- linCmt()
        cp ~ add(add.sd)
      })
    }, pars = c(
      `THETA[1]` = log(2), `THETA[2]` = log(20), `THETA[3]` = log(1),
      `THETA[4]` = log(60), `THETA[5]` = 0.5, `ETA[1]` = 0.3
    )),
    list(name = "2-cmt alpha/beta/k21", mod = function() { # nolint: indentation_linter.
      ini({
        ta <- log(0.5)
        tb <- log(0.02)
        tk21 <- log(0.05)
        tv <- log(20)
        eta.b ~ 0.1
        add.sd <- 0.5
      })
      model({
        alpha <- exp(ta)
        beta <- exp(tb) * (wt / 70)^0.5 * exp(eta.b)
        k21 <- exp(tk21)
        v <- exp(tv)
        cp <- linCmt()
        cp ~ add(add.sd)
      })
    }, pars = c(
      `THETA[1]` = log(0.5), `THETA[2]` = log(0.02), `THETA[3]` = log(0.05),
      `THETA[4]` = log(20), `THETA[5]` = 0.5, `ETA[1]` = 0.3
    )),
    list(name = "2-cmt alpha/beta/aob", mod = function() { # nolint: indentation_linter.
      ini({
        ta <- log(0.5)
        tb <- log(0.02)
        taob <- log(2)
        tv <- log(20)
        eta.aob ~ 0.1
        add.sd <- 0.5
      })
      model({
        alpha <- exp(ta)
        beta <- exp(tb)
        aob <- exp(taob) * (wt / 70) * exp(eta.aob)
        v <- exp(tv)
        cp <- linCmt()
        cp ~ add(add.sd)
      })
    }, pars = c(
      `THETA[1]` = log(0.5), `THETA[2]` = log(0.02), `THETA[3]` = log(2),
      `THETA[4]` = log(20), `THETA[5]` = 0.5, `ETA[1]` = 0.3
    )),
    list( # nolint: indentation_linter.
      name = "2-cmt A/B/alpha/beta (B in Vc)", mod = function() { # nolint: indentation_linter.
        ini({ # nolint: indentation_linter.
          ta <- log(0.5)
          tb <- log(0.02)
          tA <- log(0.04)
          tB <- log(0.01)
          eta.b ~ 0.1
          add.sd <- 0.5
        })
        model({
          alpha <- exp(ta)
          beta <- exp(tb)
          a <- exp(tA)
          b <- exp(tB) * (wt / 70)^-1 * exp(eta.b)
          cp <- linCmt()
          cp ~ add(add.sd)
        })
      }, pars = c(
        `THETA[1]` = log(0.5), `THETA[2]` = log(0.02), `THETA[3]` = log(0.04),
        `THETA[4]` = log(0.01), `THETA[5]` = 0.5, `ETA[1]` = 0.3
      )
    ),
    list( # nolint: indentation_linter.
      name = "2-cmt v/B/alpha/beta (v in Vc)", mod = function() { # nolint: indentation_linter.
        ini({ # nolint: indentation_linter.
          ta <- log(0.5)
          tb <- log(0.02)
          tv <- log(25)
          tB <- log(0.01)
          eta.v ~ 0.1
          add.sd <- 0.5
        })
        model({
          alpha <- exp(ta)
          beta <- exp(tb)
          v <- exp(tv) * (wt / 70) * exp(eta.v)
          b <- exp(tB)
          cp <- linCmt()
          cp ~ add(add.sd)
        })
      }, pars = c(
        `THETA[1]` = log(0.5), `THETA[2]` = log(0.02), `THETA[3]` = log(25),
        `THETA[4]` = log(0.01), `THETA[5]` = 0.5, `ETA[1]` = 0.3
      )
    ),
    list( # nolint: indentation_linter.
      name = "3-cmt v/B/C/alpha/beta/gamma (C in Vc)", mod = function() { # nolint: indentation_linter.
        ini({ # nolint: indentation_linter.
          ta <- log(1)
          tb <- log(0.1)
          tg <- log(0.01)
          tv <- log(25)
          tB <- log(0.01)
          tC <- log(0.005)
          eta.c ~ 0.1
          add.sd <- 0.5
        })
        model({
          alpha <- exp(ta)
          beta <- exp(tb)
          gamma <- exp(tg)
          v <- exp(tv)
          b <- exp(tB)
          c <- exp(tC) * (wt / 70)^-1 * exp(eta.c)
          cp <- linCmt()
          cp ~ add(add.sd)
        })
      }, pars = c(
        `THETA[1]` = log(1), `THETA[2]` = log(0.1), `THETA[3]` = log(0.01),
        `THETA[4]` = log(25), `THETA[5]` = log(0.01), `THETA[6]` = log(0.005),
        `THETA[7]` = 0.5, `ETA[1]` = 0.3
      )
    )
  )
  for (cs in cases) {
    rC <- .carryJumpFd(cs$mod, cs$pars, ev)
    expect_true(grepl("rx_lcCarryAdv_", rC$txt, fixed = TRUE), label = cs$name)
    expect_lt(max(rC$err), 1e-6, label = cs$name)
    rN <- .carryJumpFd(cs$mod, cs$pars, ev, carry = "none")
    expect_gt(max(rN$err), 1e-3, label = cs$name)
  }
})
