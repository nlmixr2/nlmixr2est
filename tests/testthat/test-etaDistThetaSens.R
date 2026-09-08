## Theta sensitivities built through the declared ETAS (R/etaDistPeer.R).
##
## A declared theta reaches the model only through its own random effect, so
## d(state)/d(theta) = sum_k d(state)/d(eta_k) * d(eta_k)/d(theta).  That makes
## the sensitivity system scale with the number of declared ETAS instead of the
## number of parameters the declarations carry.
##
## Two properties are worth pinning, and they pull in opposite directions:
## the construction has to be CHEAPER (or there is no reason for it) and it has
## to compute the SAME THING (or it is not a reparameterization at all).

.edTsGammaLinCmt <- function() {
  function() {
    ini({
      lclm <- 1.5; lv1m <- 1.5; tq <- 0.9; tv2 <- 4.2
      lclrv <- -3; lv1rv <- -3
      eta.cl + eta.v1 ~ c(1, 0.3, 1)
      dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv), rate = 1/(exp(lclrv)*exp(lclm)))
      dist(eta.v1) ~ dgamma(shape = 1/exp(lv1rv), rate = 1/(exp(lv1rv)*exp(lv1m)))
      eta.q + eta.v2 ~ c(0.1, 0.01, 0.1)
      prop.sd <- 0.316
    })
    model({
      cl <- eta.cl; v <- eta.v1; q <- exp(tq + eta.q); v2 <- exp(tv2 + eta.v2)
      linCmt() ~ prop(prop.sd)
    })
  }
}

.edTsGammaOde <- function() {
  function() {
    ini({
      lka <- 0.5; lclm <- 1.5; lv1m <- 1.5
      lclrv <- -1; lv1rv <- -1
      eta.cl + eta.v ~ c(1, 0.3, 1)
      dist(eta.cl) ~ dgamma(shape = 1/exp(lclrv), rate = 1/(exp(lclrv)*exp(lclm)))
      dist(eta.v) ~ dgamma(shape = 1/exp(lv1rv), rate = 1/(exp(lv1rv)*exp(lv1m)))
      eta.ka ~ 0.1
      prop.sd <- 0.3
    })
    model({
      ka <- exp(lka + eta.ka); cl <- eta.cl; v <- eta.v
      d/dt(depot) <- -ka*depot
      d/dt(central) <- ka*depot - (cl/v)*central
      cp <- central/v
      cp ~ prop(prop.sd)
    })
  }
}

## the ui the M-step actually sees: expanded, with the declaration stash the
## expansion would otherwise destroy
.edTsUi <- function(f) {
  .ui <- rxode2::rxUiDecompress(nlmixr2est::nlmixr2(f))
  .st <- nlmixr2est:::.etaDistDeclStash(.ui, rxode2::rxUiEtaDists(.ui))
  .u2 <- rxode2::rxUiDecompress(rxode2::rxEtaDistExpand(.ui))
  nlmixr2est:::.etaDistDeclSet(.u2, .st)
  .u2
}

.edTsCols <- function(txt) {
  .ln <- strsplit(txt, "\n")[[1]]
  .ln <- .ln[grepl("_BY_THETA_[0-9]+___=", .ln)]
  stats::setNames(sub("^[^=]+=", "", .ln), sub("=.*$", "", .ln))
}

test_that("the eta-routed construction computes the same columns", {
  skip_on_cran()
  ## linCmt() has no ODE states, so BOTH routes fall through to linCmtB's own
  ## parameter sensitivities and the two must agree EXACTLY.  That is the
  ## check that this is a reparameterization of one derivative rather than a
  ## different derivative that happens to look similar.
  .a <- nlmixr2est:::rxUiGet.saemThetaSens(list(.edTsUi(.edTsGammaLinCmt())))
  .b <- nlmixr2est:::rxUiGet.etaDistThetaSens(list(.edTsUi(.edTsGammaLinCmt())))
  expect_false(is.null(.a))
  expect_false(is.null(.b))
  .ca <- .edTsCols(.a$thetaSens)
  .cb <- .edTsCols(.b$thetaSens)
  expect_identical(sort(names(.ca)), sort(names(.cb)))
  for (.n in names(.ca)) expect_identical(.cb[[.n]], .ca[[.n]], label = .n)
})

test_that("the eta-routed construction needs fewer sensitivity ODEs", {
  skip_on_cran()
  ## On a real ODE model the two diverge, which is the whole point: the theta
  ## route integrates one sensitivity per STATE per structural THETA, the eta
  ## route one per state per declared ETA.  Two states, two declared etas and
  ## five structural thetas gives 4 against 10.
  .u <- .edTsUi(.edTsGammaOde())
  .a <- nlmixr2est:::rxUiGet.saemThetaSens(list(.u))
  .b <- nlmixr2est:::rxUiGet.etaDistThetaSens(list(.edTsUi(.edTsGammaOde())))
  expect_false(is.null(.a))
  expect_false(is.null(.b))
  .nOde <- function(txt) {
    length(grep("^d/dt\\(rx__sens_", strsplit(txt, "\n")[[1]]))
  }
  ## same theta columns out of both
  expect_identical(sort(names(.edTsCols(.a$thetaSens))),
                   sort(names(.edTsCols(.b$thetaSens))))
  ## and strictly cheaper to get them
  expect_lt(.nOde(.b$thetaSens), .nOde(.a$thetaSens))
  expect_gt(.nOde(.a$thetaSens), 0L)
})

test_that("it declines, rather than guessing, when nothing is declared", {
  .one <- function() {
    ini({ tka <- 0.45; tcl <- 1.0; tv <- 3.45
          eta.ka ~ 0.6; add.sd <- 0.7 })
    model({ ka <- exp(tka + eta.ka); cl <- exp(tcl); v <- exp(tv)
            linCmt() ~ add(add.sd) })
  }
  .ui <- rxode2::rxUiDecompress(nlmixr2est::nlmixr2(.one))
  expect_null(nlmixr2est:::rxUiGet.etaDistThetaSens(list(.ui)))
})
