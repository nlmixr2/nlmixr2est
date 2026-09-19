nmTest({
  # Simulate from a known truth whose likelihood parameters carry an eta or a
  # covariate, then check saem recovers it and agrees with FOCEi.
  .simModeledResid <- function(model, nSub = 60L, seed = 20260912L) {
    set.seed(seed)
    .ev <- rxode2::et(amt = 320, cmt = "depot") |>
      rxode2::et(c(0.25, 0.5, 1, 2, 3.5, 5, 7, 9, 12, 24)) |>
      rxode2::et(id = seq_len(nSub))
    .cov <- data.frame(id = seq_len(nSub), WT = round(runif(nSub, 50, 110)))
    .ev <- merge(as.data.frame(.ev), .cov, by = "id")
    .s <- rxode2::rxSolve(model, .ev, addDosing = TRUE, returnType = "data.frame")
    # WT is only carried in the solve when the model uses it
    .s <- merge(.s[, setdiff(names(.s), "WT")], .cov, by = "id")
    .d <- data.frame(
      ID = .s$id,
      TIME = .s$time,
      EVID = .s$evid,
      AMT = ifelse(is.na(.s$amt), 0, .s$amt),
      CMT = ifelse(.s$evid == 1, 1L, 2L),
      DV = ifelse(.s$evid == 1, 0, .s$sim),
      WT = .s$WT
    )
    .d[order(.d$ID, .d$TIME, -.d$EVID), ]
  }
  .noTemporaryEta <- function(fit) {
    expect_false(any(grepl(
      "^rx[.](eta|l)[.]|^rxBoundedTr",
      c(names(fit), names(fit$eta), rownames(fit$omega), names(fixef(fit)))
    )))
    expect_true(any(grepl("temporary eta for eta-less likelihood theta", fit$runInfo, fixed = TRUE)))
  }
  .fitBoth <- function(mFit, d, nu = c(2, 2, 2)) {
    list(
      saem = .nlmixr(mFit, d, est = "saem", control = saemControl(seed = 42L, print = 0L, covMethod = "", nu = nu)),
      focei = .nlmixr(mFit, d, est = "focei", control = foceiControl(print = 0L, covMethod = ""))
    )
  }
  .base <- function(extraIni, lines) {
    eval(parse(
      text = sprintf(
        "function() {
      ini({ tka <- log(1.5); tcl <- log(2.7); tv <- log(30); %s; eta.ka ~ 0.2; eta.cl ~ 0.1 })
      model({ ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv); cp <- linCmt(); %s })
    }",
        extraIni,
        lines
      )
    ))
  }
  .starts <- list(tka = log(1), tcl = log(2), tv = log(25))

  test_that("saem recovers an eta on the additive residual SD", {
    mTrue <- .base("add.sd <- 0.5; eta.sd ~ 0.2", "a <- add.sd * exp(eta.sd); cp ~ add(a)")
    .d <- .simModeledResid(mTrue)
    # an eta on the residual SD mixes slowly: at the default nu=2 most seeds stop
    # short of the optimum, at nu=4 every seed tried reaches it
    .f <- .fitBoth(do.call(rxode2::ini, c(list(mTrue), .starts, list(add.sd = 1, eta.sd = 0.1))), .d, nu = c(4, 4, 4))
    .noTemporaryEta(.f$saem)
    expect_equal(
      unname(fixef(.f$saem)[c("tka", "tcl", "tv")]),
      unname(fixef(.f$focei)[c("tka", "tcl", "tv")]),
      tolerance = 0.05
    )
    expect_equal(unname(fixef(.f$saem)[["add.sd"]]), unname(fixef(.f$focei)[["add.sd"]]), tolerance = 0.15)
  })

  test_that("saem recovers a plain + dnorm() residual SD", {
    mTrue <- .base("add.sd <- 0.5", "cp ~ add(add.sd) + dnorm()")
    .d <- .simModeledResid(mTrue)
    .f <- .fitBoth(do.call(rxode2::ini, c(list(mTrue), .starts, list(add.sd = 1))), .d)
    .noTemporaryEta(.f$saem)
    expect_equal(unname(fixef(.f$saem)[["add.sd"]]), unname(fixef(.f$focei)[["add.sd"]]), tolerance = 0.15)
  })

  test_that("saem recovers add + prop under + dnorm()", {
    mTrue <- .base("add.sd <- 0.2; prop.sd <- 0.15", "cp ~ add(add.sd) + prop(prop.sd) + dnorm()")
    .d <- .simModeledResid(mTrue)
    .f <- .fitBoth(do.call(rxode2::ini, c(list(mTrue), .starts, list(add.sd = 0.5, prop.sd = 0.3))), .d)
    .noTemporaryEta(.f$saem)
    expect_equal(
      unname(fixef(.f$saem)[c("add.sd", "prop.sd")]),
      unname(fixef(.f$focei)[c("add.sd", "prop.sd")]),
      tolerance = 0.2
    )
  })

  test_that("saem recovers a covariate on the additive residual SD", {
    mTrue <- .base("add.sd <- 0.5; cov.sd <- 0.2", "a <- add.sd + (WT - 80) / 20 * cov.sd; cp ~ add(a)")
    .d <- .simModeledResid(mTrue)
    .f <- .fitBoth(do.call(rxode2::ini, c(list(mTrue), .starts, list(add.sd = 1, cov.sd = 0.05))), .d)
    .noTemporaryEta(.f$saem)
    expect_equal(
      unname(fixef(.f$saem)[c("add.sd", "cov.sd")]),
      unname(fixef(.f$focei)[c("add.sd", "cov.sd")]),
      tolerance = 0.2
    )
  })

  test_that("saem recovers a boxCox lambda under + dnorm()", {
    mTrue <- .base("add.sd <- 0.3; lam <- 0.5", "cp ~ add(add.sd) + boxCox(lam) + dnorm()")
    .d <- .simModeledResid(mTrue)
    .d <- .d[.d$EVID == 1 | .d$DV > 0, ]
    mFit <- do.call(rxode2::ini, c(list(mTrue), .starts, list(add.sd = 0.5, lam = 1)))
    .f <- .nlmixr(mFit, .d, est = "saem", control = saemControl(seed = 42L, print = 0L, covMethod = ""))
    .noTemporaryEta(.f)
    expect_equal(unname(fixef(.f)[c("add.sd", "lam")]), c(0.3, 0.5), tolerance = 0.3)
    # the temporary etas' own mean update: their variance shrinks, and the reported
    # theta is the kernel's (refining their mu split the two and gave lam ~ 1)
    .ph <- .f$parHistData[.f$parHistData$type == "Unscaled", ]
    .last <- .ph[nrow(.ph), ]
    expect_lt(.last[["V(rx.eta.lam)"]], 0.1)
    expect_equal(unname(fixef(.f)[["lam"]]), .last[["lam"]], tolerance = 0.05)
  })

  test_that("saem recovers an ll() residual SD", {
    mTrue <- .base(
      "lsd <- log(0.5)",
      "sd <- exp(lsd); ll(err) ~ -log(sd) - 0.5 * log(2 * pi) - 0.5 * ((DV - cp) / sd)^2"
    )
    mSim <- .base("add.sd <- 0.5", "cp ~ add(add.sd)")
    .d <- .simModeledResid(mSim)
    .f <- .fitBoth(do.call(rxode2::ini, c(list(mTrue), .starts, list(lsd = log(1)))), .d)
    .noTemporaryEta(.f$saem)
    expect_equal(unname(fixef(.f$saem)[["lsd"]]), unname(fixef(.f$focei)[["lsd"]]), tolerance = 0.15)
  })

  test_that("one modeled endpoint promotes a two-endpoint fit", {
    mTrue <- function() {
      ini({
        tka <- log(1.5); tcl <- log(2.7); tv <- log(30); add.sd <- 0.5; add.sd2 <- 0.3
        eta.ka ~ 0.2; eta.cl ~ 0.1; eta.sd ~ 0.25
      })
      model({
        ka <- exp(tka + eta.ka); cl <- exp(tcl + eta.cl); v <- exp(tv)
        cp <- linCmt()
        cp2 <- 2 * cp
        a <- add.sd * exp(eta.sd)
        cp ~ add(a)
        cp2 ~ add(add.sd2)
      })
    }
    ui <- rxode2::rxode2(mTrue)
    expect_equal(.saemModeledResidualCond(ui), "cp")
    ui <- suppressWarnings(.preProcessSaemModeledResid(ui, "saem", NULL, NULL)$ui)
    expect_equal(as.character(ui$predDf$distribution), c("dnorm", "norm"))
    expect_equal(ui$saemResMod, c(0L, 0L))
    .txt <- paste(deparse(ui$saemModel0), collapse = "\n")
    expect_match(.txt, "llikNorm", fixed = TRUE)
    expect_match(.txt, "llikXNorm", fixed = TRUE)
  })
})
