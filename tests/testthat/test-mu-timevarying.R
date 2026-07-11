nmTest({
  # Shared time-varying / non-time-varying mu-referenced covariate handling
  # (.nlmixrSetMuRefTimeVarying + .nlmixrTimeVaryingCovariates), reused by saem,
  # the mu-referenced focei family and vae.

  test_that(".nlmixrSetMuRefTimeVarying stages the covariate split into the ui", {
    covm <- function() {
      ini({ tka<-0.45; tcl<-1; tv<-3.45; cl.wt<-0.5; eta.ka~0.6; eta.cl~0.3; eta.v~0.1; add.sd<-0.7 })
      model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl+cl.wt*WT); v<-exp(tv+eta.v); linCmt()~add(add.sd) })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(covm))
    # no time-varying covariates: muRefFinal keeps the full mu-ref set
    nlmixr2est:::.nlmixrSetMuRefTimeVarying(.ui, character(0))
    expect_true(exists("muRefFinal", envir = .ui))
    expect_equal(nrow(.ui$muRefFinal), nrow(.ui$muRefCovariateDataFrame))
    nlmixr2est:::.nlmixrRmMuRefTimeVarying(.ui)
    expect_false(exists("muRefFinal", envir = .ui, inherits = FALSE))

    # WT treated as time-varying: it drops out of the absorbed (mu-ref) set
    nlmixr2est:::.nlmixrSetMuRefTimeVarying(.ui, "WT")
    expect_false("WT" %in% .ui$muRefFinal$covariate)
    expect_equal(.ui$timeVaryingCovariates, "WT")
    nlmixr2est:::.nlmixrRmMuRefTimeVarying(.ui)
  })

  test_that("vae warns and excludes time-varying covariates from the search", {
    skip_if_not_installed("nlmixr2data")
    mod <- function() {
      ini({ tka<-0.45; tcl<-1; tv<-3.45; eta.ka~0.6; eta.cl~0.3; eta.v~0.1; add.sd<-0.7 })
      model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl); v<-exp(tv+eta.v); linCmt()~add(add.sd) })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(mod)); .ui$control <- vaeControl()
    .d <- nlmixr2data::theo_sd
    set.seed(1); .d$TVCOV <- rnorm(nrow(.d))          # varies within subject
    expect_warning(.p <- nlmixr2est:::.vaeDataPrep(.ui, .d),
                   "time-varying covariate.*excluded from automatic covariate search: TVCOV")
    expect_false("TVCOV" %in% .p$covNames)            # excluded
    expect_true("WT" %in% .p$covNames)                # subject-constant kept
    # subject-constant only: no warning
    expect_silent(suppressMessages(nlmixr2est:::.vaeDataPrep(.ui, nlmixr2data::theo_sd)))
  })

  test_that("saem recovers a non-time-varying covariate effect", {
    skip_if_not_installed("nlmixr2data")
    covm <- function() {
      ini({ tka<-0.45; tcl<-1; tv<-3.45; cl.wt<-0.5; eta.ka~0.6; eta.cl~0.3; eta.v~0.1; add.sd<-0.7 })
      model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl+cl.wt*log(WT/70)); v<-exp(tv+eta.v); linCmt()~add(add.sd) })
    }
    .f <- suppressMessages(nlmixr2(covm, nlmixr2data::theo_sd, est = "saem",
      control = saemControl(nBurn = 200, nEm = 100, nmc = 3, seed = 5, print = 0L, calcTables = FALSE)))
    # cl.wt is estimated (WT is time-invariant -> absorbed into the phi term)
    expect_true("cl.wt" %in% names(fixef(.f)))
    expect_true(is.finite(fixef(.f)[["cl.wt"]]))
  })

  test_that("a time-varying covariate is detected and kept as a model regressor", {
    skip_if_not_installed("nlmixr2data")
    # build a within-subject-varying covariate (time-dependent), then confirm
    # .nlmixrTimeVaryingCovariates flags it and saem fits its coefficient.
    tvm <- function() {
      ini({ tka<-0.45; tcl<-1; tv<-3.45; cl.tv<-0.1; eta.ka~0.6; eta.cl~0.3; eta.v~0.1; add.sd<-0.7 })
      model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl+cl.tv*TVC); v<-exp(tv+eta.v); linCmt()~add(add.sd) })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(tvm))
    .d <- nlmixr2data::theo_sd
    .d$TVC <- as.numeric(scale(.d$TIME))               # varies within subject
    .tv <- nlmixr2est:::.nlmixrTimeVaryingCovariates(.d, .ui, rxode2::rxControl())
    expect_true("TVC" %in% .tv)
    .f <- suppressMessages(nlmixr2(tvm, .d, est = "saem",
      control = saemControl(nBurn = 150, nEm = 80, nmc = 3, seed = 5, print = 0L, calcTables = FALSE)))
    expect_true("cl.tv" %in% names(fixef(.f)))
    expect_true(is.finite(fixef(.f)[["cl.tv"]]))
  })

  test_that(".saemDropMuRefFromModel(keepEtas=) absorbs covs but keeps etas", {
    covm <- function() {
      ini({ tka<-0.45; tcl<-1; tv<-3.45; cl.wt<-0.5; eta.ka~0.6; eta.cl~0.3; eta.v~0.1; add.sd<-0.7 })
      model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl+cl.wt*WT); v<-exp(tv+eta.v); linCmt()~add(add.sd) })
    }
    .ui <- rxode2::rxUiDecompress(rxode2::rxode2(covm))
    nlmixr2est:::.nlmixrSetMuRefTimeVarying(.ui, character(0))
    on.exit(nlmixr2est:::.nlmixrRmMuRefTimeVarying(.ui), add = TRUE)
    # saem: mu-ref etas and covariates both dropped (phi model)
    .saem <- vapply(nlmixr2est:::.saemDropMuRefFromModel(.ui, keepEtas = FALSE),
                    function(e) paste(deparse(e), collapse = ""), character(1))
    expect_true(any(grepl("cl <- exp\\(tcl\\)$", .saem)))
    # fsaem inner: covariate absorbed, but eta kept (distinct from the plain
    # focei inner where the covariate stays in the model)
    .inner <- vapply(nlmixr2est:::.saemDropMuRefFromModel(.ui, keepEtas = TRUE),
                     function(e) paste(deparse(e), collapse = ""), character(1))
    expect_true(any(grepl("cl <- exp\\(tcl \\+ eta.cl\\)$", .inner)))
    expect_false(any(grepl("cl.wt", .inner)))          # covariate absorbed
  })

  test_that("mufocei recovers non-time-varying and time-varying covariate effects", {
    skip_if_not_installed("nlmixr2data")
    .fc <- foceiControl(print = 0L, calcTables = FALSE,
                        maxInnerIterations = 30L, maxOuterIterations = 40L)
    # non-time-varying covariate (absorbed into the phi term via the mu2 hook)
    covm <- function() {
      ini({ tka<-0.45; tcl<-1; tv<-3.45; cl.wt<-0.5; eta.ka~0.6; eta.cl~0.3; eta.v~0.1; add.sd<-0.7 })
      model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl+cl.wt*log(WT/70)); v<-exp(tv+eta.v); linCmt()~add(add.sd) })
    }
    .m1 <- suppressMessages(nlmixr2(covm, nlmixr2data::theo_sd, est = "mufocei", control = .fc))
    expect_true("cl.wt" %in% names(fixef(.m1)))
    expect_true(is.finite(fixef(.m1)[["cl.wt"]]))

    # time-varying covariate (kept as a beta regressor in the model)
    tvm <- function() {
      ini({ tka<-0.45; tcl<-1; tv<-3.45; cl.tv<-0.1; eta.ka~0.6; eta.cl~0.3; eta.v~0.1; add.sd<-0.7 })
      model({ ka<-exp(tka+eta.ka); cl<-exp(tcl+eta.cl+cl.tv*TVC); v<-exp(tv+eta.v); linCmt()~add(add.sd) })
    }
    .d <- nlmixr2data::theo_sd; .d$TVC <- as.numeric(scale(.d$TIME))
    .m2 <- suppressMessages(nlmixr2(tvm, .d, est = "mufocei", control = .fc))
    expect_true("cl.tv" %in% names(fixef(.m2)))
    expect_true(is.finite(fixef(.m2)[["cl.tv"]]))
  })
})
