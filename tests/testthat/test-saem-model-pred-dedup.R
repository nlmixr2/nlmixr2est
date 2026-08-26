nmTest({
  # rxUiGet.saemModelPred() emits the mu-reference replacement block before the
  # model body and the THETA/ETA alias block after it.  The alias lines that
  # exactly duplicate a replacement line are dropped; the ones whose rhs differs
  # (mu-referenced `THETA[k] + ETA[j]` vs the split `THETA[k]`) are kept.

  .mod <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      pow <- 1.1
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv)
      d/dt(depot) <- -ka * depot
      d/dt(central) <- ka * depot - cl / v * central
      cp <- central / v * pow
      cp ~ add(add.sd)
    })
  }

  .ui <- .mod()
  .pred <- suppressMessages(rxUiGet.saemModelPred(list(.ui))$predOnly)
  .lines <- strsplit(rxode2::rxNorm(.pred), "\n")[[1]]
  # count the assignments to a single variable
  .nAssign <- function(v) {
    sum(grepl(paste0("^", gsub(".", "\\.", v, fixed=TRUE), "="), .lines))
  }

  test_that(".normAssign() normalizes only parseable single assignments", {
    expect_equal(.normAssign(c("tka<-THETA[1]", "tka <-  THETA[1]")),
                 c("tka <- THETA[1]", "tka <- THETA[1]"))
    # an empty rhs (a placeholder non-mu eta) does not parse and is left alone
    expect_equal(.normAssign("eta.cl1 <- "), "eta.cl1 <- ")
    expect_equal(.normAssign(character(0)), character(0))
  })

  test_that("saem predOnly emits no duplicated line", {
    expect_false(any(duplicated(.lines)))
  })

  test_that("a non-mu-referenced parameter is aliased once", {
    # tv/pow/add.sd carry no eta, so the replacement block and the alias block
    # emit the identical assignment and only one survives
    expect_equal(.nAssign("tv"), 1L)
    expect_equal(.nAssign("pow"), 1L)
    expect_equal(.nAssign("add.sd"), 1L)
  })

  test_that("a mu-referenced parameter keeps both the combined and split form", {
    expect_equal(.nAssign("tka"), 2L)
    expect_true(any(grepl("^tka=THETA\\[1\\]\\+ETA\\[1\\];?$", .lines)))
    expect_true(any(grepl("^tka=THETA\\[1\\];?$", .lines)))
    # the combined form comes first, so the trailing split alias still wins the
    # reported `tka` column exactly as it did before the de-duplication
    expect_lt(which(grepl("^tka=THETA\\[1\\]\\+ETA\\[1\\];?$", .lines)),
              which(grepl("^tka=THETA\\[1\\];?$", .lines)))
    expect_equal(.nAssign("eta.ka"), 1L)
    expect_equal(.nAssign("eta.cl"), 1L)
  })

  test_that("the dropped aliases were dead stores", {
    # rebuild the pre-de-duplication model: the full THETA/ETA alias block
    # re-appended after the model body, ahead of the trailing cmt()/dvid()
    .thetaEta <- vapply(.uiGetThetaEta(.ui), deparse1, character(1),
                        USE.NAMES=FALSE)
    .split <- max(which(!grepl("^(cmt|dvid)\\(", .lines)))
    .dup <- suppressMessages(rxode2::rxode2(paste(
      c(.lines[seq_len(.split)], .thetaEta, .lines[-seq_len(.split)]),
      collapse="\n")))
    # the solve column layout must not shift
    expect_equal(.pred$params, .dup$params)
    expect_equal(.pred$lhs, .dup$lhs)
    expect_equal(.pred$state, .dup$state)
    .p <- c("THETA[1]"=0.45, "THETA[2]"=1, "THETA[3]"=3.45, "THETA[4]"=1.1,
            "THETA[5]"=0.7, "ETA[1]"=0.2, "ETA[2]"=-0.3)
    .et <- rxode2::et(amt=320, cmt="depot")
    .et <- rxode2::et(.et, seq(0, 24, by=0.5))
    .s1 <- rxode2::rxSolve(.pred, .p, .et, returnType="data.frame",
                           addDosing=FALSE)
    .s2 <- rxode2::rxSolve(.dup, .p, .et, returnType="data.frame",
                           addDosing=FALSE)
    expect_equal(.s1, .s2)
  })
})
