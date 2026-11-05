test_that("generated models number string literals like the data (#1299)", {
  .strCmp <- rxode2::rxModelVars(
    "v = b0 + b1*(tnom == \"1\") + b2*(tnom == \"2\") + b24*(tnom == \"24\")"
  )$strCmpParams
  .txt <- "params(b0, b2, b24, tnom)\nv = b0 + b24*(tnom == \"24\") + b2*(tnom == \"2\")"
  expect_equal(levels(rxode2::rxModelVars(.txt)$strCmpParams$tnom), c("24", "2"))
  .seeded <- .nlmixr2estStrCmpSeed(.txt, .strCmp)
  .mv <- rxode2::rxModelVars(.seeded)
  expect_equal(levels(.mv$strCmpParams$tnom), c("1", "2", "24"))
  # the params() statement still fixes the parameter order
  expect_equal(.mv$params, c("b0", "b2", "b24", "tnom"))
  expect_length(.nlmixr2estStrCmpBad(.mv, .strCmp), 0L)
  # a params() statement spanning lines stays whole
  .multi <- .nlmixr2estStrCmpSeed(sub("b2, ", "b2,\n  ", .txt, fixed = TRUE), .strCmp)
  expect_equal(levels(rxode2::rxModelVars(.multi)$strCmpParams$tnom), c("1", "2", "24"))
  # already-matching text is left alone
  expect_identical(.nlmixr2estStrCmpSeed(.seeded, .strCmp), .seeded)
})

.issue1299Data <- function() {
  set.seed(42)
  .lv <- c("0.5", "1", "2", "3", "24")
  .eff <- c("0.5" = 0, "1" = 5, "2" = 6.5, "3" = 5.5, "24" = -2)
  .d <- expand.grid(ID = 1:20, tnom = .lv, stringsAsFactors = FALSE)
  .d <- .d[order(.d$ID), ]
  .d$TIME <- seq_len(nrow(.d))
  .d$DV <- 5 + .eff[.d$tnom] + rnorm(20, 0, 2)[.d$ID] + rnorm(nrow(.d), 0, 1)
  .d$tnom <- factor(.d$tnom, levels = .lv)
  .d
}

test_that("focei string-compared covariates are not permuted (#1299)", {
  skip_on_cran()
  .d <- .issue1299Data()
  .m <- function() {
    ini({
      b0 <- 1
      b1 <- 1
      b2 <- 1
      b3 <- 1
      b24 <- 1
      addSd <- 1
      bRe ~ 1
    })
    model({
      v <- b0 + b1 * (tnom == "1") + b2 * (tnom == "2") +
        b3 * (tnom == "3") + b24 * (tnom == "24") + bRe
      v ~ add(addSd)
    })
  }
  .fit <- suppressMessages(nlmixr2(.m, .d, est = "focei", control = foceiControl(print = 0)))
  .ref <- coef(lm(DV ~ tnom, data = .d))
  expect_equal(
    unname(.fit$theta[c("b1", "b2", "b3", "b24")]),
    unname(.ref[c("tnom1", "tnom2", "tnom3", "tnom24")]),
    tolerance = 1e-2
  )
  .lv <- levels(rxode2::rxModelVars(.fit$ui)$strCmpParams$tnom)
  for (.n in c("inner", "predOnly")) {
    expect_equal(levels(rxode2::rxModelVars(.fit$foceiModel[[.n]])$strCmpParams$tnom), .lv)
  }
  .th <- .fit$theta
  .pred <- .th["b0"] + c(0, .th[c("b1", "b2", "b3", "b24")])[as.integer(.d$tnom)]
  expect_equal(.fit$PRED, unname(.pred), tolerance = 1e-4)
})
