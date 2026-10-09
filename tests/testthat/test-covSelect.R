# every warning an expression gives, in order, and its value
.covSelectWarnings <- function(expr) {
  .acc <- new.env(parent = emptyenv())
  .acc$w <- character(0)
  .v <- withCallingHandlers(
    expr,
    warning = function(w) {
      .acc$w <- c(.acc$w, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = .v, warnings = .acc$w)
}

.covSelectEnv <- function(covR = diag(c(0.04, 0.09)), covS = diag(c(0.05, 0.1)), covRS = diag(c(0.045, 0.095))) {
  .env <- new.env(parent = emptyenv())
  .env$covR <- covR
  .env$covS <- covS
  .env$.covSinv <- covS
  .env$covRS <- covRS
  .env
}

test_that(".covSelectFocei() installs the requested covariance when its matrices are usable", {
  # "r,s": the sandwich
  .e <- .covSelectEnv()
  expect_identical(.covSelectFocei(.e, 1L, 1L, 1L, "r", "s", FALSE, FALSE, 1e-5), list(slot = 1L, label = "r,s"))
  expect_identical(.e$cov, .e$covRS)
  expect_identical(.e$covMethod, "r,s")
  # "r"
  .e <- .covSelectEnv()
  expect_identical(.covSelectFocei(.e, 2L, 1L, 0L, "r", "s", FALSE, FALSE, 1e-5), list(slot = 2L, label = "r"))
  expect_identical(.e$cov, .e$covR)
  # "s": the S inverse; a request for S alone still says it uses S
  .e <- .covSelectEnv()
  expect_warning(
    .r <- .covSelectFocei(.e, 3L, 0L, 1L, "r", "s", FALSE, FALSE, 1e-5),
    "using S matrix to calculate covariance"
  )
  expect_identical(.r, list(slot = 3L, label = "s"))
  expect_identical(.e$cov, .e$.covSinv)
})

test_that(".covSelectFocei() falls back between R and S and labels what it used", {
  # R not usable: S
  .e <- .covSelectEnv()
  expect_warning(.r <- .covSelectFocei(.e, 1L, 2L, 1L, "r", "s", FALSE, FALSE, 1e-5), "using S matrix")
  expect_identical(.r$label, "s")
  expect_identical(.e$cov, .e$.covSinv)
  # "r" whose R is not usable: S
  .e <- .covSelectEnv()
  expect_warning(.r <- .covSelectFocei(.e, 2L, 3L, 1L, "r", "s", FALSE, FALSE, 1e-5), "using S matrix")
  expect_identical(.r$slot, 3L)
  # S not positive definite: R
  .e <- .covSelectEnv()
  expect_warning(.r <- .covSelectFocei(.e, 1L, 1L, 2L, "r", "s", FALSE, FALSE, 1e-5), "using R matrix")
  expect_identical(.r, list(slot = 2L, label = "r"))
  expect_identical(.e$cov, .e$covR)
  # S not computed: R, with the console note
  .e <- .covSelectEnv()
  expect_output(
    expect_warning(.r <- .covSelectFocei(.e, 1L, 1L, 3L, "r", "s", FALSE, FALSE, 1e-5), "using R matrix"),
    "S matrix calculation failed; Switch to R-matrix covariance"
  )
  expect_identical(.r$slot, 2L)
  # neither usable: none
  .e <- .covSelectEnv()
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 2L, 2L, "r", "s", FALSE, FALSE, 1e-5))
  expect_identical(.r$value, list(slot = 0L, label = "failed"))
  expect_identical(.r$warnings, c("cannot calculate covariance", "covariance step failed"))
  expect_false(exists("cov", envir = .e, inherits = FALSE))
})

test_that(".covSelectFocei() checks a doubtful sandwich against R and S", {
  # a repaired R: S, the one not repaired
  .e <- .covSelectEnv()
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 1L, 1L, "r+", "s", TRUE, FALSE, 1e-5))
  expect_identical(.r$value, list(slot = 3L, label = "s"))
  expect_identical(.e$cov, .e$covS)
  expect_identical(.r$warnings, "using S matrix to calculate covariance, can check sandwich or R matrix with $covRS and $covR")
  # a repaired S: R, with the R repair warning only when R is in the covariance
  .e <- .covSelectEnv()
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 1L, 1L, "r", "|s|", TRUE, FALSE, 1e-5))
  expect_identical(.r$value$label, "r")
  expect_identical(.e$cov, .e$covR)
  # both repaired: the diagonal sums decide (R x 2, S x 4 against the sandwich)
  .e <- .covSelectEnv(covR = diag(c(0.04, 0.09)), covS = diag(c(0.05, 0.1)), covRS = diag(c(0.5, 0.5)))
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 1L, 1L, "r+", "s+", TRUE, FALSE, 1e-5))
  expect_identical(.r$value$label, "r+")
  expect_identical(.e$cov, .e$covR)
  expect_identical(
    .r$warnings,
    c(
      "R matrix non-positive definite but corrected (because of cholAccept)",
      "using R matrix to calculate covariance, can check sandwich or S matrix with $covRS and $covS"
    )
  )
  .e <- .covSelectEnv(covR = diag(c(0.04, 0.09)), covS = diag(c(0.05, 0.1)), covRS = diag(c(0.01, 0.01)))
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 1L, 1L, "|r|", "s+", TRUE, FALSE, 1e-5))
  expect_identical(.r$value, list(slot = 1L, label = "|r|,s+"))
  expect_identical(.e$cov, .e$covRS)
  expect_identical(
    .r$warnings,
    c(
      "R matrix non-positive definite but corrected by R = sqrtm(R%*%R)",
      "S matrix non-positive definite but corrected (because of cholAccept)",
      "since sandwich matrix is corrected, you may compare to $covR or $covS if you wish"
    )
  )
  # an unrepaired sandwich with a variance below covSmall is checked the same way: kept
  # here, its diagonal sum being the smallest of the three
  .e <- .covSelectEnv(covRS = diag(c(1e-6, 0.095)))
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 1L, 1L, "r", "s", FALSE, FALSE, 1e-5))
  expect_identical(.r$value$label, "r,s")
  expect_identical(.e$cov, .e$covRS)
  # and replaced by R when the sandwich's sum is the largest and R's is below S's
  .e <- .covSelectEnv(covRS = diag(c(1e-6, 2)))
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 1L, 1L, "r", "s", FALSE, FALSE, 1e-5))
  expect_identical(.r$value$label, "r")
  expect_identical(.e$cov, .e$covR)
  expect_identical(.covSandwichChoice(diag(2), diag(2), diag(2), "r", "s", FALSE, 1e-5), "covRS")
})

test_that(".covSelectFocei() refuses a covariance of all tiny variances and reports S problems", {
  .e <- .covSelectEnv(covR = diag(c(1e-8, 1e-9)))
  .r <- .covSelectWarnings(.covSelectFocei(.e, 2L, 1L, 0L, "r", "s", FALSE, FALSE, 1e-5))
  expect_identical(.r$value$slot, 0L)
  expect_identical(.r$warnings, c("The variance of all elements are unreasonably small, <1e-7", "covariance step failed"))
  .e <- .covSelectEnv()
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 1L, 1L, "r", "s", FALSE, TRUE, 1e-5))
  expect_identical(.r$value$label, "r,s")
  expect_identical(.r$warnings, "S matrix had problems solving for some subject and parameters")
})

test_that(".covFallbackCheck() takes a named list of methods and their fallbacks", {
  expect_identical(.covFallbackCheck(NULL), list())
  expect_identical(.covFallbackCheck(list()), list())
  expect_identical(.covFallbackCheck(list("r,s" = "s", r = NULL)), list("r,s" = "s", r = character(0)))
  expect_error(.covFallbackCheck(c("r", "s")), "must be a named list")
  expect_error(.covFallbackCheck(list("r", "s")), "must be a named list")
  expect_error(.covFallbackCheck(list(r = "s", r = "r,s")), "uniquely named")
  expect_error(.covFallbackCheck(list(sa = "r")), "names a method without fallbacks: \"sa\"")
  expect_error(.covFallbackCheck(list("r,s" = c("s", "s"))), "must be distinct method names")
  expect_error(.covFallbackCheck(list("r,s" = NA_character_)), "must be distinct method names")
  expect_error(.covFallbackCheck(list("r,s" = "r,s")), "cannot fall back to \"r,s\"")
  expect_error(.covFallbackCheck(list(r = "analytic")), "cannot fall back to \"analytic\"")
  # the default reproduces the established fallbacks
  expect_identical(
    foceiControl()$covFallback,
    list("r,s" = c("r", "s"), r = "s", s = character(0), analytic = c("r,s", "r", "s"))
  )
  expect_identical(foceiControl(covFallback = list(r = "s"))$covFallback, list(r = "s"))
})

test_that(".covSelectFocei() falls back only to the listed methods and records what it tried", {
  # R not usable and "s" not listed: none
  .e <- .covSelectEnv()
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 2L, 0L, "r", "s", FALSE, FALSE, 1e-5, fallback = "r"))
  expect_identical(.r$value$slot, 0L)
  expect_identical(.e$covTried, data.frame(method = "r,s", outcome = "R not positive definite"))
  # S not usable and "r" not listed: none
  .e <- .covSelectEnv()
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 1L, 2L, "r", "s", FALSE, FALSE, 1e-5, fallback = "s"))
  expect_identical(.r$value$slot, 0L)
  expect_identical(.r$warnings, c("cannot calculate covariance", "covariance step failed"))
  # the sandwich check may not pick an unlisted method: the sandwich stays
  .e <- .covSelectEnv()
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 1L, 1L, "r+", "s", TRUE, FALSE, 1e-5, fallback = "r"))
  expect_identical(.r$value$label, "r+,s")
  expect_identical(.e$cov, .e$covRS)
  # a fallback that is used, and the request it replaced
  .e <- .covSelectEnv()
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 2L, 1L, "r", "s", FALSE, FALSE, 1e-5))
  expect_identical(
    .e$covTried,
    data.frame(method = c("r,s", "s"), outcome = c("R not positive definite", "used"))
  )
  .e <- .covSelectEnv()
  .r <- .covSelectWarnings(.covSelectFocei(.e, 1L, 1L, 1L, "r+", "s", TRUE, FALSE, 1e-5))
  expect_identical(
    .e$covTried,
    data.frame(method = c("r,s", "s"), outcome = c("sandwich not used (covSmall check)", "used"))
  )
})

test_that("a declined analytic covariance falls back as covFallback$analytic lists", {
  skip_on_cran()
  # linCmt() is outside the analytic covariance's scope
  .lin <- function() {
    ini({
      tka <- 0.45
      tcl <- 1
      tv <- 3.45
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl)
      v <- exp(tv + eta.v)
      linCmt() ~ add(add.sd)
    })
  }
  .fit <- function(fb) {
    .ctl <- foceiControl(print = 0, calcTables = FALSE, covMethod = "analytic", covFull = FALSE, covFallback = fb)
    suppressMessages(suppressWarnings(.nlmixr(.lin, theo_sd, "focei", .ctl)))
  }
  .s <- .fit(list(analytic = "s"))
  expect_identical(.s$covMethod, "s")
  expect_identical(.s$env$covTried$method[nrow(.s$env$covTried)], "s")
  .none <- .fit(list())
  expect_identical(.none$covMethod, "failed")
  expect_true(any(grepl("covFallback lists no fallback", .none$runInfo, fixed = TRUE)))
})

test_that("covFallback round-trips through rxUiDeparse()", {
  .ctl <- foceiControl(covFallback = list("r,s" = "s"))
  .d <- rxode2::rxUiDeparse(.ctl, "ctl")
  expect_match(deparse1(.d), "covFallback", fixed = TRUE)
  expect_identical(eval(.d[[3]])$covFallback, list("r,s" = "s"))
  expect_false(grepl("covFallback", deparse1(rxode2::rxUiDeparse(foceiControl(), "ctl")), fixed = TRUE))
})

test_that("a control without covFallback gets its control function's default", {
  expect_identical(.covFallbackOf(list()), foceiControl()$covFallback)
  expect_identical(.covFallbackOf(list(), "saemControl"), saemControl()$covFallback)
  # an empty list given by the user is no fallback, not the default
  expect_identical(.covFallbackOf(list(covFallback = list())), list())
  expect_identical(.covFallbackOf(list(covFallback = list(r = "s"))), list(r = "s"))
  expect_identical(
    saemControl()$covFallback,
    list(sa = c("linFim", "Ha"), fim = c("linFim", "Ha"), analytic = c("linFim", "Ha"), linFim = "Ha")
  )
  expect_error(saemControl(covFallback = list(sa = "r")), "cannot fall back to \"r\"")
})
