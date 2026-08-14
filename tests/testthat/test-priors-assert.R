nmTest({

  ## A prior an estimation method cannot use must be refused rather than
  ## quietly dropped, otherwise the fit does something other than what
  ## the model says with nothing to tell the user.

  .hasPriors <- function() {
    ## the `prior` column arrived with a newer lotri
    "prior" %in% names(rxode2::.rxBlankIni("empty")) ||
      exists("lotriPriorDists", envir=asNamespace("lotri"), inherits=FALSE)
  }

  .hasRxAsserts <- function() {
    exists("assertRxUiNormalPriors", envir=asNamespace("rxode2"), inherits=FALSE)
  }

  .mod <- function(prior=NULL) {
    .ini <- c("tka <- 0.45", "tcl <- 1", "tv <- 3.45", "eta.ka ~ 0.6",
              "add.sd <- 0.7", prior)
    .txt <- paste0("function() {\n ini({\n", paste(.ini, collapse="\n"),
                   "\n})\n model({\n ka <- exp(tka + eta.ka)\n cl <- exp(tcl)\n",
                   " v <- exp(tv)\n linCmt() ~ add(add.sd)\n})\n}")
    eval(str2lang(.txt))()
  }

  .env <- function(cls, ui) {
    .e <- new.env(parent=emptyenv())
    assign("ui", ui, envir=.e)
    class(.e) <- c(cls, "nlmixr2Est")
    .e
  }

  test_that("a method declares what priors it supports", {
    ## an ordinary method has no declaration, which means none
    expect_equal(.nlmixr2PriorSupport(.env("focei", NULL)), "none")
    ## and so does a class with no method at all
    expect_equal(.nlmixr2PriorSupport(.env("notAMethod", NULL)), "none")
  })

  test_that("an unknown support level is an error", {
    nlmixr2Est.fakeBadLevel <- function(env, ...) TRUE
    attr(nlmixr2Est.fakeBadLevel, "nlmixr2Priors") <- "sometimes"
    registerS3method("nlmixr2Est", "fakeBadLevel", nlmixr2Est.fakeBadLevel,
                     envir=globalenv())
    expect_error(.nlmixr2PriorSupport(.env("fakeBadLevel", NULL)),
                 "nlmixr2Priors")
  })

  test_that("a model without priors is accepted by every method", {
    .ui <- .mod()
    expect_error(.nlmixr2AssertPriors(.env("focei", .ui)), NA)
    expect_error(.nlmixr2AssertPriors(.env("saem", .ui)), NA)
  })

  test_that("a prior is refused by a method that cannot use it", {
    skip_if_not(.hasPriors())
    .ui <- .mod("prior(tka) ~ dnorm(0, 10)")

    ## the message names the parameter and the method, so the user can
    ## see which prior and which est
    expect_error(.nlmixr2AssertPriors(.env("focei", .ui)), "tka")
    expect_error(.nlmixr2AssertPriors(.env("focei", .ui)), "focei")
    expect_error(.nlmixr2AssertPriors(.env("saem", .ui)), "saem")

    ## a method registered by another package gets the same treatment,
    ## which is the point of checking in the generic
    expect_error(.nlmixr2AssertPriors(.env("nonmem", .ui)), "nonmem")
  })

  test_that("a method that declares support is not blocked", {
    skip_if_not(.hasPriors())
    .ui <- .mod("prior(tka) ~ dnorm(0, 10)")

    nlmixr2Est.fakeAll <- function(env, ...) TRUE
    attr(nlmixr2Est.fakeAll, "nlmixr2Priors") <- "all"
    registerS3method("nlmixr2Est", "fakeAll", nlmixr2Est.fakeAll, envir=globalenv())

    expect_equal(.nlmixr2PriorSupport(.env("fakeAll", .ui)), "all")
    expect_error(.nlmixr2AssertPriors(.env("fakeAll", .ui)), NA)
  })

  test_that("a normal-only method takes a normal prior but not others", {
    skip_if_not(.hasPriors())
    skip_if_not(.hasRxAsserts())

    nlmixr2Est.fakeNormal <- function(env, ...) TRUE
    attr(nlmixr2Est.fakeNormal, "nlmixr2Priors") <- "normal"
    registerS3method("nlmixr2Est", "fakeNormal", nlmixr2Est.fakeNormal,
                     envir=globalenv())

    expect_error(.nlmixr2AssertPriors(
      .env("fakeNormal", .mod("prior(tka) ~ dnorm(0, 10)"))), NA)
    expect_error(.nlmixr2AssertPriors(
      .env("fakeNormal", .mod("prior(tka) ~ dgamma(2, 1)"))))
  })

  test_that("an nwpri method takes omega degrees of freedom", {
    skip_if_not(.hasPriors())
    skip_if_not(.hasRxAsserts())

    nlmixr2Est.fakeNwpri <- function(env, ...) TRUE
    attr(nlmixr2Est.fakeNwpri, "nlmixr2Priors") <- "nwpri"
    registerS3method("nlmixr2Est", "fakeNwpri", nlmixr2Est.fakeNwpri,
                     envir=globalenv())

    ## degrees of freedom on the omega are what NWPRI needs
    expect_error(.nlmixr2AssertPriors(
      .env("fakeNwpri", .mod("prior(eta.ka) ~ invWishart(2)"))), NA)
    ## a normal prior on the omega values is TNPRI, which it does not do
    expect_error(.nlmixr2AssertPriors(
      .env("fakeNwpri", .mod("om.eta.ka ~ 0.01"))))
  })

})
