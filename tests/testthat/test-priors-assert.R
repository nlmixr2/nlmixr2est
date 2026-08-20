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
    ## an ordinary method with no declaration means none; saem does not
    ## declare nlmixr2Priors (focei's family does, as of #931, so it is no
    ## longer a representative "undeclared" example)
    expect_equal(.nlmixr2PriorSupport(.env("saem", NULL)), "none")
    ## and so does a class with no method at all
    expect_equal(.nlmixr2PriorSupport(.env("notAMethod", NULL)), "none")
    ## focei's family declares "general" support (#931): a prior on a
    ## population parameter AND on an omega element, either directly
    ## (the "tnpri" convention) or via invWishart() degrees of freedom
    ## (the "nwpri" convention) -- foceiPriorOmegaGradAdd()'s chain-rule
    ## into op_focei.cholOmegaInv (src/inner.cpp) is the same regardless
    ## of which convention built the term.
    expect_equal(.nlmixr2PriorSupport(.env("focei", NULL)), "general")
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
    ## see which prior and which est.  saem declares no support at all
    ## (still "none"); focei's family now declares "general" (#931), so a
    ## refusal for it needs something outside even that: a joint block
    ## with no om.<eta> member reachable is model-shape refused elsewhere,
    ## so saem alone demonstrates this.
    expect_error(.nlmixr2AssertPriors(.env("saem", .ui)), "tka")
    expect_error(.nlmixr2AssertPriors(.env("saem", .ui)), "saem")

    ## a method registered by another package gets the same treatment,
    ## which is the point of checking in the generic
    expect_error(.nlmixr2AssertPriors(.env("nonmem", .ui)), "nonmem")
  })

  test_that("focei's family accepts a prior on omega too (#931)", {
    skip_if_not(.hasPriors())
    expect_error(.nlmixr2AssertPriors(.env("focei", .mod("om.eta.ka ~ 0.01"))), NA)
    expect_error(.nlmixr2AssertPriors(
      .env("focei", .mod("prior(eta.ka) ~ invWishart(2)"))), NA)
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

  test_that("a tnpri method takes a normal prior (incl. on omega) but not others", {
    skip_if_not(.hasPriors())
    skip_if_not(.hasRxAsserts())

    nlmixr2Est.fakeTnpri <- function(env, ...) TRUE
    attr(nlmixr2Est.fakeTnpri, "nlmixr2Priors") <- "tnpri"
    registerS3method("nlmixr2Est", "fakeTnpri", nlmixr2Est.fakeTnpri,
                     envir=globalenv())

    expect_error(.nlmixr2AssertPriors(
      .env("fakeTnpri", .mod("prior(tka) ~ dnorm(0, 10)"))), NA)
    ## a normal prior directly on omega is exactly what tnpri is for
    expect_error(.nlmixr2AssertPriors(
      .env("fakeTnpri", .mod("om.eta.ka ~ 0.01"))), NA)
    expect_error(.nlmixr2AssertPriors(
      .env("fakeTnpri", .mod("prior(tka) ~ dgamma(2, 1)"))))
    ## nwpri's own mechanism (omega degrees of freedom) is refused
    expect_error(.nlmixr2AssertPriors(
      .env("fakeTnpri", .mod("prior(eta.ka) ~ invWishart(2)"))))
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

  test_that("a theta-only method takes a theta prior but not one on omega", {
    skip_if_not(.hasPriors())

    nlmixr2Est.fakeTheta <- function(env, ...) TRUE
    attr(nlmixr2Est.fakeTheta, "nlmixr2Priors") <- "theta"
    registerS3method("nlmixr2Est", "fakeTheta", nlmixr2Est.fakeTheta,
                     envir=globalenv())

    ## a Cauchy is fine too -- "theta" restricts WHERE the prior can sit
    ## (never on omega), not the distribution family
    expect_error(.nlmixr2AssertPriors(
      .env("fakeTheta", .mod("prior(tka) ~ dnorm(0, 10)"))), NA)
    expect_error(.nlmixr2AssertPriors(
      .env("fakeTheta", .mod("prior(tka) ~ dcauchy(0, 5)"))), NA)
    expect_error(.nlmixr2AssertPriors(
      .env("fakeTheta", .mod("om.eta.ka ~ 0.01"))), "omega")
    expect_error(.nlmixr2AssertPriors(
      .env("fakeTheta", .mod("prior(eta.ka) ~ invWishart(2)"))), "omega")
  })

  test_that("a general method takes anything the kernel supports", {
    skip_if_not(.hasPriors())

    nlmixr2Est.fakeGeneral <- function(env, ...) TRUE
    attr(nlmixr2Est.fakeGeneral, "nlmixr2Priors") <- "general"
    registerS3method("nlmixr2Est", "fakeGeneral", nlmixr2Est.fakeGeneral,
                     envir=globalenv())

    expect_error(.nlmixr2AssertPriors(
      .env("fakeGeneral", .mod("prior(tka) ~ dcauchy(0, 5)"))), NA)
    expect_error(.nlmixr2AssertPriors(
      .env("fakeGeneral", .mod("om.eta.ka ~ 0.01"))), NA)
    expect_error(.nlmixr2AssertPriors(
      .env("fakeGeneral", .mod("prior(eta.ka) ~ invWishart(2)"))), NA)
  })

  test_that(".nlmixr2PriorMethod() reads the omega convention off the ini() syntax", {
    skip_if_not(.hasPriors())

    ## no omega-referencing prior at all: the three methods agree, so this
    ## is a harmless default
    expect_equal(.nlmixr2PriorMethod(.mod("prior(tka) ~ dnorm(0, 10)")), "general")
    ## a normal directly on omega only makes sense as tnpri
    expect_equal(.nlmixr2PriorMethod(.mod("om.eta.ka ~ 0.01")), "tnpri")
    ## omega degrees of freedom only makes sense as nwpri
    expect_equal(.nlmixr2PriorMethod(.mod("prior(eta.ka) ~ invWishart(2)")), "nwpri")
    ## neither nwpri nor tnpri has a Cauchy analogue
    expect_equal(.nlmixr2PriorMethod(.mod("prior(tka) ~ dcauchy(0, 5)")), "general")
  })

  test_that(".nlmixr2BuildPriorSpec() is a no-op for a model with no priors", {
    expect_null(.nlmixr2BuildPriorSpec(.mod()))
  })

  test_that(".nlmixr2BuildPriorSpec() returns a usable external pointer", {
    skip_if_not(.hasPriors())
    skip_if_not(exists("rxPriorBuildSpec", envir=asNamespace("rxode2"), inherits=FALSE))

    .spec <- .nlmixr2BuildPriorSpec(.mod("prior(tka) ~ dnorm(0, 10)"))
    expect_true(is(.spec, "externalptr"))
  })

})
