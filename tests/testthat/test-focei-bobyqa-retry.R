nmTest({
  # Issue investigation 2026-08-23: on a real 2-compartment IV infusion
  # steady-state FOCEi fit (innerOpt="trust"), the outer bobyqa() optimizer's
  # default rhobeg=0.2 made its initial quadratic model collapse -- bobyqa
  # reported normal convergence, but 4 of 5 population parameters never moved
  # from their starting values at all. rhobeg=0.25 or larger escaped it. Root
  # cause not understood (n1qn1, FOCEi's other inner optimizer, never got
  # stuck on the same model at the same rhobeg) -- .bobyqaRetryIfStuck()
  # detects the symptom (final point never left its own starting radius) and
  # retries with a wider rhobeg, capped at 0.3, as a pragmatic safeguard.
  .fn <- function(x) sum((x - c(0.5, -0.3, 0.2))^2)
  .par <- c(0, 0, 0)
  .ctl <- list(rhobeg = 0.2, rhoend = 0.001, npt = 7)
  .distFromStart <- function(x) sqrt(sum((x - .par)^2))
  .distFromTrue <- function(x) sqrt(sum((x - c(0.5, -0.3, 0.2))^2))

  test_that(".bobyqaRetryIfStuck() retries a stuck search with a wider rhobeg", {
    skip_on_cran()
    # Fabricate a "stuck" first attempt: par unchanged from the start, as if
    # bobyqa's own search never left its initial exploration radius.
    .stuckRet <- list(par = .par, fval = .fn(.par), feval = 7L, ierr = 0L,
                      msg = "Normal exit from bobyqa")
    .ret <- suppressWarnings(
      nlmixr2est:::.bobyqaRetryIfStuck(.par, .fn, -Inf, Inf, .ctl, .stuckRet))
    # Escalated and actually moved beyond the original stuck radius, closer
    # to the true minimum than the fabricated stuck point was.
    expect_gt(.distFromStart(.ret$par), .ctl$rhobeg)
    expect_lt(.distFromTrue(.ret$par), .distFromTrue(.par))
  })

  test_that(".bobyqaRetryIfStuck() warns when it retries", {
    skip_on_cran()
    .stuckRet <- list(par = .par, fval = .fn(.par), feval = 7L, ierr = 0L,
                      msg = "Normal exit from bobyqa")
    expect_warning(
      nlmixr2est:::.bobyqaRetryIfStuck(.par, .fn, -Inf, Inf, .ctl, .stuckRet),
      "stalled")
  })

  test_that(".bobyqaRetryIfStuck() leaves an already-moved result alone", {
    skip_on_cran()
    .movedRet <- list(par = c(0.5, -0.3, 0.2), fval = 0, feval = 20L,
                      ierr = 0L, msg = "Normal exit from bobyqa")
    .ret <- nlmixr2est:::.bobyqaRetryIfStuck(.par, .fn, -Inf, Inf, .ctl, .movedRet)
    expect_identical(.ret, .movedRet)
  })

  test_that(".bobyqaRetryIfStuck() is a no-op when rhobeg was not explicitly set", {
    skip_on_cran()
    # .boundedResidOpt() (npag/saem's phi0 step) calls .bobyqa() without a
    # rhobeg, letting minqa::bobyqa() pick its own internal default -- there
    # is no known starting radius to compare against or widen, so this must
    # not engage.
    .ctlNoRhobeg <- list(rhoend = 0.001, npt = 7)
    .stuckRet <- list(par = .par, fval = .fn(.par), feval = 7L, ierr = 0L,
                      msg = "Normal exit from bobyqa")
    .ret <- nlmixr2est:::.bobyqaRetryIfStuck(.par, .fn, -Inf, Inf, .ctlNoRhobeg, .stuckRet)
    expect_identical(.ret, .stuckRet)
  })

  test_that(".bobyqaRetryIfStuck() stops escalating once rhobeg=0.3 is reached", {
    skip_on_cran()
    .stuckRet <- list(par = .par, fval = .fn(.par), feval = 7L, ierr = 0L,
                      msg = "Normal exit from bobyqa")
    .ctl03 <- list(rhobeg = 0.3, rhoend = 0.001, npt = 7)
    # Already at the escalation ceiling -- no larger candidate exists, so this
    # must return the (still-stuck) input unchanged rather than loop forever.
    .ret <- suppressWarnings(
      nlmixr2est:::.bobyqaRetryIfStuck(.par, .fn, -Inf, Inf, .ctl03, .stuckRet))
    expect_identical(.ret, .stuckRet)
  })
})
