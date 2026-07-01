nmTest({
  test_that("SAEM does not clamp non-mu-referenced omega below 1.0 (regression for accidental 1.0 floor)", {
    set.seed(100)
    n_subj <- 20
    d <- do.call(rbind, lapply(1:n_subj, function(i) {
      times <- c(0.5, 1, 2, 4, 8, 12, 24)
      data.frame(
        ID = i,
        TIME = c(0, times),
        AMT = c(100, rep(0, length(times))),
        EVID = c(1, rep(0, length(times))),
        DV = c(0, rep(1, length(times))),
        CMT = c(1, rep(2, length(times)))
      )
    }))

    one.compartment.i0 <- function() {
      ini({
        tka <- log(1.5)
        tcl <- log(2.0)
        tv <- log(20)
        eta.cl ~ 0.01
        # non-mu-referenced ETA (add.err on tka, not linear in the parameter):
        eta.ka2 ~ 0.09
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka) + eta.ka2
        cl <- exp(tcl + eta.cl)
        v <- exp(tv)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    fit <- suppressWarnings(.nlmixr(
      one.compartment.i0, d, est = "saem",
      saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                  calcTables = FALSE, covMethod = 0L)
    ))
    # eta.ka2 is not mu-referenced (added, not multiplied); its estimated
    # omega must not be silently floored up to 1.0 by the mixture-branch
    # omega-floor regression.
    .om <- fit$omega
    expect_true(.om["eta.ka2", "eta.ka2"] < 0.99)
  })

  test_that("test SAEM mixture model estimation", {

    set.seed(42)
    n_subj <- 30
    sub_pop <- rbinom(n_subj, 1, 0.6) + 1 # 1 or 2
    cl_sim <- ifelse(sub_pop == 1, 1.2, 6.0)

    sim_data <- do.call(rbind, lapply(1:n_subj, function(i) {
      subj_cl <- cl_sim[i]
      times <- c(0.5, 1, 2, 4, 8, 12, 24)
      ka_val <- 1.5
      v_val <- 24.0
      k_val <- subj_cl / v_val
      cp <- 100 * ka_val / (v_val * (ka_val - k_val)) * (exp(-k_val * times) - exp(-ka_val * times)) + rnorm(length(times), 0, 0.05)
      cp[cp < 0] <- 0
      data.frame(
        ID = i,
        TIME = c(0, times),
        AMT = c(100, rep(0, length(times))),
        EVID = c(1, rep(0, length(times))),
        DV = c(0, cp),
        CMT = c(1, rep(2, length(times)))
      )
    }))

    one.compartment.mix <- function() {
      ini({
        tka <- log(1.5)
        tcl1 <- log(1.0)
        tcl2 <- log(5.0)
        tv <- log(20)
        p1 <- 0.5
        eta.cl ~ 0.01
        eta.v ~ 0.01
        eta.ka ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- mix(exp(tcl1 + eta.cl), p1, exp(tcl2 + eta.cl))
        v <- exp(tv + eta.v)
        d/dt(depot) <- -ka * depot
        d/dt(center) <- ka * depot - cl / v * center
        cp <- center / v
        cp ~ add(add.sd)
      })
    }

    # Run SAEM estimation with mixture
    fit_saem <- expect_error(
      .nlmixr(one.compartment.mix, sim_data, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = 0L)),
      NA
    )

    # Check that estimated mixture parameter p1 exists
    expect_true("p1" %in% names(fit_saem$theta))

    # Check that estimated parameters are accessible
    expect_true(is.numeric(fit_saem$theta["p1"]))

    # Test SAEM mixture model with split ETAs and covariance calculation
    twoPopSplit <- function() {
      ini({
        tka   <- log(1.5)
        tcl1  <- log(1.0)
        tcl2  <- log(5.0)
        tv    <- log(20)
        p1    <- 0.5
        eta.cl1 ~ 0.01
        eta.cl2 ~ 0.01
        eta.v  ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl1), p1, exp(tcl2 + eta.cl2))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit_saem_split <- expect_error(
      .nlmixr(twoPopSplit, sim_data, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = "linFim")),
      NA
    )
    expect_true("p1" %in% names(fit_saem_split$theta))
    expect_true(!is.null(fit_saem_split$cov))
    expect_true("eta.cl" %in% rownames(fit_saem_split$omega))
    expect_true(!"eta.cl1" %in% rownames(fit_saem_split$omega))
    expect_true("eta.cl" %in% names(fit_saem_split$ranef))
    expect_true(!"eta.cl1" %in% names(fit_saem_split$ranef))

    # Test SAEM mixture model with 3 split ETAs and covariance calculation
    threePopSplit <- function() {
      ini({
        tka   <- log(1.5)
        tcl1  <- log(1.0)
        tcl2  <- log(3.0)
        tcl3  <- log(6.0)
        tv    <- log(20)
        p1    <- 0.3
        p2    <- 0.4
        eta.cl1 ~ 0.01
        eta.cl2 ~ 0.01
        eta.cl3 ~ 0.01
        eta.v  ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl1), p1, exp(tcl2 + eta.cl2), p2, exp(tcl3 + eta.cl3))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit_saem_split3 <- expect_error(
      .nlmixr(threePopSplit, sim_data, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = "linFim")),
      NA
    )
    expect_true("p1" %in% names(fit_saem_split3$theta))
    expect_true("p2" %in% names(fit_saem_split3$theta))
    expect_true(!is.null(fit_saem_split3$cov))
    expect_true("eta.cl" %in% rownames(fit_saem_split3$omega))
    expect_true(!"eta.cl1" %in% rownames(fit_saem_split3$omega))
    expect_true("eta.cl" %in% names(fit_saem_split3$ranef))
    expect_true(!"eta.cl1" %in% names(fit_saem_split3$ranef))

    # Test SAEM mixture model with split ETAs and bounded parameters
    # to verify that back-transformations are correctly applied (not NaN)
    twoPopSplitBounded <- function() {
      ini({
        tka   <- log(1.5)
        tcl1  <- log(c(0.01, 1.0, 100))  # Bounded!
        tcl2  <- log(c(0.01, 5.0, 100))  # Bounded!
        tv    <- log(20)
        p1    <- 0.5
        eta.cl1 ~ 0.01
        eta.cl2 ~ 0.01
        eta.v  ~ 0.01
        add.sd <- 0.05
      })
      model({
        ka <- exp(tka)
        cl <- mix(exp(tcl1 + eta.cl1), p1, exp(tcl2 + eta.cl2))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit_saem_split_bounded <- expect_error(
      .nlmixr(twoPopSplitBounded, sim_data, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = "linFim")),
      NA
    )
    # Test SAEM mixture model with nested expit/logit functions under mix()
    # (non-split ETAs, so covariance cannot be calculated; covMethod=0L)
    twoPopNestedExpit <- function() {
      ini({
        tka   <- log(1.2)
        tcl1  <- logit(1.0, 0.1, 200)
        tcl2  <- logit(5.0, 0.1, 200)
        tv    <- log(30)
        p1    <- 0.67
        eta.cl ~ 0.09
        eta.v  ~ 0.04
        add.sd <- 0.2
      })
      model({
        ka <- exp(tka)
        cl <- mix(expit(tcl1 + eta.cl, 0.1, 200), p1, expit(tcl2 + eta.cl, 0.1, 200))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit_saem_nested <- expect_error(
      .nlmixr(twoPopNestedExpit, sim_data, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = 0L)),
      NA
    )
    expect_true("p1" %in% names(fit_saem_nested$theta))

    # Test SAEM mixture model with nested expit/logit functions under mix() and split ETAs
    # (split ETAs, so covariance is supported; covMethod="linFim")
    twoPopNestedExpitSplit <- function() {
      ini({
        tka   <- log(1.2)
        tcl1  <- logit(1.0, 0.1, 200)
        tcl2  <- logit(5.0, 0.1, 200)
        tv    <- log(30)
        p1    <- 0.67
        eta.cl1 ~ 0.09
        eta.cl2 ~ 0.09
        eta.v  ~ 0.04
        add.sd <- 0.2
      })
      model({
        ka <- exp(tka)
        cl <- mix(expit(tcl1 + eta.cl1, 0.1, 200), p1, expit(tcl2 + eta.cl2, 0.1, 200))
        v  <- exp(tv + eta.v)
        linCmt() ~ add(add.sd)
      })
    }

    fit_saem_nested_split <- expect_error(
      .nlmixr(twoPopNestedExpitSplit, sim_data, est="saem",
              saemControl(print = 0, seed = 1234, nBurn = 5, nEm = 5,
                          calcTables = FALSE, covMethod = "linFim")),
      NA
    )
    expect_true("p1" %in% names(fit_saem_nested_split$theta))
    expect_true(!is.null(fit_saem_nested_split$cov))
  })
})


