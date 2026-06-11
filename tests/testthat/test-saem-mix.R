nmTest({
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
  })
})
