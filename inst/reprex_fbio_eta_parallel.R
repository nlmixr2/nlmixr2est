# Reprex: parallel FOCEi crashes when bioavailability f() depends on ETA
#
# Setup: nlmixr2est on parallel-v2, rxode2 on parfocei-eta-fix (or main —
# crash also reproduces on parallel-v2 alone, before the parfocei-eta-fix
# changes), rxode2ll >= 2.0.14.
#
# Symptom on cores=2:
#   "EE:[lsoda] ewt[1] = 0 <= 0." (repeated)
#   "double free or corruption (!prev)" / "corrupted size vs. prev_size"
#
# Cores=1 runs fine.  The bug is parallel-specific and predates the
# changes on parfocei-eta-fix; Matt Fidler flagged it on PR #621
# (2026-05-01): "your code in rxode2 ind_solve will corrupt any etas
# where they are used in f(), dur() or the inner problem has an issue
# with the solve and falls back to numerical differences."
#
# Run: NOT_CRAN=true Rscript repro_fbio_eta.R

suppressMessages({
  devtools::load_all(
    "/home/bill/github/nlmixr2/_worktrees-parfocei/nlmixr2est",
    quiet = TRUE, helpers = FALSE
  )
})

# 1-cmt PK with bioavailability that depends on eta.f
mod_fbio <- function() {
  ini({
    tka   <- 0.45
    tcl   <- 1
    tv    <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v  ~ 0.1
    eta.f  ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v  <- exp(tv  + eta.v)
    f(depot) <- exp(eta.f)        # <-- f() reads eta.f
    d/dt(depot)  <- -ka * depot
    d/dt(center) <-  ka * depot - cl / v * center
    cp <- center / v
    cp ~ add(add.sd)
  })
}

cat("=== cores=1 (works) ===\n")
fit1 <- nlmixr2(
  mod_fbio, nlmixr2data::theo_sd, est = "focei",
  control = foceiControl(
    print = 0, covMethod = "",
    rxControl = rxode2::rxControl(cores = 1L)
  )
)
cat("cores=1 OK, obj =", fit1$objective, "\n")

cat("=== cores=2 (crashes) ===\n")
fit2 <- nlmixr2(
  mod_fbio, nlmixr2data::theo_sd, est = "focei",
  control = foceiControl(
    print = 0, covMethod = "",
    rxControl = rxode2::rxControl(cores = 2L)
  )
)
cat("cores=2 OK, obj =", fit2$objective, "\n")
