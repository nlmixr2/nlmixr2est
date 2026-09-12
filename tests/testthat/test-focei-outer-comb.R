# The combined build: a fast=TRUE fit carries the outer-gradient sensitivities on
# the inner model, so the solve pool holds ONE model and the outer gradient and
# Hessian read the inner solve instead of swapping a peer model in.
.combModel <- function() {
  ini({ tka <- log(1.57); tcl <- log(2.72); tv <- log(31.5); add.sd <- 0.7
        eta.ka ~ 0.6; eta.cl ~ 0.3; eta.v ~ 0.1 })
  model({ ka <- exp(tka+eta.ka); cl <- exp(tcl+eta.cl); v <- exp(tv+eta.v)
          d/dt(depot) <- -ka*depot; d/dt(center) <- ka*depot - cl/v*center
          cp <- center/v; cp ~ add(add.sd) })
}

test_that("the combined build makes the inner model the outer model", {
  skip_on_cran()
  skip_if_not(.rxode2DydtCompact(), "rxode2 cannot compact a model to its leading block")
  .ui <- rxode2::rxode2(.combModel)
  .m <- function(...) {
    .u <- rxode2::rxUiDecompress(.ui)
    .u$control <- foceiControl(fast = TRUE, ...)
    suppressMessages(rxUiGet.foceiModel(list(.u)))
  }
  .comb <- .m()
  expect_true(.comb$outerComb)
  .md5 <- function(x) rxode2::rxModelVars(x)$md5[["parsed_md5"]]
  expect_identical(.md5(.comb$inner), .md5(.comb$outer))
  .lhs <- rxode2::rxModelVars(.comb$inner)$lhs
  # the FOCEi block is untouched at the front ...
  expect_identical(.lhs[1:8], c("rx_pred_", paste0("rx__sens_rx_pred__BY_ETA_", 1:3, "___"),
                                "rx_r_", paste0("rx__sens_rx_r__BY_ETA_", 1:3, "___")))
  # ... and the outer columns follow it
  expect_true(all(c("rx_predf_", "rx_f1_ETA_1_", "rx_f2_ETA_1__ETA_1_", "rx_rvarf_",
                    "rx_rsig_4_") %in% .lhs[-(1:8)]))
  expect_equal(length(rxode2::rxModelVars(.comb$inner)$state), 26L)
  # the column map resolves by name in the combined model
  .cols <- .vaeOuterCols(c(list(augMod = .comb$inner), .comb$outerMeta))
  expect_equal(.lhs[.cols$predf + 1L], "rx_predf_")
  expect_identical(.cols$lhsNames, .lhs)
  # opting out keeps the separate outer model and the 8-state inner
  .sep <- .m(outerCombine = FALSE)
  expect_false(.sep$outerComb)
  expect_false(identical(.md5(.sep$inner), .md5(.sep$outer)))
  expect_equal(length(rxode2::rxModelVars(.sep$inner)$state), 8L)
})

test_that("the outer gradient and Hessian read the inner solve, and agree with the separate model", {
  skip_on_cran()
  skip_if_not(.rxode2DydtCompact(), "rxode2 cannot compact a model to its leading block")
  .run <- function(outerCombine) {
    observed <- NULL
    optimizer <- function(par, fn, gr, lower, upper, control) {
      fn(par)
      g <- gr(par)
      h <- control$hessian(par)
      i <- .odeSwapInfo()
      # central differences of the objective, the reference both builds must match
      ref <- vapply(seq_along(par), function(j) {
        step <- 1e-4 * max(1, abs(par[j]))
        plus <- minus <- par; plus[j] <- plus[j] + step; minus[j] <- minus[j] - step
        (fn(plus) - fn(minus)) / (2 * step)
      }, numeric(1))
      fn(par)
      observed <<- list(gradient = g, hessian = h, reference = ref, poolName = i$poolName,
                        loaded = i$models$name[i$models$loaded])
      list(x = par, convergence = 0L, message = "check")
    }
    fit <- .nlmixr(.combModel, nlmixr2data::theo_sd, "focei", foceiControl(
      fast = TRUE, outerCombine = outerCombine, outerOpt = optimizer, innerOpt = "n1qn1",
      print = 0, covMethod = "", calcTables = FALSE, maxInnerIterations = 1000L, epsilon = 1e-10,
      rxControl = rxode2::rxControl(atol = 1e-10, rtol = 1e-10)))
    list(fit = fit, obs = observed)
  }
  .c <- .run(TRUE)
  .s <- .run(FALSE)
  nsub <- 12L; neta <- 3L
  expect_true(.c$fit$env$outerComb)
  expect_false(.s$fit$env$outerComb)
  expect_equal(.c$obs$poolName, "inner")
  expect_equal(.s$obs$poolName, "outer")
  # combined: the gradient's outer solve and the Hessian's base solve are the inner
  # solve; only the 4*neta third-order probes integrate
  expect_equal(.c$fit$env$nOuterSolveReused, 2L * nsub)
  expect_equal(.c$fit$env$nOuterSolveRun, 4L * neta * nsub)
  # ... and the inner iterations ran compacted to the inner block; only the two
  # derivative passes (gradient, Hessian) solved every subject at full width
  expect_equal(.c$fit$env$combInnerNeq, 8L)
  expect_equal(.c$fit$env$nInnerSolveFull, 2L * nsub)
  expect_gt(.c$fit$env$nInnerSolveCompact, .c$fit$env$nInnerSolveFull)
  expect_equal(.s$fit$env$nOuterSolveReused, 0L)
  expect_gt(.c$fit$env$nAnalyticGradDirect, 0L)
  # both builds differentiate the same objective; the two solves are of different ODE
  # systems (26 vs 8 + 18 states) so they agree only to solver tolerance
  .rel <- function(a, b) max(abs(a - b)) / max(abs(b))
  expect_lt(.rel(.c$obs$gradient, .c$obs$reference), 1e-3)
  expect_lt(.rel(.s$obs$gradient, .s$obs$reference), 1e-3)
  expect_lt(.rel(.c$obs$gradient, .s$obs$gradient), 1e-3)
  expect_lt(norm(.c$obs$hessian - .s$obs$hessian, "F") / norm(.s$obs$hessian, "F"), 1e-2)
})

test_that("a combined fast fit reaches the separate build's optimum", {
  skip_on_cran()
  skip_if_not(.rxode2DydtCompact(), "rxode2 cannot compact a model to its leading block")
  .fit <- function(...) .nlmixr(.combModel, nlmixr2data::theo_sd, "focei",
    foceiControl(fast = TRUE, print = 0, covMethod = "", calcTables = FALSE, ...))
  .c <- .fit(outerCombine = TRUE)
  .s <- .fit(outerCombine = FALSE)
  expect_true(.c$env$outerComb)
  expect_equal(.c$objf, .s$objf, tolerance = 1e-4)
  expect_equal(unname(fixef(.c)), unname(fixef(.s)), tolerance = 1e-2)
  expect_gt(.c$env$nOuterSolveReused, 0L)
})
