# saem two-level IOV with subjects that have different numbers of occasions
# (#1504).  An absent occasion is an uninformative eta, so it must not count
# as a zero deviation in the pooled occasion variance.

# Simulated one-compartment oral data, occasion = dose number, true IIV and
# IOV variance 0.1 on V; `nOcc` is the number of occasions per subject.
.iovUnbalancedData <- function(nOcc, seed) {
  withr::local_seed(seed)
  .n <- length(nOcc)
  .etaCl <- rnorm(.n, 0, sqrt(0.1))
  .etaV <- rnorm(.n, 0, sqrt(0.1))
  .tad <- c(0.5, 1, 2, 4, 8, 12, 23.5)
  .rows <- lapply(seq_len(.n), function(i) {
    .occ <- seq_len(nOcc[i])
    .iov <- rnorm(length(.occ), 0, sqrt(0.1))
    .dose <- data.frame(
      ID = i,
      TIME = (.occ - 1) * 24,
      AMT = 320,
      EVID = 1,
      DV = NA,
      occ = .occ,
      IOVV = .iov
    )
    .obs <- data.frame(
      ID = i,
      TIME = rep((.occ - 1) * 24, each = length(.tad)) + .tad,
      AMT = 0,
      EVID = 0,
      DV = NA,
      occ = rep(.occ, each = length(.tad)),
      IOVV = rep(.iov, each = length(.tad))
    )
    .d <- rbind(.dose, .obs)
    .d$ETACL <- .etaCl[i]
    .d$ETAV <- .etaV[i]
    .d[order(.d$TIME, -.d$EVID), ]
  })
  .d <- do.call(rbind, .rows)
  # nolint start: object_usage_linter. rxode2 model block (NSE).
  .sim <- rxode2::rxode2({
    ka <- 1.5
    cl <- 2.5 * exp(ETACL)
    v <- 30 * exp(ETAV + IOVV)
    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - cl / v * central
    cp <- central / v
  })
  # nolint end
  .s <- rxode2::rxSolve(.sim, .d, returnType = "data.frame")
  .o <- .d$EVID == 0
  .d$DV[.o] <- .s$cp * (1 + rnorm(sum(.o), 0, 0.1)) + rnorm(sum(.o), 0, 0.1)
  .d[, c("ID", "TIME", "AMT", "EVID", "DV", "occ")]
}

# nolint start: object_usage_linter. rxode2 ini()/model() are NSE blocks.
.iovUnbalancedModel <- function() {
  ini({
    tka <- log(1.5)
    tcl <- log(2.5)
    tv <- log(30)
    eta.cl ~ 0.1
    eta.v ~ 0.1
    iov.v ~ 0.1 | occ
    prop.sd <- 0.1
    add.sd <- 0.1
  })
  model({
    ka <- exp(tka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v + iov.v)
    linCmt() ~ add(add.sd) + prop(prop.sd)
  })
}
# nolint end

.iovUnbalancedOcc <- function() {
  withr::local_seed(2)
  sample(1:8, 60, replace = TRUE, prob = 0.6^(0:7))
}

test_that("an unbalanced occasion design does not shrink the two-level IOV variance (#1504)", {
  skip_on_cran()
  .nOcc <- .iovUnbalancedOcc()
  # most occasion slots are empty: that is the design the bias grew with
  expect_lt(sum(.nOcc) / (length(.nOcc) * max(.nOcc)), 0.5)
  .d <- .iovUnbalancedData(.nOcc, seed = 1)
  .f <- suppressWarnings(nlmixr2(
    .iovUnbalancedModel(),
    .d,
    est = "saem",
    control = saemControl(
      nBurn = 150,
      nEm = 100,
      seed = 42L,
      print = 0L,
      covMethod = "",
      calcTables = FALSE,
      iovMethod = "twoLevel"
    )
  ))
  # truth 0.1; pooling over all subjects gave about 0.023 here
  .iov <- .f$omega$occ["iov.v", "iov.v"]
  expect_gt(.iov, 0.06)
  expect_lt(.iov, 0.16)
})

test_that("the collapsed form samples absent occasions instead of freezing them (#1504)", {
  # phi columns: tka (not pooled), then three collapsed occasion columns
  .ue <- cbind(a = c(1, 1), b = c(1, 0), c = c(0, 1), d = c(1, 1))
  .i1 <- 1:4
  .pool <- c(0L, 1L, 1L, 1L)
  .r <- .saemUeCollapsedPool(.ue, c(1L, 2L), .i1, .pool, 1L)
  expect_true(all(.r$ue[, 2:4] == 1))
  expect_equal(.r$revisitCols, integer(0))
  # the two-level form keeps the mask: its occasion columns are independent
  .r <- .saemUeCollapsedPool(.ue, c(1L, 2L), .i1, .pool, 0L)
  expect_equal(.r$ue, .ue)
  expect_equal(.r$revisitCols, c(1L, 2L))
  # a non-pooled column keeps its own verdict
  .ue[1, 1] <- 0
  .r <- .saemUeCollapsedPool(.ue, 0L, .i1, .pool, 1L)
  expect_equal(unname(.r$ue[1, 1]), 0)
  expect_equal(.r$revisitCols, 0L)
})
