test_that("the outer Hessian native interface reports scope and restores probes", {
  skip_on_cran()
  api <- new.env(parent = globalenv())
  header <- system.file("include", "nlmixr2estFoceiPtr.h", package = "nlmixr2est")
  Rcpp::sourceCpp(code = paste0('
#include <Rcpp.h>
#include "', gsub("\\\\", "/", header), '"
iniNlmixr2estFoceiGlobals
// [[Rcpp::export]]
Rcpp::List probeOuterHessian(Rcpp::List ptr, Rcpp::NumericVector x, double step) {
  iniNlmixr2estFocei(ptr);
  if (!nlmixr2FoceiOuterHessianP) return Rcpp::List::create(Rcpp::Named("status")=-99);
  Rcpp::NumericMatrix H(x.size(), x.size());
  int rc = nlmixr2FoceiOuterHessianP(x.begin(), x.size(), step, H.begin());
  return Rcpp::List::create(Rcpp::Named("status")=rc, Rcpp::Named("hessian")=H);
}
'), env = api, showOutput = FALSE)
  ptr <- .nlmixr2estFoceiPtrs()
  expect_equal(api$probeOuterHessian(ptr[1:10], c(0, 0), 1e-3)$status, -99L)
  expect_equal(api$probeOuterHessian(ptr, c(0, 0), 1e-3)$status, -1L)
  model <- function() {
    ini({ level <- 0.2; error <- fix(0.5); etaLevel ~ 0.2 })
    model({ prediction <- exp(level+etaLevel); prediction ~ add(error) })
  }
  data <- data.frame(ID = c(1, 1, 2, 2), TIME = c(1, 2, 1, 2),
                     DV = c(1.5, 1.7, 0.8, 1), AMT = 0, EVID = 0)
  for (fast in c(TRUE, FALSE)) {
    observed <- NULL
    optimizer <- function(par, fn, gr, lower, upper, control) {
      fn(par); before <- fn(par)
      result <- api$probeOuterHessian(ptr, par, 1e-3)
      badShape <- api$probeOuterHessian(ptr, par[-1], 1e-3)$status
      badStep <- api$probeOuterHessian(ptr, par, 0)$status
      failed <- api$probeOuterHessian(ptr, par, 1e100)$status
      observed <<- list(result = result, badShape = badShape, badStep = badStep,
                        failed = failed, before = before, after = fn(par),
                        wrapper = if (fast) control$hessian(par) else NULL)
      list(x = par, convergence = 0L, message = "Outer Hessian API check")
    }
    .nlmixr(model, data, "focei", control = foceiControl(
      fast = fast, outerOpt = optimizer, innerOpt = "n1qn1", epsilon = 1e-10,
      maxInnerIterations = 1000L, print = 0, covMethod = "", calcTables = FALSE))
    expect_equal(observed$result$status, if (fast) 0L else -4L)
    expect_equal(observed$badShape, -2L)
    expect_equal(observed$badStep, -2L)
    expect_true(observed$failed %in% c(-3L, -4L))
    expect_equal(observed$before, observed$after, tolerance = 1e-7)
    if (fast) expect_equal(observed$result$hessian, observed$wrapper, tolerance = 1e-6)
  }
  expect_equal(api$probeOuterHessian(ptr, c(0, 0), 1e-3)$status, -1L)
})
