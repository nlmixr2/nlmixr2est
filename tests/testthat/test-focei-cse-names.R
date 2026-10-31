test_that('repeated outer CSE passes preserve existing temporaries', {
  cse <- getFromNamespace('.rxOptExprC', 'rxode2')
  x <- 'rx_expr_0~a+b\ny=rx_expr_0*(c+d)+(c+d)'
  evaluate <- function(txt) {
    env <- list2env(list(a = 2, b = 3, c = 7, d = 11))
    eval(parse(text = gsub('~', '<-', txt, fixed = TRUE)), env)
    env$y
  }
  for (pass in 1:3) {
    out <- cse(.foceiCseRename(x, pass))
    if (is.na(out)) {
      expect_gt(pass, 1L)
      break
    }
    x <- out
    expect_equal(evaluate(x), 108)
  }
})

test_that('outer CSE preserves a small ODE and its sensitivity', {
  cse <- getFromNamespace('.rxOptExprC', 'rxode2')
  original <- paste(
    'rx_expr_0~a+b',
    'k=rx_expr_0*(c+d)+(c+d)',
    'd/dt(A)=-k*A',
    'd/dt(S)=-k*S-A',
    'A(0)=10',
    'S(0)=0',
    'pred=A',
    sep = '\n'
  )
  optimized <- original
  for (pass in 1:3) {
    out <- cse(.foceiCseRename(optimized, pass))
    if (is.na(out)) {
      break
    }
    optimized <- out
  }
  params <- c(a = .02, b = .03, c = .07, d = .11)
  events <- data.frame(time = 0:24, evid = 0, amt = 0)
  solve <- function(txt) {
    as.data.frame(rxode2::rxSolve(
      rxode2::rxode2(txt),
      params = params,
      events = events
    ))
  }
  ref <- solve(original)
  actual <- solve(optimized)
  expect_true(all(is.finite(as.matrix(actual))))
  expect_equal(actual[c('A', 'S', 'pred')], ref[c('A', 'S', 'pred')], tolerance = 1e-8)
})


test_that('outer CSE renaming also avoids existing target prefixes', {
  text <- 'rx_focei_cse1_0~17\nrx_expr_0~a+b\ny=rx_expr_0*(c+d)+(c+d)+rx_focei_cse1_0'
  env <- list2env(list(a = 2, b = 3, c = 7, d = 11))
  renamed <- .foceiCseRename(text, 1L)
  eval(parse(text = gsub('~', '<-', renamed, fixed = TRUE)), env)
  expect_equal(env$y, 125)
})

test_that('outer CSE renames only defined temporary identifiers', {
  text <- 'rx_expr_1~a+b\ny=rx_expr_1+rx_expr_10+my_rx_expr_1+rx_expr_1_extra+my.rx_expr_1'
  expect_identical(
    .foceiCseRename(text, 1L),
    'rx_focei_cse1_1~a+b\ny=rx_focei_cse1_1+rx_expr_10+my_rx_expr_1+rx_expr_1_extra+my.rx_expr_1'
  )
  expect_identical(.foceiCseRename('y=rx_expr_1+my_rx_expr_1', 1L), 'y=rx_expr_1+my_rx_expr_1')
})

test_that('large-model CSE preserves covariate names in a fast fit', {
  withr::local_options(rxode2.optExprChunkChars = 1L)
  model <- function() {
    ini({ level <- .211123; error <- .2; eta ~ .3 })
    model({ d/dt(dummy) <- 0; dummy(0) <- 1
            prediction <- exp(level + eta) * dummy + my_rx_expr_testcov
            prediction ~ add(error) })
  }
  data <- data.frame(ID = 1, TIME = 1:3, DV = c(1.6, 1.7, 1.8), my_rx_expr_testcov = c(.1, .2, .3), AMT = 0, EVID = 0)
  fit <- function(optimize) {
    suppressWarnings(suppressMessages(nlmixr2(
      model,
      data,
      'focei',
      foceiControl(
        fast = TRUE,
        optExpression = optimize,
        maxOuterIterations = 0,
        print = 0,
        covMethod = '',
        calcTables = FALSE
      )
    )))
  }
  reference <- fit(FALSE)
  optimized <- fit(TRUE)
  expect_equal(as.numeric(optimized$objf), as.numeric(reference$objf), tolerance = 1e-8)
  expect_equal(as.matrix(optimized$eta[-1]), as.matrix(reference$eta[-1]), tolerance = 1e-8)
  expect_gt(optimized$env$nAnalyticGradDirect, 0)
  expect_equal(optimized$env$nFDGradFast, 0)
})
