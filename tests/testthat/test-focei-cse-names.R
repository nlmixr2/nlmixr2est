test_that('repeated outer CSE passes preserve existing temporaries', {
  x <- 'rx_expr_0~a+b\ny=rx_expr_0*(c+d)+(c+d)'
  evaluate <- function(txt) {
    env <- list2env(list(a=2, b=3, c=7, d=11))
    eval(parse(text=gsub('~', '<-', txt, fixed=TRUE)), env)
    env$y
  }
  for (pass in 1:3) {
    out <- rxode2:::.rxOptExprC(.foceiCseRename(x, pass))
    if (is.na(out)) {
      expect_gt(pass, 1L)
      break
    }
    x <- out
    expect_equal(evaluate(x), 108)
  }
})

test_that('outer CSE preserves a small ODE and its sensitivity', {
  original <- paste(
    'rx_expr_0~a+b', 'k=rx_expr_0*(c+d)+(c+d)',
    'd/dt(A)=-k*A', 'd/dt(S)=-k*S-A',
    'A(0)=10', 'S(0)=0', 'pred=A', sep='\n'
  )
  optimized <- original
  for (pass in 1:3) {
    out <- rxode2:::.rxOptExprC(.foceiCseRename(optimized, pass))
    if (is.na(out)) break
    optimized <- out
  }
  params <- c(a=.02, b=.03, c=.07, d=.11)
  events <- data.frame(time=0:24, evid=0, amt=0)
  solve <- function(txt) as.data.frame(rxode2::rxSolve(
    rxode2::rxode2(txt), params=params, events=events
  ))
  ref <- solve(original)
  actual <- solve(optimized)
  expect_true(all(is.finite(as.matrix(actual))))
  expect_equal(actual[c('A', 'S', 'pred')], ref[c('A', 'S', 'pred')], tolerance=1e-8)
})

