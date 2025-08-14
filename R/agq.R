.agq <- function(neta=2, nAQD=3) {
  .gh <- fastGHQuad::gaussHermiteData(nAQD)
  .x <- .gh$x
  .w <- .gh$w
  if (nAQD %% 2 == 1) {
    # If nAQD is odd, have the zero weight at the beginning so that
    # the F value is cached and it doesn't need to evaluate it twice
    .zero <- (nAQD+1L)/2L
    .x <- c(0, .x[-.zero])
    .w <- c(.w[.zero], .w[-.zero])
  }
  .x <-   as.matrix(do.call("expand.grid", lapply(1:neta, function(x) .x)))
  .w <-   as.matrix(do.call("expand.grid", lapply(1:neta, function(x) .w)))
  list(
    x = .x,
    w = .w,
    n = nrow(.x),
    neta = neta,
    nAQD = nAQD
  )
}
