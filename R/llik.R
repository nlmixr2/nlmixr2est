#' llik Binomial exported from stan
#'
#' @param x Number of successes observed
#' @param size number of trials (zero or more).
#' @param prob probability of success on each trial.
#' @return list with `fx` and the dbinom distribution and `J` the
#'   Jacobian. This is what is used internally for estimation methods.
#' @author Wenping Wang, Matthew L. Fidler, Stan team
#' @export
#' @examples
#'
#' llikBinomial(46:54, 100, 0.5)
#'
#' n <- 2000
#' k <- seq(0, n, by = 20)
#'
#' llikBinomial(k, n, pi/10)
#'
llikBinomial <- function(x, size, prob) {
  checkmate::assertIntegerish(x, lower=0, any.missing=FALSE)
  checkmate::assertIntegerish(size, lower=0, any.missing=FALSE)
  checkmate::assertNumeric(prob, lower=0, upper=1, any.missing=FALSE)
  .df <- data.frame(x=as.numeric(x), size=as.numeric(size), prob=as.numeric(prob))
  if (any(.df$x > .df$size)) {
    stop("some binomial observations 'x' are more than the number of experiments 'size'",
         call.=FALSE)
  }
  .Call(`_nlmixr2est_llikBinomialInternal`, .df$x, .df$size, .df$prob)
}
#' llik Normal exported from stan
#'
#' @param x observed values
#' @param mean vector of means.
#' @param sd vector of standard deviations.
#' @return list with `fx` and the dbinom distribution and `J` the
#'   Jacobian. This is what is used internally for estimation methods.
#' @author Matthew L. Fidler
#' @export
#' @examples
#'
#' llikNormal(c(-0.2, 0, 0.2), mean=0, sd=1)
#'
llikNormal <- function(x, mean=0, sd=1) {
  checkmate::assertNumeric(x, any.missing=FALSE)
  checkmate::assertNumeric(mean, any.missing=FALSE)
  checkmate::assertNumeric(sd, any.missing=FALSE)
  .df <- data.frame(x=as.numeric(x), mean=as.numeric(mean), sd=as.numeric(sd))
  .mat <- as.matrix(.df[-1, ])
  .Call(`_nlmixr2est_llikNormalInternal`, .df$x, .mat)
}
