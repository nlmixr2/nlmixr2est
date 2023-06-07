#' UI modify covariates with reps
#'
#' @param expr expression to change
#' @param old old expression to change
#' @param new new expression to change
#' @return new expression with replacement
#' @author Matthew L. Fidler
#' @noRd
.uiModifyForCovsRep <- function(expr, old, new) {
  if (identical(expr, old)) return(new)
  if (is.call(expr)) {
    as.call(c(expr[[1]],lapply(expr[-1], .uiModifyForCovsRep, old=old, new=new)))
  } else {
    expr
  }
}
#' This function handles mu2 covariates
#'
#' In general the dataset is modified with nlmixrMuDerCov# and the mu2
#' expressions are changed to traditional mu-expressions
#'  
#' @param data input dataset
#' @param ui input ui
#' @return a list with list(ui=mu referenced ui, data=mu referenced dataset)
#' @author Matthew L. Fidler
#' @noRd
.uiModifyForCovs <- function(ui, data) {
  .datEnv <- new.env(parent=emptyenv())
  .datEnv$data <- data
  .datEnv$model <- rxode2::as.model(ui)
  lapply(seq_along(ui$mu2RefCovariateReplaceDataFrame$covariate),
         function(i) {
           .datEnv$data[[paste0("nlmixrMuDerCov", i)]] <-
             with(.datEnv$data,
                  eval(str2lang(ui$mu2RefCovariateReplaceDataFrame$covariate[i])))
           .new <- str2lang(paste0("nlmixrMuDerCov", i, "*",
                                   ui$mu2RefCovariateReplaceDataFrame$covariateParameter[i]))
           .old <- str2lang(ui$mu2RefCovariateReplaceDataFrame$modelExpression[i])
           .datEnv$model <- .uiModifyForCovsRep(.datEnv$model, .old, .new)
           invisible()
         })
  ui2 <- ui
  rxode2::model(ui2) <- .datEnv$model
  ui2 <- rxode2::rxUiDecompress(ui2)
  list(ui=ui2, data=.datEnv$data)
}
