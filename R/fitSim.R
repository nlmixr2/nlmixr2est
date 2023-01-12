getBaseSimModel.nlmixr2.focei <- function(obj) {
  .base <-   rxode2::rxCombineErrorLines(obj$ui)
  .base[[1]] <- quote(`rxModelVars`)
  .base <- eval(.base)
  # focei models are always pruned
  .base <- rxode2::rxPrune(.base)
  .ctl <- obj$foceiControl
  if (.ctl$sumProd) {
    .malert("stabilizing round off errors in simulation problem...")
    .base <- rxode2::rxSumProdModel(.base)
    .msuccess("done")
  }
  if (.ctl$optExpression) {
    .base <- rxode2::rxOptExpr(.base, ifelse(.getRxPredLlikOption(),
                                             "simulation llik model",
                                             "simulation model"))
  }
  eval(parse(text=paste0("quote(rxode2({\n", .base, "\n}))")))
}

getBaseSimModel.nlmixr2FitCoreSilent <- function(obj) {
  .est <- obj$est
  if (.est %in% c("focei", "foce", "foi", "fo", "posthoc")) {
    if (!exists(".foceiSim", obj)) {
      assign(".foceiSim", getBaseSimModel.nlmixr2.focei(obj), envir=obj)
    }
    return(obj$.foceiSim)
  }
  NextMethod()
}
