# This acts like the global environment for nlmixr2est
nlmixr2global <- new.env(parent = emptyenv())

.nlmixr2globalReset <- function() {
  nlmixr2global$finalUiCompressed <- TRUE # is the final UI compressed?
  nlmixr2global$rxPredLlik <- TRUE # is this a log-likelihood?
  nlmixr2global$nlmeFitDataAll <- NULL # data for nlme fit
}
