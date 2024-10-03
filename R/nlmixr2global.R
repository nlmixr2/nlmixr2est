# This acts like the global environment for nlmixr2est
nlmixr2global <- new.env(parent = emptyenv())

.nlmixr2globalReset <- function(ini=FALSE) {
  nlmixr2global$finalUiCompressed <- TRUE # is the final UI compressed?
  nlmixr2global$rxPredLlik <- TRUE # is this a log-likelihood?
  nlmixr2global$nlmeFitDataAll <- NULL # data for nlme fit
  nlmixr2global$nlmeFitRxModel <- NULL # rx model for nlme fit
  nlmixr2global$nlmeFitRxControl <- NULL # rx control for nlme fit
  nlmixr2global$nlmixr2SimInfo <- NULL # sim info for nlmixr2

  ## Timer information
  nlmixr2global$nlmixr2Time <- NULL # timer for nlmixr2 steps
  nlmixr2global$currentTimingEnvironment <- new.env(parent=emptyenv()) # current timing environment
  nlmixr2global$timingStackNlmixr <- NULL
  if (ini) {
    nlmixr2global$timingStack <- NULL
  }
}
