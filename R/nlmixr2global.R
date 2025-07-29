# This acts like the global environment for nlmixr2est
nlmixr2global <- new.env(parent = emptyenv())

.nlmixr2globalReset <- function(ini=FALSE) {
  nlmixr2global$finalUiCompressed <- TRUE # is the final UI compressed?
  nlmixr2global$rxPredLlik <- FALSE # is this a log-likelihood?
  nlmixr2global$nlmeFitDataAll <- NULL # data for nlme fit
  nlmixr2global$nlmeFitRxModel <- NULL # rx model for nlme fit
  nlmixr2global$nlmeFitRxControl <- NULL # rx control for nlme fit
  nlmixr2global$nlmixr2SimInfo <- NULL # sim info for nlmixr2

  ## Timer information
  nlmixr2global$nlmixr2Time <- NULL # timer for nlmixr2 steps
  nlmixr2global$currentTimingEnvironment <- NULL # current timing environment
  nlmixr2global$timingStackNlmixr <- NULL
  nlmixr2global$timingStack <- NULL
  nlmixr2global$extraTimingTable <- NULL

  nlmixr2global$toRxParam <- "" # string of the params() + cmt() for rx model in focei
  nlmixr2global$toRxDvidCmt <- "" # string of dvid() spec in rxode2 focei model

  nlmixr2global$nlmixr2objectName <- NULL # Allows external methods
                                          # (like those in nlmixr2) to
                                          # assign object name

  nlmixr2global$lastPredSimulationInfo <- NULL # to get observation dataset with pred attached for pred_corr

  nlmixr2global$nlmixrEvalEnv <- new.env(parent=emptyenv()) # evaluate environment for udf

  nlmixr2global$nlmEnv <- new.env(parent=emptyenv()) # nlmEnv data etc for nlm related methods

  nlmixr2global$nlmixr2EstEnv <- new.env(parent=emptyenv())
  nlmixr2global$nlmixr2EstEnv$uiUnfix <- NULL
  nlmixr2global$nlmixr2EstEnv$nlmixrPureInputUi <- NULL

  nlmixr2global$nlsEnv <- new.env(parent=emptyenv())

  nlmixr2global$etaMat <- NULL # eta matrix for foceiControl

  if (ini) {
    nlmixr2global$nlmixr2pipeData <- NULL
    nlmixr2global$nlmixr2pipeControl <- NULL
    nlmixr2global$nlmixr2pipeTable <- NULL
    nlmixr2global$nlmixr2pipeEst <- NULL
  }
}
