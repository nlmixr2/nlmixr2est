# This acts like the global environment for nlmixr2est
nlmixr2global <- new.env(parent = emptyenv())

.nlmixr2globalReset <- function() {
  nlmixr2global$finalUiCompressed <- TRUE
  nlmixr2global$rxPredLlik <- TRUE
}
