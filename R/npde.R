
.npdeSim <- function(object, nsim = 300, ties = TRUE, seed = 1009, updateObject = TRUE,
                     cholSEtol = (.Machine$double.eps)^(1 / 3), ..., addDosing=FALSE, subsetNonmem=TRUE) {
  set.seed(seed)
  .si <- object$simInfo
  .si$rx <- eval(.getSimModel(object, hideIpred=TRUE))
  .si$object <- object
  .si$returnType <- "data.frame.TBS"
  .si$nsim <- nsim
  .si <- c(.si, list(...))
  .pt <- proc.time()
  .si$dfObs <- 0
  .si$dfSub <- 0
  .si$thetaMat <- NA
  .si$addDosing <- addDosing
  .si$subsetNonmem <- subsetNonmem
  do.call(rxode2::rxSolve, .si)
}
