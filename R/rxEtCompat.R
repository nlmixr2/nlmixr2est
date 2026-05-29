# Only define dim.rxEt if rxode2 doesn't already register one — avoids
# overriding a future rxode2-native implementation.
if (!isS3method("dim", "rxEt")) {
  dim.rxEt <- function(x) {
    dim(as.data.frame(x, all = TRUE))
  }
}
