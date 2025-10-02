#' @export
nmObjGet.etaMat <- function(x, ...) {
  .ui <- x[[1]]
  if (is.null(.ui$eta)) return(NULL)
  if (is.null(.ui$iov)) {
    as.matrix(.ui$eta[-1])
  } else {
    .eta <- as.matrix(.ui$eta[-1])
    .n <- names(.ui$iov)
    as.matrix(do.call(`cbind`,
                      c(list(.eta),
                        lapply(.n, function(n) {
                          .dt <- data.table::as.data.table(.ui$iov[[n]])
                          .frm <- eval(str2lang(paste0("ID ~ ", n)))
                          .nr <- names(.dt)[-(1:2)]
                          do.call(`cbind`,
                                  lapply(.nr, function(nr) {
                                    .dt0 <- .dt[,c("ID", n, nr)]
                                    .df <- as.data.frame(data.table::dcast(.dt,
                                                                           formula=.frm,
                                                                           value.var=nr)[, -1])
                                    names(.df) <- paste0("rx.",nr, ".", names(.df))
                                    .df
                                  }))
                        }))))
  }
}
