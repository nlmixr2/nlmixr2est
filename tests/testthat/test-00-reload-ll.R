test_that("between session saem ll works", {

  src <-
    deparse(quote({
      library(nlmixr2est) #nolint
      rxode2::rxClean()
      one.cmt <- function() {
        ini({
          tka <- 0.45
          tcl <- log(c(0, 2.7, 100))
          tv <- 3.45
          eta.ka ~ 0.6
          eta.cl ~ 0.3
          eta.v ~ 0.1
          add.sd <- 0.7
        })
        model({
          ka <- exp(tka + eta.ka)
          cl <- exp(tcl + eta.cl)
          v <- exp(tv + eta.v)
          linCmt() ~ add(add.sd)
        })
      }
      fit <- suppressMessages(nlmixr(one.cmt, theo_sd, est="saem", control = saemControl(print = 0, nBurn = 1, nEm = 1, nmc = 1, nu = c(1, 1, 1))))
      saveRDS(fit, "fit.rds")
    }))

  src <- src[seq(2, length(src) - 1)]

  rds <- withr::with_tempdir({
    writeLines(src, "000.R")
    .cmd <- file.path(R.home("bin"), "R")
    .args <- c("CMD", "BATCH", "000.R")
    .out <- sys::exec_internal(cmd = .cmd, args = .args, error = FALSE)
    message(paste(readLines(paste0("000.Rout")), collapse="\n"))
    readRDS("fit.rds")
  })

  expect_true(is.numeric(rds$objf))
})
