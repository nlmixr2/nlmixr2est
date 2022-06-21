test_that("between session saem ll works", {

  deparse(quote({
    library(nlmixr2est) #nolint
    rxode2::rxClean()
    one.cmt <- function() {
      ini({
        ## You may label each parameter with a comment
        tka <- 0.45 # Ka
        tcl <- log(c(0, 2.7, 100)) # Log Cl
        ## This works with interactive models
        ## You may also label the preceding line with label("label text")
        tv <- 3.45; label("log V")
        ## the label("Label name") works with all models
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
    fit <- nlmixr(one.cmt, theo_sd, est="saem", control=saemControl(print=0))
    saveRDS(fit, "fit.rds")
  })) -> src

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
