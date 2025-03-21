if (FALSE) {

  nmTest({

    skip_if_not(file.exists(test_path("test-saem-theo_sd.qs")))

    mod <- function() {
      ini({
        tka <- 0.45 ; label("Log Ka")
        tcl <- 1 ; label("Log Cl")
        tv <- 3.45 ; label("Log V")
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        add.sd <- 0.7
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        ipre <- center / v
        ipre ~ add(add.sd)
      })
    }

    mod2 <- function() {
      ini({
        tka <- 0.45 ; label("Log Ka")
        tcl <- 1 ; label("Log Cl")
        tv <- 3.45 ; label("Log V")
        eta.ka ~ 0.6
        eta.cl ~ 0.3
        eta.v ~ 0.1
        lnorm.sd <- 0.1
      })
      model({
        ka <- exp(tka + eta.ka)
        cl <- exp(tcl + eta.cl)
        v <- exp(tv + eta.v)
        d / dt(depot) <- -ka * depot
        d / dt(center) <- ka * depot - cl / v * center
        ipre <- log(center / v)
        ipre ~ add(lnorm.sd)
      })
    }

    dat <- theo_sd
    dat <- dat[!(dat$TIME == 0 & dat$EVID == 0), ]

    estVal <- function(i, n = tot) {
      .cur <- unlist(ops[i, ])
      .est <- c()
      .mod <- c()
      .doIt <- FALSE
      .addProp <- 0
      .lnorm <- FALSE
      .logitNorm <- FALSE
      .probitNorm <- FALSE
      .trans <- FALSE
      .propT <- FALSE
      .powT <- FALSE
      if (.cur["add"] == "add") {
        .mod <- c(.mod, "add(add.sd)")
        .est <- c(.est, c(add.sd = 0.1))
        .doIt <- TRUE
        .addProp <- 1
      }
      if (.cur["add"] == "lnorm") {
        .mod <- c(.mod, "lnorm(lnorm.sd)")
        .est <- c(.est, c(lnorm.sd = 0.1))
        .doIt <- TRUE
        .addProp <- 1
        .trans <- TRUE
        .lnorm <- TRUE
      }
      if (.cur["add"] == "logitNorm") {
        .mod <- c(.mod, "logitNorm(logit.sd, -0.5, 14)")
        .est <- c(.est, c(logit.sd = 0.1))
        .trans <- TRUE
        .doIt <- TRUE
        .addProp <- 1
        .logitNorm <- TRUE
      }
      if (.cur["add"] == "probitNorm") {
        .mod <- c(.mod, "probitNorm(probit.sd, -0.5, 14)")
        .est <- c(.est, c(probit.sd = 0.1))
        .doIt <- TRUE
        .trans <- TRUE
        .addProp <- 1
        .probitNorm <- TRUE
      }
      if (.cur["prop"] == "prop") {
        .mod <- c(.mod, "prop(prop.sd)")
        .est <- c(.est, c(prop.sd = 0.1))
        .doIt <- TRUE
        if (.addProp == 1) .addProp <- 2
      }
      if (.cur["prop"] == "propT") {
        .mod <- c(.mod, "propT(prop.sd)")
        .est <- c(.est, c(prop.sd = 0.1))
        .doIt <- TRUE
        .propT <- TRUE
        if (.addProp == 1) .addProp <- 2
      }
      if (.cur["prop"] == "pow") {
        .mod <- c(.mod, "pow(pow.sd, pw)")
        .est <- c(.est, c(pow.sd = 0.1, pw = 1))
        .doIt <- TRUE
        if (.addProp == 1) .addProp <- 2
      }
      if (.cur["prop"] == "powT") {
        .mod <- c(.mod, "powT(pow.sd, pw)")
        .est <- c(.est, c(pow.sd = 0.1, pw = 1))
        .doIt <- TRUE
        .powT <- TRUE
        if (.addProp == 1) .addProp <- 2
      }
      if (.cur["tbs"] != "") {
        .mod <- c(.mod, paste0(.cur["tbs"], "(lambda)"))
        .est <- c(.est, c(lambda = 1))
        if (.lnorm) {
          .doIt <- FALSE
        } else if ((.logitNorm | .probitNorm) & .cur["tbs"] == "boxCox") {
          .doIt <- FALSE
        }
        .trans <- TRUE
      }
      if (.addProp <= 1 & .cur["addProp"] == "combined1") {
        .doIt <- FALSE
      }
      if (.doIt && !(.trans) && (.propT || .powT)) {
        .doIt <- FALSE
      }
      if (.doIt) {
        ctl1 <- saemControl(print=0, nEm = n, nBurn = n, logLik = TRUE, addProp = .cur["addProp"])
        mod2 <- eval(parse(text = paste0(
          "mod %>% model(ipre~", paste(.mod, collapse = "+"), ") %>% ",
          gsub("c[(]", "ini(", deparse1(.est))
        )))
        v <- .nlmixr(mod2, dat, est = "saem", control = ctl1)
        assign("mod2", mod2, globalenv())
        if (!inherits(v, "nlmixr2FitCore")) {
          message(sprintf("bad model at %d", i))
          print(mod2)
        }
        # saveRDS(v, paste0("test-saem-theo_sd-", i, "-", n, ".rds"))
        .sum <- c(objf = v$objective, v$theta)
        return(invisible(setNames(sapply(.nm, function(x) {
          .sum[x]
        }), .nm)))
      }
      sapply(.nm, function(x) {
        NA_real_
      })
    }

    ctl1 <- saemControl(nEm = 5, nBurn = 5, logLik = TRUE, print = 0, addProp = "combined1")
    ctl2 <- saemControl(nEm = 5, nBurn = 5, logLik = TRUE, print = 0, addProp = "combined2")

    ## tot <- 250
    tot <- 15

    ops <- expand.grid(
      add = c("", "add", "lnorm", "logitNorm", "probitNorm"), prop = c("", "prop", "pow", "powT", "propT"),
      tbs = c("", "yeoJohnson", "boxCox"), addProp = c("combined1", "combined2"),
      stringsAsFactors = FALSE
    )
    ops$id <- seq_along(ops$add)

    .nm <- c("objf", "tka", "tcl", "tv", "lnorm.sd", "logit.sd", "probit.sd", "add.sd", "pow.sd", "pw", "lambda")

    ## estVal(117)

    v <- suppressMessages(suppressWarnings(lapply(seq_along(ops$add), estVal)))

    m <- t(matrix(unlist(v), length(.nm), length(v)))
    dimnames(m) <- list(NULL, .nm)
    val <- cbind(ops, m)
    val <- val[!is.na(val$objf), ]

    for (.n in .nm) {
      val[, .n] <- round(val[[.n]], 2)
    }

    ##qs::qsave(val, file=test_path("test-saem-theo_sd.qs"))

    .test <- qs::qread(test_path("test-saem-theo_sd.qs"))

    for (i in seq_along(.test$add)) {
      test_that(
        with(
          as.list(.test[3,]),
          paste0(
            "add: ", add,
            " prop: ", prop,
            " tbs: ", tbs,
            " addProp: ", addProp
          )
        ), {
          expect_equal(as.list(val[i, ]),
                       as.list(.test[i, ]),
                       tolerance=1e-3)
        })
    }

  })
}
