test_that("saem mu reference 1", {
  skip_if_not(utils::packageVersion("rxode2") >= "4.0.0")
  theo_sd2 <- nlmixr2data::theo_sd

  theo_sd2$lwt<-log(theo_sd2$WT/70)

  # The basic model consists of an ini block that has initial estimates
  # Original mu-referencing
  one.compartment <- function() {
    ini({
      tka <- log(1.57); label("Ka")
      tcl <- log(2.72); label("Cl")
      tv <- log(31.5); label("V")
      covwt<- 0.01
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + covwt*lwt)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  fit1 <- .nlmixr(one.compartment, theo_sd2, est="saem",
                  control=saemControl(print = 0,seed = 1234, nBurn = 1, nEm = 1,
                                      calcTables = FALSE))

  # true mu expression should not have information in $runInfo
  expect_true(is.null(fit1$runInfo))

  # mu2-referencing
  theo_sd2 <- nlmixr2data::theo_sd
  one.compartment <- function() {
    ini({
      tka <- log(1.57); label("Ka")
      tcl <- log(2.72); label("Cl")
      tv <- log(31.5); label("V")
      covwt <- 0.01
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      cl <- exp(tcl + eta.cl + log(WT/70)*covwt)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }


  fit2 <-
    withr::with_options(list(cli.unicode=FALSE),{
      .nlmixr(one.compartment, theo_sd2, est="saem",
              control=saemControl(print = 0,seed = 1234, nBurn = 1, nEm = 1,
                                  calcTables = FALSE))
    })

  expect_true(grepl("mu2 item:", fit2$runInfo[1]))

  # mu3 referencing
  one.compartment <- function() {
    ini({
      tka <- log(1.57); label("Ka")
      tcl <- log(2.72); label("Cl")
      tv <- log(31.5); label("V")
      covwt <- 0.01
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      wt70 <- WT/70
      cl <- exp(tcl + eta.cl + log(wt70)*covwt)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  fit3 <-   withr::with_options(list(cli.unicode=FALSE), {
    .nlmixr(one.compartment, theo_sd2, est="saem",
            control=saemControl(print = 0, seed = 1234, nBurn = 1, nEm = 1,
                                calcTables = FALSE))
  })

  expect_true(grepl("mu3 item", fit3$runInfo[1]))


  # mu4 referencing
  one.compartment <- function() {
    ini({
      tka <- log(1.57); label("Ka")
      tcl <- log(2.72); label("Cl")
      tv <- log(31.5); label("V")
      covwt <- 0.01
      eta.ka ~ 0.6
      eta.cl ~ 0.3
      eta.v ~ 0.1
      add.sd <- 0.7
    })
    model({
      ka <- exp(tka + eta.ka)
      wt70 <- WT/70
      if (wt70 > 1) {
        wt70s <- "larger than 70"
      } else {
        wt70s <- "smaller than 70"
      }
      cl <- exp(tcl + eta.cl + (wt70s == "larger than 70")*covwt)
      v <- exp(tv + eta.v)
      d/dt(depot) <- -ka * depot
      d/dt(center) <- ka * depot - cl / v * center
      cp <- center / v
      cp ~ add(add.sd)
    })
  }

  fit4 <- withr::with_options(list(cli.unicode=FALSE), {
    .nlmixr(one.compartment, theo_sd2, est="saem",
            control=saemControl(print=0,seed = 1234, nBurn = 1, nEm = 1,
                                calcTables = FALSE))
  })

  expect_true(grepl("mu4 item", fit4$runInfo[1]))
})
