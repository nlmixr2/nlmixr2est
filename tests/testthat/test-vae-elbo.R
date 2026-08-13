## Phase 4a: full-pipeline ELBO wiring -- encoder (C++) forward -> rxode2 decoder
## -> p(x|z)+KL -> encoder (C++) backward. Verifies the analytic gradient of the
## ELBO w.r.t. every encoder parameter block against finite differences, which
## exercises the entire encoder<->decoder composition end to end.

nmTest({
  test_that("vae ELBO gradient w.r.t. encoder params matches finite differences", {
    theo <- function() {
      ini({
        lka <- log(1.8); lke <- log(0.086); lV <- log(32)
        eta.ka ~ 0.3; eta.ke ~ 0.03; eta.V ~ 0.03
        add.err <- 0.7
      })
      model({
        ka <- exp(lka + eta.ka); ke <- exp(lke + eta.ke); V <- exp(lV + eta.V)
        d/dt(depot) = -ka * depot
        d/dt(central) = ka * depot - ke * central
        cp <- central / V
        cp ~ add(add.err)
      })
    }
    ui <- rxode2::assertRxUi(theo)
    prep <- .vaeDataPrep(ui, nlmixr2data::theo_sd)
    expect_equal(prep$N, 12L)
    expect_equal(prep$neta, 3L)
    am <- .vaeDecoderModel(ui)

    ## nCov MUST come from the prep -- the encoder head is [hDim + nCov] wide and
    ## the encoder is conditioned on the covariates (as in the reference)
    zDim <- prep$zDim; hDim <- 6L; nCov <- ncol(prep$covIn)
    .testSeed(7)
    params <- .vaeEncoderInitParams(zDim, hDim, nCov, prep$zPop, rep(0.1, zDim))
    .testSeed(123); eps <- matrix(rnorm(prep$N * zDim), prep$N, zDim)
    alphaKL <- 0.7

    st <- .vaeElboStep(params, prep, am, prep$zPop, prep$omega, prep$a, alphaKL, eps)
    expect_false(is.null(st$grads))
    expect_true(is.finite(st$loss))

    elboOnly <- function(p) .vaeElboStep(p, prep, am, prep$zPop, prep$omega, prep$a,
                                         alphaKL, eps, withGrad = FALSE)$loss
    h <- 1e-5
    fdMax <- function(name, ncoord = 6) {
      tmpl <- params[[name]]; g <- as.numeric(st$grads[[name]]); flat <- as.numeric(tmpl)
      idx <- if (length(flat) > ncoord) sort(sample(length(flat), ncoord)) else seq_along(flat)
      d <- 0
      for (k in idx) {
        up <- flat; up[k] <- up[k] + h; dn <- flat; dn[k] <- dn[k] - h
        pu <- params; pu[[name]] <- if (is.null(dim(tmpl))) up else array(up, dim(tmpl))
        pd <- params; pd[[name]] <- if (is.null(dim(tmpl))) dn else array(dn, dim(tmpl))
        d <- max(d, abs(g[k] - (elboOnly(pu) - elboOnly(pd)) / (2 * h)))
      }
      d
    }
    .testSeed(1)
    for (nm in c("fcB", "fcW", "Wih", "Whh", "bih")) expect_lt(fdMax(nm), 1e-4)
  })
})
