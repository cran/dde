## ----echo = FALSE, results = "hide"-------------------------------------------
output_lang <- function(x, lang) {
  writeLines(c(sprintf("```%s", lang), x, "```"))
}
output_r <- function(x) output_lang(x, "r")
output_c <- function(x) output_lang(x, "c")
knitr::opts_chunk$set(
  error = FALSE,
  fig.width = 7,
  fig.height = 5)
set.seed(1)

## -----------------------------------------------------------------------------
lorenz_dde <- function(t, y, p) {
  sigma <- p$sigma
  R     <- p$R
  b     <- p$b
  y1 <- y[[1L]]
  y2 <- y[[2L]]
  y3 <- y[[3L]]
  c(sigma * (y2 - y1),
    R * y1 - y2 - y1 * y3,
    -b * y3 + y1 * y2)
}

## -----------------------------------------------------------------------------
p <- list(sigma = 10.0,
          R     = 28.0,
          b     =  8.0 / 3.0)

tt <- seq(0, 100, length.out = 50001)
y0 <- c(1, 1, 1)
yy <- dde::dopri(y0, tt, lorenz_dde, p)

## -----------------------------------------------------------------------------
par(mar=rep(.5, 4))
plot(yy[, c(2, 4)], type = "l", lwd = 0.5, col = "#00000066",
     axes = FALSE, xlab = "", ylab = "")

## ----eval = FALSE-------------------------------------------------------------
#  lorenz_ds <- function(t, y, p) {
#    sigma <- 10.0
#    R     <- 28.0
#    b     <-  8.0 / 3.0
#    y1 <- y[[1L]]
#    y2 <- y[[2L]]
#    y3 <- y[[3L]]
#    list(c(sigma * (y2 - y1),
#           R * y1 - y2 - y1 * y3,
#           -b * y3 + y1 * y2))
#  }
#  yy_ds <- deSolve::ode(y0, tt, lorenz_ds, p)

## -----------------------------------------------------------------------------
yy <- dde::dopri(y0, tt, lorenz_dde, p, return_statistics = TRUE)
attr(yy, "statistics")

## -----------------------------------------------------------------------------
yy2 <- dde::dopri(y0, range(tt), lorenz_dde, p, return_minimal = TRUE,
                  n_history = 5000, return_history = TRUE)

## -----------------------------------------------------------------------------
dim(yy2)
h <- attr(yy2, "history")
dim(h)

## -----------------------------------------------------------------------------
yy2 <- dde::dopri_interpolate(h, tt)
all.equal(yy2, yy[, -1], check.attributes = FALSE)

## -----------------------------------------------------------------------------
seir <- function(t, y, p) {
  b <- 0.1
  N <- 1e7
  beta <- 10.0
  sigma <- 1.0 / 3.0
  delta <- 1.0 / 21.0
  t_latent <- 14.0
  Births <- N * b
  surv <- exp(-b * t_latent)

  S <- y[[1L]]
  E <- y[[2L]]
  I <- y[[3L]]
  R <- y[[4L]]

  tau <- t - t_latent
  y_lag <- dde::ylag(tau, c(1L, 3L)) # Here is ylag!
  S_lag <- y_lag[[1L]]
  I_lag <- y_lag[[2L]]

  new_inf <- beta * S * I / N
  lag_inf <- beta * S_lag * I_lag * surv / N

  c(Births  - b * S - new_inf + delta * R,
    new_inf - lag_inf - b * E,
    lag_inf - (b + sigma) * I,
    sigma * I - b * R - delta * R)
}

## -----------------------------------------------------------------------------
y0 <-  y0 <- c(1e7 - 1, 0, 1, 0)
tt <- seq(0, 365, length.out = 100)
yy <- dde::dopri(y0, tt, seir, NULL, n_history = 1000L, return_history = FALSE)
matplot(tt, yy[, 2:5], type="l")

## ----echo=FALSE, results="asis"-----------------------------------------------
output_r(capture.output(print(seir)))

## -----------------------------------------------------------------------------
seir_deSolve <- function(t, y, parms) {
  b <- 0.1
  N <- 1e7
  beta <- 10
  sigma <- 1 / 3
  delta <- 1 / 21
  t_latent <- 14.0
  I0 <- 1
  Births <- N * b
  surv <- exp(-b * t_latent)

  S <- y[[1L]]
  E <- y[[2L]]
  I <- y[[3L]]
  R <- y[[4L]]

  tau <- t - t_latent
  if (tau < 0.0) { # NOTE: assuming that t0 is always zero
    S_lag <- parms$S0
    I_lag <- parms$I0
  } else {
    y_lag <- deSolve::lagvalue(tau, c(1L, 3L))
    S_lag <- y_lag[[1L]]
    I_lag <- y_lag[[2L]]
  }

  new_inf <- beta * S * I / N
  lag_inf <- beta * S_lag * I_lag * surv / N

  list(c(Births  - b * S - new_inf + delta * R,
         new_inf - lag_inf - b * E,
         lag_inf - (b + sigma) * I,
         sigma * I - b * R - delta * R))
}

## -----------------------------------------------------------------------------
y0 <-  y0 <- c(1e7 - 1, 0, 1, 0)
tt <- seq(0, 365, length.out = 100)
initial <- list(S0 = y0[[1]], I0 = y0[[3]])
yy_ds <- deSolve::dede(y0, tt, seir_deSolve, initial)

## -----------------------------------------------------------------------------
yy_dde <- dde::dopri(y0, tt, seir, NULL, n_history = 1000L,
                     return_history = FALSE)
op <- par(mfrow=c(1, 2), mar=c(4, .5, 1.4, .5), oma=c(0, 2, 0, 0))
matplot(tt, yy_dde[, -1], type="l", main = "dde")
matplot(tt, yy_ds[, -1], type="l", main = "deSolve", yaxt="n")

## -----------------------------------------------------------------------------
tR <- microbenchmark::microbenchmark(times = 30,
  deSolve = deSolve::dede(y0, tt, seir_deSolve, initial),
  dde = dde::dopri(y0, tt, seir, NULL, n_history = 1000L,
                   return_history = FALSE))
tR

## ----include = FALSE----------------------------------------------------------
.dlls <- local({
  build <- sprintf("%s.c", c("seir", "seir_ds"))
  files <- file.path(dde:::dde_example_path(), build)
  lapply(files, dde:::shlib, "dde_")
})

## ----echo = FALSE, results = "asis"-------------------------------------------
output_c(readLines(system.file("examples/seir_ds.c", package = "dde")))

## -----------------------------------------------------------------------------
initial <- c(0.0, y0[[1]], y0[[3]])

zz_ds <- deSolve::dede(y0, tt, "seir_deSolve", initial,
                       initfunc = "seir_initmod", dllname = "dde_seir_ds")
zz_dde <- dde::dopri(y0, tt, "seir", numeric(), dllname = "dde_seir",
                     n_history = 1000L, return_history = FALSE)

## -----------------------------------------------------------------------------
all.equal(zz_ds, yy_ds, check.attributes = FALSE)
all.equal(zz_dde, yy_dde, check.attributes = FALSE)

## -----------------------------------------------------------------------------
tC <- microbenchmark::microbenchmark(
  deSolve = deSolve::dede(y0, tt, "seir_deSolve", initial,
                          initfunc = "seir_initmod", dllname = "dde_seir_ds"),
  dde = dde::dopri(y0, tt, "seir", numeric(), dllname = "dde_seir",
                   n_history = 1000L, return_history = FALSE))
tC

## -----------------------------------------------------------------------------
ptr <- getNativeSymbolInfo("seir")
tC2 <- microbenchmark::microbenchmark(
  deSolve = deSolve::dede(y0, tt, "seir_deSolve", initial,
                          initfunc = "seir_initmod", dllname = "dde_seir_ds"),
  dde = dde::dopri(y0, tt, "seir", numeric(), dllname = "dde_seir",
                   n_history = 1000L, return_history = FALSE),
  dde2 = dde::dopri(y0, tt, ptr, numeric(), n_history = 1000L,
                    return_history = FALSE, return_minimal = TRUE))
tC2

## ----include = FALSE----------------------------------------------------------
for (x in .dlls) {
  tryCatch({
    dyn.unload(x$dll)
    unlink(x$tmp, recursive = TRUE)
  }, error = function(e) NULL)
}

