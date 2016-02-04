context("CAL3INTG")

Y <- 1.
SGY <- 1.
LIM <- abs(Y/SGY)
stopifnot(LIM <= 5)
ARG0 <- LIM^2 / (-2.0)
T0 <- 1.0 / (1.0 + 0.23164 * LIM)

CAL3INTG_ <- function (Y, SGY) {
  stopifnot(all(lengths(list(Y, SGY)) == 1))
  call_ <- .Fortran("CAL3INTG", Y = as.single(Y), SGY = as.single(SGY), INTG = as.double(0.0))
  return(as.numeric(call_$INTG))
}

expect_equal(
  CAL3INTG_(Y = 1., SGY = 1.),
  0.3989 * exp(ARG0) * (0.3194 * T0 - 0.3566 * T0 ^ 2 + 1.7815 * T0 ^ 3 - 1.8213 * T0 ^ 4 + 1.3303 * T0 ^ 5),
  tolerance = 0.000001)
