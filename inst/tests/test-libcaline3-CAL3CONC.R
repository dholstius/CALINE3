context("CAL3CONC")

args <- list(
  U = as.single(1.0),
  W = as.single(30.0),
  W2 = as.single(15.0),
  H = as.single(0.0),
  D = as.double(0.0),
  Z = as.single(1.8),
  UWL = as.single(5000.),
  DWL = as.single(-5000.),
  MIXH = as.single(1000.),
  BASE = as.single(4.0),
  PHI = as.single(1.57),
  Q1 = as.single(38835.),
  PY1 = as.single(0.132025152),
  PY2 = as.single(0.883044183),
  PZ1 = as.single(1.000000000),
  PZ2 = as.single(0.000000000),
  DSTR = as.single(1.0),
  CONC = as.single(0.0))

expect_equal(
  do.call(.Fortran, c("CAL3CONC", args))$CONC,
  as.single(3065.9436),
  tolerance = 1e-6)
