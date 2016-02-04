context("rural curved")

NR <- 4
XR <- c(400, 100, 200, 100)
YR <- c(1700, 1500, 1300, 350)
ZR <- rep(1.8, NR)

NL <- 10
XL1 <- c(-707, 0, 120, 150, 150, 175, 265, 350, 475, 650)
YL1 <- c(-707, 0, 175, 350, 1350, 1510, 1640, 1760, 1830, 1850)
XL2 <- c(0, 120, 150, 150, 175, 265, 350, 475, 650, 1650)
YL2 <- c(0, 175, 350, 1350, 1510, 1640, 1760, 1830, 1850, 1850)
NTYP <- rep(1, NL)
HL <- rep(0, NL)
WL <- rep(28, NL)
VPHL <- rep(8500, NL)
EFL <- rep(30, NL)

test_that('array lengths', {
	expect_equal(length(XR), NR)
	expect_equal(length(YR), NR)
	expect_equal(length(ZR), NR)
	expect_equal(length(XL1), NL)
	expect_equal(length(YL1), NL)
	expect_equal(length(XL2), NL)
	expect_equal(length(YL2), NL)
	expect_equal(length(HL), NL)
	expect_equal(length(WL), NL)
	expect_equal(length(VPHL), NL)
	expect_equal(length(EFL), NL)
})

UM <- 1
BRGM <- 45
CLASM <- 6
MIXHM <- 1000

ATIM <- 60.0
Z0 <- 50.0
VS <- 0.0
VD <- 0.0

expected_result <- structure(
  c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.8, 0, 0, 0, 0, 0,
    1.5, 0, 0, 0, 3.7, 0, 0, 0, 2.1, 0, 0, 3.1, 0.4, 0.1, 0, 0, 0, 1.3, 0.5), .Dim = c(4L, 10L))

test_that('CAL3RXL', {
  C_gm3 <- CAL3RXL(
    XR, YR, ZR,
    XL1, YL1, XL2, YL2, WL, HL, NTYP, VPHL, EFL,
    UM, BRGM, CLASM, MIXHM,
    ATIM, Z0, VS, VD
  )
  C_ppm <- C_gm3 * 0.0245 / 28.0
	expect_that(round(C_ppm, digits=1), equals(expected_result))
})

test_that('CAL3RXM', {
  C_gm3 <- CAL3RXM(
    XR, YR, ZR,
    XL1, YL1, XL2, YL2, WL, HL, NTYP, VPHL, EFL,
    UM, BRGM, CLASM, MIXHM,
    ATIM, Z0, VS, VD
  )
  C_ppm <- C_gm3 * 0.0245 / 28.0
  expect_that(abs(mean(as.vector(round(C_ppm, digits=1)) - rowSums(expected_result))), is_less_than(0.1))
})
