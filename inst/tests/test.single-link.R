context("single link")

XR <- 30.0
YR <- 0.0
ZR <- 1.8

XL1 <- 0.0
YL1 <- -5000.0
XL2 <- 0.0
YL2 <- 5000.0
WL <- 30.0
HL <- 0.0
NTYP <- 1
VPHL <- 7500.0
EFL <- 30.0

UM <- 1.0
BRGM <- 270.0
CLASM <- 6
MIXHM <- 1000.0

ATIM <- 60.0
Z0 <- 10.0
VS <- 0.0
VD <- 0.0

test_that('CALINE3_LINK_CONTRIBUTIONS', {

	C_ugm3 <- CALINE3_LINK_CONTRIBUTIONS(
		XR, YR, ZR, XL1, YL1, XL2, YL2, WL, HL, NTYP, VPHL, EFL,
		UM, BRGM, CLASM, MIXHM, ATIM, Z0, VS, VD
	)

	C_ppm <- C_ugm3 * 0.0245 / 28.0
	expect_equal(
		round(C_ppm, digits=1),
		array(4.6, c(1, 1))
	)

})

test_that('CALINE3_RECEPTOR_TOTALS', {

  C_ugm3 <- CALINE3_RECEPTOR_TOTALS(
    XR, YR, ZR, XL1, YL1, XL2, YL2, WL, HL, NTYP, VPHL, EFL,
    UM, BRGM, CLASM, MIXHM,
    ATIM, Z0, VS, VD
  )

  C_ppm <- C_ugm3 * 0.0245 / 28.0
  expect_equal(
    round(C_ppm, digits=1),
    array(4.6, c(1, 1))
  )

})
