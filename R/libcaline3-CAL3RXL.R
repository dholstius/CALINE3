#' CAL3RXL
#'
#' CAL3RXL predicts the incremental contribution from each link
#' to the concentration at each receptor under a single set of meteorological conditions.
#'
#' @return CAL3RXL returns a matrix of concentrations of size NR x NL,
#'         where NR is the number of receptors and NL is the number of links.
#'
#' @rdname CALINE3
#' @export
#'
#' @examples
#' CAL3RXL(
#'   XR = 30., YR = 0., ZR = 1.8,
#'   XL1 = 0., YL1 = -5000., XL2 = 0., YL2 = 5000.,
#'   WL = 30., HL = 0., NTYP = 1, VPHL = 7500., EFL = 30.,
#'   UM = 1.0, BRGM = 270., CLASM = 6, MIXHM = 1000.,
#'   ATIM = 60., Z0 = 10., VS = 0., VD = 0.)
#'
CAL3RXL <- function(
	XR, YR, ZR,
	XL1, YL1, XL2, YL2, WL, HL, NTYP, VPHL, EFL,
	UM, BRGM, CLASM, MIXHM,
	ATIM, Z0, VS, VD,
	LXR,
	.coerce = TRUE
) {

	if (.coerce) {
		XR    <- as.single(XR)
		YR    <- as.single(YR)
		ZR    <- as.single(ZR)
		XL1   <- as.single(XL1)
		YL1   <- as.single(YL1)
		XL2   <- as.single(XL2)
		YL2   <- as.single(YL2)
		WL    <- as.single(WL)
		HL    <- as.single(HL)
		NTYP  <- as.integer(NTYP)
		VPHL  <- as.single(VPHL)
		EFL   <- as.single(EFL)
		UM    <- as.single(UM)
		BRGM  <- as.single(BRGM)
		CLASM <- as.integer(CLASM)
		MIXHM <- as.single(MIXHM)
		ATIM  <- as.single(ATIM)
		Z0    <- as.single(Z0)
		VS    <- as.single(VS)
		VD    <- as.single(VD)
	}

	NR <- length(XR)
	NL <- length(XL1)

	# Link-receptor pairs
	if (missing(LXR)) {
	  message("Processing all ", NR * NL, " possible link-receptor pairs")
	  LXR <- array(TRUE, dim = c(NR, NL))
	} else {
	  message("Processing ", sum(LXR), " of ", NR * NL, " possible link-receptor pairs")
	  stopifnot(dim(LXR) == c(NR, NL))
	  stopifnot(is.logical(LXR))
	}

	# Call native code, using array C for results
	shape <- c(NR, NL)
	C <- as.single(array(0.0, dim = shape))
	retval <- .Fortran(
	  'CAL3RXL',
	  NR, XR, YR, ZR,
	  NL, XL1, YL1, XL2, YL2, WL, HL, NTYP, VPHL, EFL,
	  UM, BRGM, CLASM, MIXHM,
	  ATIM, Z0, VS, VD,
	  LXR = LXR,
	  C = C,
	  PACKAGE = "CALINE3"
	)

  # Returned values are in g/m^3
	with(retval, array(as.double(C), dim = shape))

}

