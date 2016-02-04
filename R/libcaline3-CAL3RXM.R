#' CAL3RXM
#'
#' Given a sample of representative meteorological conditions, CALINE3_RECEPTOR_TOTALS
#' predicts cumulative concentrations at each receptor (from incremental concentrations
#' contributed by each link)
#'
#' All coordinates are in meters unless otherwise specified. By default, predicted
#' concentrations are returned in units of grams per cubic meter (µg/m^3).
#'
#' @param XR x-coordinates of the receptors
#' @param YR y-coordinates of the receptors
#' @param ZR z-coordinates of the receptors (height above ground level, usually 1.8m)
#' @param XL1 starting x-coordinates of the links
#' @param YL1 starting y-coordinates of the links
#' @param XL2 ending x-coordinates of the links
#' @param YL2 ending y-coordinates of the links
#' @param WL widths of the links
#' @param HL heights of the links (above ground level)
#' @param NTYP link classifications (1=at grade, 2=bridge, 3=fill, 4=depressed)
#' @param VPHL link-level traffic volumes, in vehicles per hour
#' @param EFL link-level emission factors, in grams per vehicle-mile per hour
#' @param UM wind speeds, in meters per second (not less than 1.0)
#' @param BRGM wind bearings, in degrees (direction wind is blowing from)
#' @param CLASM stability classes (1, 2, 3, 4, 5, or 6)
#' @param MIXHM mixing heights, in meters (over 1000 skips mixing height calculations)
#' @param ATIM averaging time, in minutes (usually 60)
#' @param Z0 surface roughness, in centimeters
#' @param VS settling velocity, in cm/sec
#' @param VD deposition velocity, in cm/sec
#' @param LXR whether to process (or skip) a given receptor-link pair
#' @param .coerce force arguments to be cast to correct type
#'
#' @return CAL3RXM returns a matrix of concentrations of size NR x NM,
#'         where NR is the number of receptors and NM is the number of meteorological conditions
#'
#' @useDynLib CALINE3
#' @rdname CALINE3
#' @export
#'
#' @examples
#' CAL3RXM(
#'   XR = 30., YR = 0., ZR = 1.8,
#'   XL1 = 0., YL1 = -5000., XL2 = 0., YL2 = 5000.,
#'   WL = 30., HL = 0., NTYP = 1, VPHL = 7500., EFL = 30.,
#'   UM = 1.0, BRGM = 270., CLASM = 6, MIXHM = 1000.,
#'   ATIM = 60., Z0 = 10., VS = 0., VD = 0.)
#'
CAL3RXM <- function(
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

  # Receptor specifications
  NR <- as.integer(length(XR))
  stopifnot(all.equal(NR, length(YR), length(ZR)))
  if(any(is.na(c(XR, YR, ZR))))
    stop("Receptor coordinates cannot include NA values.")

  # Link specifications
  NL <- as.integer(length(XL1))
  stopifnot(all.equal(NL, length(YL1), length(XL2), length(YL2),
                      length(WL), length(HL), length(NTYP), length(VPHL), length(EFL)))
  stopifnot(lapply(list(XR, YR, ZR), is.numeric) == TRUE)
  stopifnot(lapply(list(XL1, YL1, XL2, YL2, WL, HL, NTYP, VPHL, EFL), is.numeric) == TRUE)
  stopifnot(lapply(list(UM, BRGM, CLASM, MIXHM), is.numeric) == TRUE)
  if(any(is.na(c(XL1, YL1, XL2, YL2))))
    stop("Link coordinates cannot include NA values.")
  if(any(is.na(WL)))
    stop("Link widths cannot include NA values.")
  if(any(is.na(HL)))
    stop("Link heights cannot include NA values.")
  if(any(is.na(NTYP)))
    stop("Link classifications cannot include NA values.")
  if(any(is.na(VPHL)))
    stop("Link flows cannot include NA values.")
  if(any(is.na(EFL)))
    stop("Link emission factors cannot include NA values.")
  if (!is.integer(NTYP)) {
    if (is.character(NTYP)) {
      # Can't pass characters to .Fortran() with DUP = FALSE.
      clas.lookup <- list(AG=0, BR=1, FL=2, DP=3, `At Grade`=0, `Bridge`=1, `Fill`=2, `Depressed`=3)
      NTYP <- as.integer(clas.lookup[NTYP])
    } else {
      stop('NTYP argument must be character or integer')
    }
  }

  # Meteorology specifications.
  NM <- as.integer(length(UM))
  stopifnot(all.equal(NM, length(BRGM), length(CLASM), length(MIXHM)))
  if(any(is.na(UM)))
    stop("Wind speeds cannot include NA values. See ?read.ISC for more.")
  if(any(is.na(BRGM)) || any(BRGM < 0) || any(BRGM > 360))
    stop("Wind bearings must be between 0 and 360 degrees, and cannot include NA values. See ?read.ISC for more.")
  if(any(is.na(CLASM)))
    stop("Stability classes must not include NA values. See ?read.ISC for more.")
  if(any(CLASM > 6)) {
    warning("Replacing stability class 7 with class 6.")
    CLASM <- pmin(CLASM, 6)
  }
  if(any(MIXHM < 0))
    stop("Mixing heights must not be negative")
  if(any(is.na(MIXHM)))
    stop("Mixing heights cannot include NA values.")

  # Link-receptor pairs
  if (missing(LXR)) {
    message("Processing all ", NR * NL, " possible link-receptor pairs")
    LXR <- as.logical(array(TRUE, dim = c(NR, NL)))
  } else {
    message("Processing ", sum(LXR), " of ", NR * NL, " possible link-receptor pairs")
    stopifnot(is.logical(LXR))
    stopifnot(identical(dim(LXR), c(NR, NL)))
  }

  # Call native code, using array C for results
  shape <- c(NR, NM)
  C <- as.single(array(0.0, dim = shape))
  retval <- .Fortran(
    'CAL3RXM',
    NR, XR, YR, ZR,
    NL, XL1, YL1, XL2, YL2, WL, HL, NTYP, VPHL, EFL,
    NM,
    UM, BRGM, CLASM, MIXHM,
    ATIM, Z0, VS, VD,
    LXR = LXR,
    C = C,
    PACKAGE = "CALINE3"
  )

  # Returned values are in g/m^3
  with(retval, array(as.double(C), dim = shape))

}
