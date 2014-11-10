This is a rarely-updated package that breaks out Fortran routines used by another package (Rcaline),
so that users of the more frequently-updated Rcaline need not have access to a Fortran compiler.

It passes devtools::check() on OS X.

devtools::build_win() says:

    * checking for code which exercises the package ... WARNING
    No examples, no tests, no vignettes

However, there are several `testthat` cases in inst/tests/, all of which pass.
