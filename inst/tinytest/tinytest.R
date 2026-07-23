# --------------------------------------------------------------------
# Triggering tinytests
# --------------------------------------------------------------------


x <- getOption("distributions3_runtests", "0")
if (isTRUE(as.logical(x))) {
    if (requireNamespace("tinytest", quietly = TRUE)) {
      tinytest::test_package("topmodels")
    }
} else {
    message("[distributions3] Not running tinytests")
}

