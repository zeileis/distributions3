

rm(list = objects())
library("yaml")
source("autotest-functions.R")

###
config <- at_load_config()
config <- lapply(config, at_check_config)

#for (i in seq_along(config)) {
for (i in 1L) {
    ## Name of the output file
    OUTFILE <- file.path("testthat", paste0("test-auto-", names(config)[i]))
    cat("--- generating automatic test set \"", OUTFILE, "\"\n")

    x <- config[[i]]

    ## List to store all the tests required which will be written to the
    ## test file defined by OUTFILE above.
    res <- list()
    res <- append(res, at_setup(x))
    source("autotest-functions.R")

    writeLines(unlist(res), con = OUTFILE)
}


library("testthat")
test_file(OUTFILE)
