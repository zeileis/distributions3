

rm(list = objects())
library("yaml")
source("autotest-functions.R")

## TODO: Remove, just for testing
library("testthat")

###
config <- autotest_load_config()
#str(config)
config <- lapply(config, autotest_check_config)

#for (i in seq_along(config)) {
for (i in 1L) {
    ## Name of the output file
    OUTFILE <- file.path("testthat", paste0("test-auto-", names(config)[i]))
    OUTFILE <- gsub("\\.yaml$", ".R", OUTFILE)
    cat("--- generating automatic test set \"", OUTFILE, "\"\n")

    x <- config[[i]]

    ## List to store all the tests required which will be written to the
    ## test file defined by OUTFILE above.
    res <- list()
    res <- append(res, autotest_setup(x))
    res <- append(res, autotest_constructor(x))
    #res <- append(res, autotest_basic_methods(x))
    #res <- append(res, autotest_pdf(x))
    #res <- append(res, autotest_log_pdf(x))
    #res <- append(res, autotest_constructor_return(x))



    autotest_setup(x)

    writeLines(unlist(res), con = OUTFILE)

    ## TODO: Remove, just for testing
    test_file(OUTFILE)
}



