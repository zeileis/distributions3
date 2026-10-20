
#' Reading/Parsing YAML File
#'
#' @param dir character. Name of the directory containing the `.yaml` files.
#'
#' @reutrn Named list with the content of the YAML files. The names of the
#' list elements corresponds to the YAML file name.
at_load_config <- function(dir = "autotest-config") {
    stopifnot(
        "argument 'dir' must be name/path to the directory with the YAML files" =
            is.character(dir) && length(dir) == 1L && isTRUE(dir.exists(dir))
    )

    ## Searching for yaml files
    yamlfiles <- list.files(dir)[grepl("\\.yaml$", list.files(dir))]
    if (length(yamlfiles) == 0L) stop("No yaml files found in ", dir)

    ## Reading yaml files via a helper function for verbose error messages
    parse_yaml <- function(file, dir) {
        tryCatch(read_yaml(file.path(dir, file)),
            error = function(e) stop("Problems reading \"", file, "\": ", e)) |>
        structure(class = "at_config")
    }

    # Load and return as named list of at_config objects.
    lapply(yamlfiles, parse_yaml, dir = dir) |> setNames(yamlfiles)
}

#' Checking YAML Config
#'
#' Checks if the list (content of the YAML file) is as expected.
#'
#' @param x named list.
#'
#' @details
#' Each YAML file is expected to contain the 'name' of the distribution as well
#' as a section 'parameters > valid' which not only specifies the names of the
#' parameters of the distribution, but also provides a series of numeric values
#' used for testing. If multiple parameters are specified, the numeric vectors
#' must have the same length.
#'
#' Section 'parameters > invalid' (optional) can be specified, defining numeric
#' values which are seen as invalid for specific parameters. Can only be defined
#' for parameters which are included in 'parameters > valid'.
#'
#' The 'support' section (optional) can be used to specify the support of the
#' distribution if constant. When included it must contain 'lower' and 'upper'.
#'
#' A 'tests' section with 'tests > include' or 'tests > exclude' can be specified
#' to only include certain tests, or exclude specific tests. Must evaluate to
#' character vectors.
at_check_config <- function(x) {
    # Basic required elements
    stopifnot(
        "'name' must be string" = is.character(x$name) && length(x$name) == 1L,
        "missing 'parameters' in yaml file" = !is.null(x$parameters),
        "missing 'valid' parameters in yaml file" = !is.null(x$parameters$valid),
        "'parameters > valid' empty in yaml file" = length(x$parameters$valid) > 0L
    )

    ## --------------------------------------
    ## Parameters
    ## --------------------------------------

    ## Valid parameters (required)
    x$parameters$valid <- lapply(x$parameters$valid, function(k) as.numeric(unlist(k)))
    lapply(names(x$parameters$valid), function(n) {
        if (!is.numeric(x$parameters$valid[[n]]) || any(is.na(x$parameters$valid[[n]])))
            stop("'parameters > valid > ", n, "' must be numeric and not NA")
    })
    tmp <- sapply(x$parameters$valid, length)
    if (!all(tmp == tmp[[1L]])) stop("'parameters > valid' must all have the same length")

    ## Invvalid parameters (optional)
    if (!is.null(x$parameters$invalid)) {
        x$parameters$invalid <- lapply(x$parameters$invalid, function(k) as.numeric(unlist(k)))
        lapply(names(x$parameters$invalid), function(n) {
            if (!n %in% names(x$parameters$valid))
                stop("parameter ", n, " in 'parameters > invalid' unknown (not one of the ones listed in 'valid'")
            if (!is.numeric(x$parameters$invalid[[n]]) || any(is.na(x$parameters$invalid[[n]])))
                stop("'parameters > invalid > ", n, "' must be numeric and not NA")
        })
    }

    ## --------------------------------------
    ## Support (optional)
    ## --------------------------------------
    if (!is.null(x$support)) {
        stopifnot(
            "missing 'support > lower'" = !is.null(x$support$lower),
            "missing 'support > upper'" = !is.null(x$support$upper),
            "'support > lower' must be single numeric" =
                is.numeric(x$support$lower) && length(x$support$lower) == 1L,
            "'support > upper' must be single numeric" =
                is.numeric(x$support$upper) && length(x$support$upper) == 1L
        )
    }

    ## --------------------------------------
    ## Enabled/disabled tests (optional)
    ## --------------------------------------
    if (!is.null(x$tests$include))
        stopifnot("'tests > include' incorrectly specified" = is.character(x$tests$include))
    if (!is.null(x$tests$exclude))
        stopifnot("'tests > exclude' incorrectly specified" = is.character(x$tests$exclude))
    return(x)
}


## ----------------------------------------------------------------------------
## Series of functions which create the strings for the automated tests.
## Function `at_setup()` is required to set up a series of tests for the
## constructor functions (required) and will create a series of objects
## placed in the global env to be used by the different tests that follow.
## E.g., two distribution objects with one and multiple distributions, used
## for testing the methods.
##
## All other tests must are added below and can re-use these objects.
## ----------------------------------------------------------------------------

at_clean_cmd <- function(x) {
    ## Stripping first four spaces (if present)
    gsub("^[ ]{4}", "", strsplit(x, "\n")[[1L]])
}

#' Setting up Autotest 
#'
#' @param x object of class `at_config`.
at_setup <- function(x) {
    stopifnot("argument 'x' must be of class 'at_config'" = inherits(x, "at_config"))


    ## -------------------------
    ## Template 1: Tests general functionability of the
    ## constructor function and creates two objects 'd1' and 'dx'
    ## used for testing.
    ## -------------------------

    template1 <- '
    ## --- this file is auto-generated, do not change
    ##
    ## Automated tests for distribution "%1$s"
    ## ... TODO(R)

    ## For interactive use (testing the tests)
    if (interactive()) { library(\"testthat\"); library(\"distributions3\"); }

    ## Test that the constructor function with no input arguments works
    test_that("constructor function %1$s() with default arguments runs silently", {
        ## Constructor function w/ default arguments should run silently
        expect_silent(x <- %1$s())

        ## Returns has correct class attributes
        expect_s3_class(x, "%1$s")
        expect_s3_class(x, "distribution")

        ## Return is either of length 0 or 1
        expect_true(length(x) == 0L || length(x) == 1L, info = "object not of length 0 or 1")

        ## Correct parameters
        expect_identical(sort(names(as.list(x))), sort(%2$s), info = "parameters not as expected/defined")
    });

    ## Test that all parameters defined are arguments to the constructor function
    test_that("constructor %1$s() has arguments %3$s", {
        args <- names(formals(getFromNamespace("%1$s", "distributions3")))
        expect_true(all(!is.na(match(%2$s, args))), info = "not all parameters are arguments to constructor function")
    });

    ## Call constructor function with the pre-defined arguments in parameters$valid
    ## to create two objects for testing: d1 (distribution object with one distribution),
    ## dx (distributions object with as many distributions as defined in parameters$valid).
    ## They are written to the global env the tests run in.
    test_that("constructor %1$s runs silent with parmaeter inputs", {
        expect_silent(d1 <<- %4$s)
        expect_silent(dx <<- %5$s)
    });
    '

    ## Commands to create the test objects d1, dx
    valid <- x$parameters$valid
    cmd_d1 <- sprintf("%s(%s)", x$name,
                      paste(paste(names(valid), "=", sapply(valid, head, n = 1)), collapse = ", "))
    cmd_d1 <<- cmd_d1
    cmd_dx <- sprintf("%s(%s)", x$name,
                      paste(paste(names(valid), "=", sapply(valid, deparse)), collapse = ", "))

    res <- sprintf(
                template1, # template
                x$name,                                             # %1$s - name of distribution
                deparse(names(x$parameters$valid)),                 # %2$s - vector with parameter names
                paste(names(x$parameters$valid), collapse = ", "),  # %3$s
                cmd_d1,                                             # %4$s - cmd to create 'd1'
                cmd_dx                                              # %5$s - cmd to create 'dx'
            )

    return(at_clean_cmd(res))
}




