
#' Reading/Parsing YAML File
#'
#' @param dir character. Name of the directory containing the `.yaml` files.
#'
#' @reutrn Named list with the content of the YAML files. The names of the
#' list elements corresponds to the YAML file name.
autotest_load_config <- function(dir = "autotest-config") {
    stopifnot(
        "argument 'dir' must be name/path to the directory with the YAML files" =
            is.character(dir) && length(dir) == 1L && isTRUE(dir.exists(dir))
    )

    ## Searching for yaml files
    yamlfiles <- list.files(dir)[grepl("\\.yaml$", list.files(dir))]
    if (length(yamlfiles) == 0L) stop("No yaml files found in ", dir)

    ## Reading yaml files via a helper function for verbose error messages
    parse_yaml <- function(file, dir) {
        res <- tryCatch(read_yaml(file.path(dir, file)),
            error = function(e) stop("Problems reading \"", file, "\": ", e))
        ## Appending 'testfile' (the output file name), add class, return
        res$testfile <- sprintf("test-autotest-%s", gsub("\\.yaml$", ".R", file))
        res |> structure(class = "autotest_config")
    }

    # Load and return as named list of autotest_config objects.
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
autotest_check_config <- function(x) {
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
    stopifnot(
        "'parameters > valid' must all have length 1L or N" = all(tmp == 1L | tmp == max(tmp)),
        "'parameters > valid' must at least have N = 2" = max(tmp) > 1L
    )
    ## Extending parameter vectors if needed
    x$parameters$valid <- lapply(x$parameters$valid, function(x) rep_len(x, max(tmp)))

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
    if (!is.null(x$tests)) {
        ## Allows to have 'include' or 'exclude' but mutually exclusive.
        stopifnot(
            "'tests' only allow to include one include or exclude rule" =
                all(names(x$tests) %in% c("include", "exclude")) && length(x$tests) == 1L,
            "'tests > include' or 'tests > exclude' must be character vector of length > 1" =
                is.character(x$tests[[1L]]) && length(x$tests[[1L]]) > 0L
        )
    }

    return(x)
}


## ----------------------------------------------------------------------------
## Series of functions which create the strings for the automated tests.
## Function `autotest_setup()` is required to set up a series of tests for the
## constructor functions (required) and will create a series of objects
## placed in the global env to be used by the different tests that follow.
## E.g., two distribution objects with one and multiple distributions, used
## for testing the methods.
##
## All other tests must are added below and can re-use these objects.
## ----------------------------------------------------------------------------

clean_cmd <- function(x) {
    ## Stripping first four spaces (if present)
    gsub("^[ ]{4}", "", strsplit(x, "\n")[[1L]])
}

replace_placeholders <- function(template, x, indentation = 8) {
    indentation <- as.integer(indentation)[[1L]]
    stopifnot(is.character(template), is.list(x), !is.null(names(x)),
              is.integer(indentation) && indentation > 0L)

    ## Search and replace placeholders
    for (i in seq_along(x)) {
        template <- gsub(names(x)[i], x[[i]], template)
    }

    ## Adjusting indentation
    template <- gsub("[ ]{0,}<INDENT>", paste0(rep(" ", indentation), collapse = ""), template)
    return(template)
}

#' Setting up Autotest 
#'
#' @param x object of class `autotest_config`.
#'
#' @return Character vector with commands for test set.
autotest_setup <- function(config) {
    stopifnot("argument 'config' must be of class 'autotest_config'" = inherits(x, "autotest_config"))

    template <- '
    ## --- this file is auto-generated, do not change
    ##
    ## Automated tests for distribution "%1$s"
    ## ... TODO(R)

    ## For interactive use (testing the tests)
    if (interactive()) { library(\"testthat\"); library(\"distributions3\"); }

    ## Test that the constructor function with no input arguments works
    test_that("constructor function <DIST>() with default arguments runs silently", {
        ## Constructor function w/ default arguments should run silently
        expect_silent(x <- <DIST>())

        ## Returns has correct class attributes
        expect_s3_class(x, c("<DIST>", "distribution"), exact = TRUE)

        ## Return is either of length 0 or 1
        expect_true(length(x) == 0L || length(x) == 1L, info = "object not of length 0 or 1")

        ## Correct parameters
        expect_identical(sort(names(as.list(x))), sort(<PARAMVEC>), info = "parameters not as expected/defined")
    });

    ## Test that all parameters defined are arguments to the constructor function
    test_that("constructor <DIST>() has arguments <PARAMTXT>", {
        args <- names(formals(getFromNamespace("<DIST>", "distributions3")))
        expect_true(all(!is.na(match(<PARAMVEC>, args))), info = "not all parameters are arguments to constructor function")
    });

    ## Call constructor function with the pre-defined arguments in
    ## parameters$valid to create four objects for testing: d0 (default
    ## return), d1 (unnamed distribution object with one distribution), dx
    ## (unnamed distributions object with as many distributions as defined in
    ## parameters$valid), and dn (named version of dx). They are written to the
    ## global env the tests run in.
    test_that("constructor %1$s runs silent with parameter inputs", {
        expect_silent(d0 <<- <DIST>())
        expect_silent(d1 <<- <CMD_D1>)
        expect_silent(dx <<- <CMD_DX>)
        expect_silent(dn <<- <CMD_DN>)
    });
    '

    ## Commands to create the test objects d1, dx, dn
    valid <- x$parameters$valid
    cmd_d1 <- sprintf("%s(%s)", config$name,
                      paste(paste(names(valid), "=", sapply(valid, head, n = 1)), collapse = ", "))
    cmd_d1 <- sprintf("%s(%s)", config$name,
                      paste(paste(names(valid), "=", sapply(valid, head, n = 1)), collapse = ", "))
    cmd_dx <- sprintf("%s(%s)", config$name,
                      paste(paste(names(valid), "=", sapply(valid, deparse)), collapse = ", "))

    ## Same as cmd_dx but with named parameters, which should result in a named distributions object
    names(valid[[1L]]) <- LETTERS[seq_along(valid[[1L]])] ## Naming first parameter set
    cmd_dn <- sprintf("%s(%s)", config$name,
                      paste(paste(names(valid), "=", sapply(valid, deparse)), collapse = ", "))

    ## List of values/placeholders to replace
    p <- list("<DIST>"       = config$name,
              "<PARAMVEC>"   = deparse(names(config$parameters$valid)),
              "<PARAMTXT>"   = paste0(names(config$parameters$valid), collapse = ", "),
              "<CMD_D1>"     = cmd_d1,
              "<CMD_DX>"     = cmd_dx,
              "<CMD_DN>"     = cmd_dn)

    return(replace_placeholders(template, p) |> clean_cmd())
}


#' Setting up Autotests for Constructor Return
#'
#' Testing correctness of the return of the constructor function. This
#' test set expects that we already have the two named distribution objects
#' `d1` and `dx` in the global environment as created by `autotest_setup()` which
#' must be added before this test set.
#'
#' @param x object of class `autotest_config`.
#'
#' @return Returns `NULL` if the test is excluded, or a character vector
#' containing the test sets (R commands).
autotest_constructor <- function(config) {
    ## Sanity check + check if this test should be included
    stopifnot("argument 'config' must be of class 'autotest_config'" = inherits(x, "autotest_config"))

    template <- '
    ## Testing object `d1`: unnamed object of length 1
    test_that("constructor returns expected object: unnamed, length 1", {
        ## Class
        expect_s3_class(d0, c("<DIST>", "distribution"), exact = TRUE)
        expect_s3_class(d1, c("<DIST>", "distribution"), exact = TRUE)
        expect_s3_class(dn, c("<DIST>", "distribution"), exact = TRUE)
        expect_s3_class(dx, c("<DIST>", "distribution"), exact = TRUE)

        ## Converting to data frame for testing
        dd0 <- as.data.frame(as.matrix(d0))
        dd1 <- as.data.frame(as.matrix(d1))
        ddx <- as.data.frame(as.matrix(dx))
        ddn <- as.data.frame(as.matrix(dn))

        ## Dimension and names/dimension names
        expect_identical(dim(dd1),      c(1L, NCOL(as.matrix(d0))))
        expect_identical(names(dd1),    names(dd0))
        expect_null(names(d1))

        expect_identical(dim(ddx),      c(<N>, NCOL(as.matrix(d0))))
        expect_identical(names(ddx),    names(dd0))
        expect_null(names(dx))

        expect_identical(dim(ddn),      c(<N>, NCOL(as.matrix(d0))))
        expect_identical(names(ddn),    names(ddn))
        expect_identical(names(dn),     LETTERS[seq_len(<N>)])
        expect_identical(rownames(ddn), LETTERS[seq_len(<N>)])

        ## Checking values
        <CMD_VALUES>
    });
    '

    ## Generating commands to test the values of the parameters
    cmd_values <- with(config$parameters,
        sprintf("<INDENT>expect_equal(ddn$%s, %s)", names(valid), sapply(valid, deparse)))

    ## List of values/placeholders to replace
    p <- list("<DIST>"       = x$name,
              "<N>"          = deparse(length(config$parameters$valid[[1L]])),
              "<CMD_VALUES>" = paste(cmd_values, collapse = "\n"))

    return(replace_placeholders(template, p) |> clean_cmd())
}


#' Setting up Autotests for S3 Method `as.list`
#'
#' Testing function of the S3 method 'as.list()'.
#' Expects object `dn` to exist (see `autotest_setup()`).
#'
#' @param x object of class `autotest_config`.
#'
#' @return Character vector with commands for test set.
autotest_support <- function(x) {
    stopifnot("argument 'x' must be of class 'autotest_config'" = inherits(x, "autotest_config"))

    template <- '
    test_that("as.list method works as expected", {
        ## Names of parameters
        expect_silent(ln <- as.list(dn))

        ## Comparing against expected list
        expect_identical(ln, <DIST>)
    });
    '
    ## List of values/placeholders to replace
    p <- list("<DIST>"       = x$name,
              "<N>"          = deparse(length(config$parameters$valid[[1L]])))

    return(replace_placeholders(template, p) |> clean_cmd())
}


