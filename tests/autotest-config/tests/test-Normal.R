
## "Support" methods
test_that("print.Normal works correctly", {
    d_test_print(Normal(-10, +2))
    d_test_print(Normal(-100.123456, 1e-10))
})

test_that("support.Normal works correctly", {
    ## Generic method tests
    d_test_support(unname(dd), expected_min = -Inf, expected_max = Inf) # unnamed
    d_test_support(dd,         expected_min = -Inf, expected_max = Inf) # named
})

