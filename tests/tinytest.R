### "lazy" transition from "testthat" to "tinytest"
### Copyright (C) 2020-2021 Sebastian Meyer

if (!requireNamespace("tinytest", quietly = TRUE)
    || packageVersion("tinytest") < "1.2.4") {
    message("this test suite requires package 'tinytest' (>= 1.2.4)")
    q("no")
}

## provide simple replacement for test_that() expectation bundles
test_that <- function (desc, code) {
    eval(substitute(code), new.env(parent = parent.frame()))
    invisible()
}

## we use verbose = 1 to print_status() only after each test file,
## not after each expression (verbose = 2)
tinytest::test_package("surveillance", testdir = "testthat", verbose = 1)
