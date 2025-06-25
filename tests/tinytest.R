### "lazy" transition from "testthat" to "tinytest"
### Copyright (C) 2020-2022 Sebastian Meyer

if (!requireNamespace("tinytest", quietly = TRUE)
    || packageVersion("tinytest") < "1.4.1") {
    message("this test suite requires package 'tinytest' (>= 1.4.1)")
    q("no")
}

## provide simple replacement for test_that() expectation bundles
## WARNING: this wrapper doesn't print test results, not even failures!
test_that <- function (desc, code) {
    eval(substitute(code), new.env(parent = parent.frame()))
    invisible()
}

## show warnings as they appear
options(warn = 1)

## use verbose = 1 to print_status() only after each test file,
## not after each expression (verbose = 2),
## and omit ANSI escapes for a cleaner log
tinytest::test_package("surveillance", testdir = "testthat",
                       verbose = 1, color = FALSE)
