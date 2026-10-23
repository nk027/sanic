# Run by default, including on CRAN. Set RUN_TESTS=false to opt out locally.
# Tests are skipped when the suggested package tinytest is unavailable.

run_tests <- tolower(Sys.getenv("RUN_TESTS", "true")) %in%
  c("true", "t", "1", "yes", "y")

if (run_tests) {
  if (requireNamespace("tinytest", quietly = TRUE)) {
    tinytest::test_package("sanic")
  } else {
    message("Skipping tests: package 'tinytest' is not installed.")
  }
} else {
  message("Skipping tests: RUN_TESTS is disabled.")
}
