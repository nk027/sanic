# Run the tinytest suite when RUN_TESTS is set to a true value, e.g.
#
#   RUN_TESTS=true R CMD check sanic_0.0.2.tar.gz
#
# The variable can also be set in a project-local .Renviron file:
#
#   RUN_TESTS=true
#
# Tests are skipped (rather than failed) when tinytest is unavailable.

run_tests <- tolower(Sys.getenv("RUN_TESTS", "")) %in%
  c("true", "t", "1", "yes", "y")

if (run_tests) {
  if (requireNamespace("tinytest", quietly = TRUE)) {
    tinytest::test_package("sanic")
  } else {
    message("Skipping tests: package 'tinytest' is not installed.")
  }
} else {
  message("Skipping tests: set RUN_TESTS=true to run them.")
}
