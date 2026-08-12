suppressMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
})

## testthat::test_dir() runs each test file with its own directory as the
## working directory, so anchor back up to the repo root.
source(file.path("..", "..", "..", "workflow", "scripts", "model_application", "get_fill_values.R"))

test_that("get_fill_values uses config value for non-mean features and computes mean otherwise", {
  config <- tribble(
    ~feature, ~fill_value,
    "A", "0",
    "B", "mean",
    "C", "5"
  )
  features <- data.table(A = c(1, 2, 3), B = c(10, 20, 30), C = c(NA, 5, 7))

  result <- get_fill_values(features, config)

  expect_equal(result, list(A = 0, B = 20, C = 5))
})

test_that("get_fill_values excludes NA and Inf values when computing the mean", {
  config <- tribble(
    ~feature, ~fill_value,
    "A", "mean"
  )
  features <- data.table(A = c(1, NA, Inf, -Inf, 3))

  result <- get_fill_values(features, config)

  expect_equal(result, list(A = 2))
})
