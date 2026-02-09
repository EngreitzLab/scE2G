# get fill values from config file and compute mean across observations if specified
get_fill_values <- function(features, config) {
  # get fill values for each feature
  fill_values <- config %>%
    filter(feature %in% colnames(features)) %>%
    select(feature, fill_value) %>%
    distinct() %>%
    deframe()

  # get features that use mean as fill value
  mean_fill_values <- fill_values == "mean"
  mean_fill_values[is.na(mean_fill_values)] <- FALSE

  # if there are any, compute mean for these features (if statement technically not needed)
  if (any(mean_fill_values)) {
    # get feature values for features with mean as fill value
    mean_features <- names(fill_values[mean_fill_values])
    mean_features <- features[, ..mean_features]

    # compute mean for all required features while excluding NA, Inf and -Inf values
    mean_values <- apply(mean_features, MARGIN = 2, FUN = function(x) {
      mean(x[!is.infinite(x)], na.rm = TRUE)
    })

    # add mean values to pre-defined fill values
    fill_values[names(mean_values)] <- mean_values
  }

  # convert fill values to list, one element per feature
  fill_values <- as.list(structure(as.numeric(fill_values), names = names(fill_values)))

  return(fill_values)
}
