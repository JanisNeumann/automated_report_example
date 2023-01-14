#
# This file contains functions for splitting data sets
#

split_data_classifier <- function(data, dev_ratio, holdout_ratio, y_var, match_vars = NULL, r_seed = 16) {
  require(caret)
  set.seed(r_seed)

  # Use caret package createDataPartition function if there are no variables
  # to keep similar between training, dev and holdout set besides y_var.
  # Otherwise, go through a more involved process including clustering.
  if (is.null(match_vars)) {
    # Create dev and holdout set iteratively. The createDataPartition can do this in one go
    # but that would enforce equally sized dev and holdout sets which may not always be wanted.
    dev_index <- caret::createDataPartition(data[, y_var],
      p = dev_ratio,
      list = FALSE
    )
    train_set <- data[-dev_index, ]

    holdout_index <- caret::createDataPartition(train_set[, y_var],
      p = holdout_ratio * nrow(data) / nrow(train_set),
      list = FALSE
    )

    holdout_set <- train_set[holdout_index, ]
    train_set <- train_set[-holdout_index, ]
    dev_set <- data[dev_index, ]

    return(list(
      "train_set" = train_set,
      "dev_set" = dev_set,
      "holdout_set" = holdout_set
    ))
  } else {
    # Predefine vectors to be appended to by following loop.
    # Eventually, they contain the row names belonging to each subset of the data.
    train_rows <- c()
    holdout_rows <- c()
    dev_rows <- c()

    for (i in unique(data[, y_var])) {
      # Subset data based on main variable of interest.
      # Retain only match_vars to then cluster by.
      i_data <- data[data[, y_var] == i, match_vars]

      # Clustering followed by creating groups based on the cluster tree.
      i_clust <- hclust(dist(i_data))
      i_groups <- cutree(i_clust, h = max(i_clust$height / 2))

      # Create subsets based on cluster grouping using caret.
      i_dev_index <- caret::createDataPartition(i_groups,
        p = dev_ratio,
        list = FALSE
      )

      i_train <- i_groups[-i_dev_index]

      i_holdout_index <- caret::createDataPartition(i_train,
        p = holdout_ratio * nrow(i_data) / nrow(i_train),
        list = FALSE
      )

      holdout_rows <- c(holdout_rows, names(i_train[i_holdout_index]))
      train_rows <- c(train_rows, names(i_train[-i_holdout_index]))
      dev_rows <- c(dev_rows, names(i_groups[i_dev_index]))
    }

    # Subset input data based on row indices defined in for loop.
    # Return as list containing train, dev and holdout set.
    return(list(
      "train_set" = data[train_rows, ],
      "dev_set" = data[dev_rows, ],
      "holdout_set" = data[holdout_rows, ]
    ))
  }
}
