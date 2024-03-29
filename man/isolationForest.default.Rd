% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iforest.R
\name{isolationForest.default}
\alias{isolationForest.default}
\title{Create a simple isolation forest}
\usage{
\method{isolationForest}{default}(
  X,
  num_trees = 100,
  sample_size = 256,
  ncols_per_tree = ncol(X),
  max_depth = ceiling(log2(sample_size)),
  ...
)
}
\arguments{
\item{X}{Dataframe of predictors}

\item{num_trees}{Number of trees in forest}

\item{sample_size}{Number of sampled rows to create one tree}

\item{ncols_per_tree}{Number of sampled columns to create one tree}

\item{max_depth}{Maximal depth of tree}
}
\value{
isolationForest object from which 'predict' method can be called on new data
}
\description{
This creates a simple isolation forest.
}
