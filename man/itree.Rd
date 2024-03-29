% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/iforest.R
\name{itree}
\alias{itree}
\title{Builds isolation tree}
\usage{
itree(X, max_depth, curr_depth = 0)
}
\arguments{
\item{X}{data to create tree}

\item{max_depth}{maximal depth of tree}

\item{curr_depth}{current depth of tree}
}
\value{
Node of isolation tree
}
\description{
This builds isolation tree based on parameters. If chosen variable is continual, split value is random value from min to max.
If chosen variable is discrete, split value is a subset of possible values.
Node becomes leaf if size of samples is <= 1, every sample is identical or node's depth equals max depth.
}
