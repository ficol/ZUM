##' @export
isolationForest <- function(X, ...) UseMethod("isolationForest")

##' Create a simple isolation forest
##' 
##' This creates a simple isolation forest.
##' 
##' @param X  Dataframe of predictors
##' @param num_trees  Number of trees in forest
##' @param sample_size Number of sampled rows to create one tree
##' @param ncols_per_tree Number of sampled columns to create one tree
##' @param max_depth Maximal depth of tree
##' @return isolationForest object from which 'predict' method can be called on new data
##' @export
isolationForest.default <- function(X, num_trees=100, sample_size=256, ncols_per_tree=ncol(X), max_depth=ceiling(log2(sample_size)), ...) {
    sample_size=min(sample_size, nrow(X))
    ncols_per_tree=min(ncols_per_tree, ncol(X))
    iforest <- isolationForestBuild(X, num_trees, sample_size, ncols_per_tree, max_depth)
    iforest$categorical <- sum(sapply(X, is.factor))
    iforest$numerical <- ncol(X) - iforest$categorical
    iforest$call <- match.call()
    class(iforest) <- "isolationForest"
    return(iforest)
}

##' Builds trees for isolation forest
##' 
##' This builds trees for isolation forest based on parameters.
##' 
##' @param X  Dataframe of predictors
##' @param num_trees  Number of trees in forest
##' @param sample_size Number of sampled rows to create one tree
##' @param ncols_per_tree Number of sampled columns to create one tree
##' @param max_depth Maximal depth of tree
##' @return list with vector of isolation trees and parameters
isolationForestBuild <- function(X, num_trees, sample_size, ncols_per_tree, max_depth) {
    X <- data.frame(X)
    forest <- vector("list", length=num_trees)
    for (i in 1:num_trees) {
		forest[i] <- list(itree(sample(X[sample(1:nrow(X), sample_size),], ncols_per_tree), max_depth))
	}
    return(list(forest=forest,
                num_trees=num_trees,
                sample_size=sample_size,
                ncols_per_tree=ncols_per_tree,
                max_depth=max_depth))
}

##' Prints summary information of isolation forest
##'
##' This prints parameters of created isolation forest.
##' @param X Isolation forest object
##' @export
print.isolationForest <- function(X, ...) {
    cat("Call:\n")
    print(X$call)
    cat("Consisting of", X$num_trees, "trees\n")
    cat("Categorical columns:", X$categorical, "\n")
    cat("Numerical columns:", X$numerical, "\n")
}

##' Prints summary information of isolation forest
##'
##' This prints parameters of created isolation forest.
##' @param X Isolation forest object
##' @export
summary.isolationForest <- function(X, ...) {
    print.isolationForest(X, ...)
}

##' Isolation forest formula
##'
##' This implements formula S3 method.
##' @param formula Formula
##' @param data input data
##' @return Isolation forest object
##' @export
isolationForest.formula <- function(formula, data=list(), ...) {
    mf <- model.frame(formula=formula, data=data)
    X <- data.frame(model.matrix(attr(mf, "terms"), data=mf))
    iforest <- isolationForest.default(X, ...)
    iforest$call <- match.call()
    iforest$formula <- formula
    return(iforest)
}

##' Predicts method for isolation forest
##'
##' This predicts values on new data using isolation forest.
##' @param object isolation forest model
##' @param newdata data to predict
##' @return array of predicted values from range (0,1)
##' @export
predict.isolationForest <- function(object, newdata, ...) {
    newdata <- data.frame(newdata)
    path_lengths <- apply(newdata, 1, function(row){
        mean(sapply(object$forest, function(tree){
            path_length(row, tree)
        }))
    })
    return(2^(-path_lengths/c_value(object$sample_size)))
}

##' Builds isolation tree
##'
##' This builds isolation tree based on parameters. If chosen variable is continual, split value is random value from min to max.
##' If chosen variable is discrete, split value is a subset of possible values.
##' Node becomes leaf if size of samples is <= 1, every sample is identical or node's depth equals max depth.
##' @param X data to create tree
##' @param max_depth maximal depth of tree
##' @param curr_depth current depth of tree
##' @return Node of isolation tree
itree <- function(X, max_depth, curr_depth=0) {
    if (curr_depth >= max_depth || nrow(X) <= 1 || nrow(unique(X)) == 1) {
        return(list(is_leaf=TRUE, size=nrow(X)))
    }
    q <- names(X)[sample(1:ncol(X), 1)]
    if (sapply(X[q], is.factor)) {
        is_categorical <- TRUE
        p <- sample(unique(X[[q]]), sample(1:nrow(unique(X[q])), 1))
        X_l <- X[!(X[[q]] %in% p),]
        X_r <- X[X[[q]] %in% p,]
    }
    else {
        is_categorical <- FALSE
        p <- runif(1, min(X[q]), max(X[q]))
        X_l <- X[X[q] < p,]
        X_r <- X[X[q] >= p,]
    }
    return(list(is_leaf=FALSE,
                left=itree(X_l, max_depth, curr_depth + 1),
                right=itree(X_r, max_depth, curr_depth + 1),
                is_categorical=is_categorical,
                split_attr=q,
                split_val=p))
}

##' Calculates the path length of sample
##'
##' This calculates the path length of sample in tree 
##' @param x data sample
##' @param T isolation tree
##' @param curr_length current length of path
##' @return length of path from root to external node
path_length <- function(x, T, curr_length=0)
{
    if (T$is_leaf) {
        return(curr_length + c_value(T$size))
    }
    if ((T$is_categorical && !(x[T$split_attr] %in% T$split_val)) || (!T$is_categorical && x[T$split_attr] < T$split_val)) {
        return(path_length(x, T$left, curr_length + 1))
    }
    return(path_length(x, T$right, curr_length + 1))
}

##' Calculates c value
##'
##' This calculates c value based on size of node
##' @param n size of node
##' @return c value of node
c_value <- function(n) {
    if (n < 2) {
        value <- 0
    }
    else if (n == 2) {
        value <- 1
    }
    else {
        value <- 2 * (log(n - 1) + 0.5772156649) - (2 * (n - 1) / n)
    }
    return(value)
}
