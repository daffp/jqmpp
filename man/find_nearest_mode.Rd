% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_mode.R
\name{find_nearest_mode}
\alias{find_nearest_mode}
\title{Find the nearest maximum in an array}
\usage{
find_nearest_mode(A, ...)
}
\arguments{
\item{A}{An array.}

\item{...}{The coordinates of a point within the dimensions of the array.}
}
\value{
Returns a list of the array positions of the nearest maximum that was found.
}
\description{
Find the nearest maximum in an array to a given point by hillclimbing.
}
\details{
A hillclimbing search from the input point across the array to find the nearest local maximum. 
The search is through each unit n-cube surrounding the next input point in the search (accounting for the array boundaries).
Ties are broken by taking the nearest euclidean distance (i.e. point which are vertical or horizontal from the input rather than diagonal)
and then random sampling. To make this reproducible the seed should be set.
}
\examples{
set.seed(444)
msize <- 20
A <- array(sample(seq_len(msize), msize^3, replace=TRUE), dim = rep(msize, 3))
find_nearest_mode(A, 3,3,3)
find_nearest_mode(A, 4,1,5)
}
