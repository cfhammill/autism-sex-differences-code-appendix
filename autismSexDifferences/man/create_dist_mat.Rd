% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matching.R
\name{create_dist_mat}
\alias{create_dist_mat}
\title{Construct a distance matrix from a data frame}
\usage{
create_dist_mat(data, w_iq = 1/2, w_bv = 1/4)
}
\arguments{
\item{data}{The data frame of interest}

\item{w_iq}{The weight for iq (as a fraction of the weight for age}

\item{w_bv}{The weight for brain-volume as a fraction of the weight
for age}
}
\value{
An n x n matrix of distances where n is the number of rows
in data. Row and colnames are set to the rownames of data
}
\description{
Apply \link{subj_dist} to each pair of rows in
a data.frame, yeilds results compatible with matchit
}
