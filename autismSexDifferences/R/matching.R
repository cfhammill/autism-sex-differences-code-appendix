#' Custom matching distance function
#'
#' Used for the OHBM analysis, exact matching on sex
#' euclidean distance for age after scaling, half
#' weight for IQ, quarter weight for brain volume by
#' default.
#'
#' @param data A data frame containing Dx, sex, bv, and
#' iq columns
#' @param w_iq The weight for iq (as a fraction of the weight for age
#' @param w_bv The weight for brain-volume as a fraction of the weight
#' for age
#' @return The distance between subjects x and y
#' @export
subj_dist <- function(i,j,data, w_iq = 1/2, w_bv = 1/4){
    with(data, {
        if(sex[i] != sex[j]) return(Inf)
        return(sqrt((age_scan[i] - age_scan[j])^2 +
                    (w_iq * (iq[i] - iq[j]))^2 +
                    (w_bv * (bv[i] - bv[j]))^2
                    ))
    })
}

#' Construct a distance matrix from a data frame
#'
#' Apply \link{subj_dist} to each pair of rows in
#' a data.frame, yeilds results compatible with matchit
#'
#' @param data The data frame of interest
#' @param w_iq The weight for iq (as a fraction of the weight for age
#' @param w_bv The weight for brain-volume as a fraction of the weight
#' for age
#' @return An n x n matrix of distances where n is the number of rows
#' in data. Row and colnames are set to the rownames of data
#' @export
create_dist_mat <-
    function(data, w_iq = 1/2, w_bv = 1/4){
        dmat <- sapply(seq_len(nrow(data))
                     , function(i) sapply(seq_len(nrow(data))
                                        , function(j) subj_dist(i,j,data, w_iq, w_bv)))

        rownames(dmat) <- rownames(data)
        colnames(dmat) <- rownames(data)

        dmat
    }

#' Perform full matching
#'
#' Given a formula, a distance matrix and a data.frame
#' perform full matching on the rows of data. Returns
#' the data frame augmented with \code{weight} and
#' \code{subclass} columns
#' @param formula The formula for the matching
#' @param data The data.frame of interest
#' @param dist_mat The n x n distance matrix
#' @param ... Additional arguments to matchit and optmatch below.
#' @return data augmented with \code{weight} and
#' \code{subclass} columns
#' @export
full_match <-
    function(data
           , dist_mat
           , formula = Dx ~ age_scan + bv + iq + sex
           , ...){
        matches <-
            matchit(Dx ~ sex + age_scan + bv + iq
                  , data = data
                  , method = "full"
                  , distance = dist_mat
                  , ...)

        data %>%
            mutate(weight = matches[["weights"]]
                 , subclass = matches[["subclass"]])
    }
