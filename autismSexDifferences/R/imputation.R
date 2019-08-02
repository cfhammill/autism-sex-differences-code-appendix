#' Convert mids to list of data frames
#'
#' Utility function extract a \code{mids} object into a
#' list of \code{data.frames}
#' @param mids The multiply imputed data set from \link[mice]{mice}
#' @return A list of \code{data.frames}
#' @export
mids_to_dfs <-
    function(mids){
        lapply(seq_len(mids$m), function(i){
            Reduce(function(df,vn){
                df[[vn]][is.na(df[[vn]])] <- mids$imp[[vn]][[i]]
                df
            }
          , names(mids$imp)[!vapply(mids$imp, is.null, logical(1))]
          , init = mids$data)
        })
    }
                           
