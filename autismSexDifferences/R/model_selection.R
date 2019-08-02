#' Fit a common suite of anatLms for model selection
#'
#' Fit all models to a common anat object and covariate data.frame
#' for use in model selection
#'
#' @param form_list A named list of one sided formulae for model selection
#' @param anat A matrix of anatomy values
#' @param model_data A data frame of covariates
#' @return A list of \code{anatModel} objects for model selection
#' @export
fit_models <-
    function(form_list =
                 c(nint = ~ sex + Dx + bv + age_scan + Dx:sex + new_scanner
                 , nint_crb = ~ sex + Dx + crbv + age_scan + Dx:sex + new_scanner
                 , nint_ttrb = ~ sex + Dx + ttrbv + age_scan + Dx:sex + new_scanner
                 , pnint = ~ sex + Dx + bv + poly(age_scan,2) + Dx:sex + new_scanner
                 , pnint_crb = ~ sex + Dx + crbv + poly(age_scan,2) + Dx:sex + new_scanner
                 , pnint_ttrb = ~ sex + Dx + ttrbv + poly(age_scan,2) + Dx:sex + new_scanner
                 , nex = ~ sex + Dx + bv + age_scan + Dx:sex + age_scan:sex + new_scanner
                 , nex_crb = ~ sex + Dx + crbv + age_scan + Dx:sex + age_scan:sex + new_scanner
                 , nex_ttrb = ~ sex + Dx + ttrbv + age_scan + Dx:sex + age_scan:sex + new_scanner
                 , nage = ~ sex + Dx + bv + age_scan + Dx:sex + bv:sex + new_scanner
                 , nage_crb = ~ sex + Dx + crbv + age_scan + Dx:sex + crbv:sex + new_scanner
                 , nage_ttrb = ~ sex + Dx + ttrbv + age_scan + Dx:sex + ttrbv:sex + new_scanner
                 , st = ~ sex + Dx + bv + age_scan + Dx:sex + bv:sex + age_scan:sex + new_scanner
                 , crb = ~ sex + Dx + crbv + age_scan + Dx:sex + crbv:sex + age_scan:sex + new_scanner
                 , ttrb = ~ sex + Dx + ttrbv + age_scan + Dx:sex + ttrbv:sex + age_scan:sex + new_scanner
                 , page = ~ sex + Dx + bv + poly(age_scan,2) + Dx:sex + bv:sex + poly(age_scan,2):sex + new_scanner
                 , page_crb = ~ sex + Dx + crbv + poly(age_scan,2) + Dx:sex + crbv:sex + poly(age_scan,2):sex + new_scanner
                 , page_ttrb = ~ sex + Dx + ttrbv + poly(age_scan,2) + Dx:sex + ttrbv:sex + poly(age_scan,2):sex + new_scanner)
           , anat, model_data, ...){

        # Fit all the models in the model list
        capture.output(
            models <- lapply(form_list
                           , function(form) anatLm(form
                                                 , data = model_data
                                                 , anat = anat
                                                 , ...))
        )
        attr(models, "form_list") <- form_list

        models
    }


#' Compute evidence ratios
#'
#' Compute the ratio of evidence for the best model relative to each competing model
#'
#' @param aics Vector of AIC values
#' @return A vector of evidence ratios
#' @export
evidence_ratio <- function(aics) exp(1/2 * (aics - min(aics)))

#' Model ranking
#'
#' Take a list of anatLms and rank the models based on minimum median evidence ratio
#' @param model_list A list of anatLm models
#' @param form_list A named list of formulae, extracted from the anatLm list if not supplied and
#' possible
#' @return A \code{tibble} with the model names, formulae, and evidence ratios
#' @export
rank_models <-
    function(model_list, form_list = NULL){
        if(is.null(form_list))
            form_list <- attr(model_list, "form_list")

        if(is.null(form_list))
            stop("form_list not supplied and not present as an attribute of model_list")

        aics <- sapply(model_list, AIC)
        rank_AICs(aics, form_list) 
    }

#' Rank AICs
#'
#' Take a set of AICS and produce a nice ranking tibble
#' @inheritParams rank_models
#' @return A \code{tibble} with the model names, formulae, and evidence ratios
#' @export
rank_AICs <-
    function(aics, form_list){
        evr <- apply(aics, 1, evidence_ratio)
        med_evr <- apply(evr, 1, median)

        data_frame(model = names(form_list)
                 , median_evidence_ratio = med_evr
                 , formula = form_list) %>%
            arrange(median_evidence_ratio)
    }
