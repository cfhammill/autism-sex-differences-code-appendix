#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(RMINC)
  library(dplyr)
  library(purrr)
  library(readxl)
  library(tidyr)
  library(ggplot2)
})

# Read latest summaries
civet_summaries <- 
  read.csv("../data/combined_neuroanatomy20170929.csv"
           , stringsAsFactors = FALSE)

civet_files_CIHR <-
  civet.getAllFilenames(
    filter(civet_summaries, src == "CIHR")
    , idvar = "scan", prefix = "d8_", basedir = "margot_subjs2.1/"
    , civetVersion = "2.1.0"
    , cnf = yaml::yaml.load_file("margot_subjs2.1/d8_0001_01/CBRAIN.params.yml")
  )



civet_paths <- c("../data/civetOutputs/pond_2-1_20170227"
               , "../data/civetOutputs/pond_addition_20170908")

civet_pathed <-
  filter(civet_summaries, src == "POND") %>%
  mutate(scan_np = gsub("-", "_", sub("MR160-", "", scan))
       , civet_path =
           sapply(scan_np, function(subj){
             paths <- file.path(civet_paths, subj)
             for(path in paths)
               if(file.exists(path)) return(dirname(path))

             return(NA)
           })
         )  

civet_files_PND <-
  civet.getAllFilenames(
    filter(civet_pathed, civet_path == civet_paths[1])
  , idvar = "scan_np", prefix = "MR160_", basedir = "../data/civetOutputs/pond_2-1_20170227/"
  , civetVersion = "2.1.0"
  , cnf = yaml::yaml.load_file("../data/civetOutputs/pond_2-1_20170227/088_0002_01_002/CBRAIN_Colosse-384672-1.params.yml")
  )

civet_files_PND_update <-
  civet.getAllFilenames(
    filter(civet_pathed, civet_path == civet_paths[2])
  , idvar = "scan_np", prefix = "MR160_", basedir = "../data/civetOutputs/pond_addition_20170908/"
  , civetVersion = "2.1.0"
  , cnf = yaml::yaml.load_file("../data/civetOutputs/pond_addition_20170908/088_0405_01_002/CBRAIN_Mammouth-453451-1.params.yml")
  )
  

civet_files <- bind_rows(civet_files_CIHR, civet_files_PND, civet_files_PND_update)
civet_files$new_scanner <- civet_files$scan_date > as.Date("2016-06-01")

civet_files <- filter(civet_files, sapply(RSL_mean_curvature_left, file.exists))

iqs <-
    read.csv("../data/combined_iqs20170620.csv") %>%
    filter(!is.na(iq)) %>%
    group_by(visit) %>%
    arrange(desc(scan_date)) %>%
    slice(1) %>%
    ungroup

## Filter
self_intersect_limit <- 150
civet_deduped <- 
  civet_files %>%
  filter(Dx %in% c("ASD", "CTRL")) %>%
  filter(LEFT_SURF_SURF < self_intersect_limit
         , RIGHT_SURF_SURF < self_intersect_limit
         , RIGHT_INTER < self_intersect_limit
         , LEFT_INTER < self_intersect_limit
         , !is.na(age_scan)
         , age_scan < 50) %>%
#  group_by(subject) %>% #remove this stanza to skip delongi
#  slice(`if`(any(best_of_subject), which(best_of_subject), 1)) %>%
#  ungroup %>%
  mutate(passed_old_qc = QC_PASS & best_of_subject) %>%
  left_join(iqs %>% select(visit, iq), by = "visit") %>%
  filter(!duplicated(scan)) %>% ## when merging CIVET dupes can happen
  mutate(RSL_mean_curvature_left = sub("left", "left_abs", RSL_mean_curvature_left)
         , RSL_mean_curvature_right = sub("right", "right_abs", RSL_mean_curvature_right))

model_data <-
  civet_deduped %>%
  mutate(age_scan = as.numeric(scale(age_scan, center = TRUE))
         , bv = as.numeric(scale(BRAIN_VOL, center = TRUE))
         , crbv = as.numeric(scale(BRAIN_VOL^(1/3), center = TRUE))
         , ttrbv = as.numeric(scale(BRAIN_VOL^(2/3), center = TRUE))
         , sex = relevel(factor(sex), "M")
         , iq = as.numeric(scale(iq))
         , Dx = factor(Dx, c("CTRL","ASD")))


## Bring in surfaces
surface_left <- read_obj("../data/CIVET_2.0/CIVET_2.0_icbm_avg_mid_sym_mc_left.obj")
surface_right <- read_obj("../data/CIVET_2.0_icbm_avg_mid_sym_mc_right.obj")

aal_left_mask <-
  readLines("../data/CIVET_2.0/CIVET_2.0_AAL_left.txt") %>%
  as.numeric %>%
  .[1:40962] %>%
  `!=`(0)

aal_right_mask <-
  readLines("../data/CIVET_2.0/CIVET_2.0_AAL_right.txt") %>%
  as.numeric %>%
  .[1:40962] %>%
  `!=`(0)

## Setup vertex tables
lt <- vertexTable(civet_deduped$nativeRMS_RSLtlink_left)
rt <- vertexTable(civet_deduped$nativeRMS_RSLtlink_right)
la <- vertexTable(civet_deduped$midSurfaceleftNativeArea)
ra <- vertexTable(civet_deduped$midSurfacerightNativeArea)
lv <- vertexTable(civet_deduped$SurfaceleftNativeVolume)
rv <- vertexTable(civet_deduped$SurfacerightNativeVolume)

lmc <- vertexTable(civet_deduped$RSL_mean_curvature_left)
rmc <- vertexTable(civet_deduped$RSL_mean_curvature_right)

thickness <- rbind(rt[aal_right_mask,], lt[aal_left_mask,])
area <- rbind(ra[aal_right_mask,], la[aal_left_mask,])
volume <- rbind(rv[aal_right_mask,], lv[aal_left_mask,])
mean_curv <- rbind(rmc[aal_right_mask,], lmc[aal_left_mask,])

## Generate random shuffles so that each imputated set gets the same randomizations
shuffles <- sapply(seq_len(5000), function(i) sample(seq_len(nrow(civet_deduped))))

#del for space
rm(lt, rt, la, ra, lv, rv, lmc, rmc)

## Write out objects
to_date <- format(Sys.Date(), "%Y%m%d")
save.image(
  file.path(".."
            , paste0("cortical_objects_longi"
                     , to_date
                     , ".rda")))
