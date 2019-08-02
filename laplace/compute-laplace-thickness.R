


suppressPackageStartupMessages({
    library(dplyr)
    library(lenses)
    library(purrr)
    library(batchtools)
    library(tidyr)
})

load("../setup/cortical_objects_longi20180115.rda")

files <-
  civet_deduped %>%
  select(subject:age_scan
       , nativeRMStlink_left
       , nativeRMStlink_right
         ) %>%
  mutate(prefix = sub("thick.*", "", nativeRMStlink_left)
       , subject =
           sub(".*thickness/", "", nativeRMStlink_left) %>%
           sub("_native.*", "", .)
       , left_w_surface = paste0(prefix, "surfaces/", subject, "_white_surface_left_81920.obj")
       , left_g_surface = paste0(prefix, "surfaces/", subject, "_gray_surface_left_81920.obj")
       , left_m_surface = paste0(prefix, "surfaces/", subject, "_mid_surface_left_81920.obj")
       , right_w_surface = paste0(prefix, "surfaces/", subject, "_white_surface_right_81920.obj")
       , right_g_surface = paste0(prefix, "surfaces/", subject, "_gray_surface_right_81920.obj")
       , right_m_surface = paste0(prefix, "surfaces/", subject, "_mid_surface_right_81920.obj") 
       , left_surfmap = paste0(prefix, "transforms/surfreg/", subject, "_left_surfmap.sm")
       , right_surfmap = paste0(prefix, "transforms/surfreg/", subject, "_right_surfmap.sm")
       , left_tlap = paste0("thickness/", subject, "_left_tlaplace_30mm.txt")
       , right_tlap = paste0("thickness/", subject, "_right_tlaplace_30mm.txt")
       , left_tlap_rs = paste0("thickness/", subject, "_left_tlaplace_resampled_30mm.txt")
       , right_tlap_rs = paste0("thickness/", subject, "_right_tlaplace_resampled_30mm.txt")
         ) 

cmds <-
  files %>%
  mutate(calc_right_thick = paste0("cortical_thickness -fwhm 30 -tlaplace ", right_w_surface, " ", right_g_surface, " ", right_tlap)
       , calc_left_thick = paste0("cortical_thickness -fwhm 30 -tlaplace ", left_w_surface, " ", left_g_surface, " ", left_tlap)
       , resample_right_thick = paste0("surface-resample ../data/CIVET_2.1/Linux-x86_64/CIVET-2.1.0/models/icbm/icbm_avg_mid_sym_mc_right.obj  ", right_m_surface, " ", right_tlap, " ", right_surfmap, " ", right_tlap_rs)
       , resample_left_thick = paste0("surface-resample ../data/CIVET_2.1/Linux-x86_64/CIVET-2.1.0/models/icbm/icbm_avg_mid_sym_mc_left.obj  ", left_m_surface, " ", left_tlap, " ", left_surfmap, " ", left_tlap_rs)         
        ) 

transposed_commands <-
  cmds %>% purrr::transpose()

reg <- makeRegistry("tlap-reg")

reg$cluster.functions <-
  makeClusterFunctionsTORQUE("torque.tmpl")

reg$default.resources <-
  list(nodes = 1
     , memory = "4G"
     , walltime = "4:00:00")

reg$packages <- "dplyr"

jobs <- batchMap(function(d){
  d$calc_right_thick[[1]] %>% system
  d$calc_left_thick[[1]] %>% system
  d$resample_right_thick[[1]] %>% system
  d$resample_left_thick[[1]] %>% system
}, transposed_commands)

chunks <- chunk(jobs$job.id, chunk.size = 5)
jobs$chunk <- chunks

submitJobs(jobs)
