suppressPackageStartupMessages({
    library(dplyr)
    library(lenses)
    library(purrr)
    library(batchtools)
    library(tidyr)
})

files <- read.csv("civet_deduped_with_tlaplace.csv", stringsAsFactors = FALSE)

get_vals <- function(x) as.numeric(readLines(x))

process_pair <- function(link, lap){
    link <- get_vals(link)
    lap <- get_vals(lap)    
    list(cor = cor(link, lap)
       , mean_link = mean(link)
       , mean_lap = mean(lap)
       , sd_link = sd(link)
       , sd_lap = sd(lap)
       , mean_dif = mean(abs(link - lap))
       , sd_dif = sd(abs(link - lap)))
}

correlations <-
  files %>%
  mutate(left_cor = map2(nativeRMStlink_left, left_tlap_rs, process_pair)
       , right_cor = map2(nativeRMStlink_right, right_tlap_rs, process_pair))


summary_frame <- 
    correlations %>%
    select(subject:age_scan, left_cor, right_cor) %>%
    gather(side, comp, left_cor, right_cor) %>%
    mutate(comp = map(comp, as_data_frame)) %>%
    unnest

summary_frame %>%
    summarize_each(mean, cor:sd_dif)

tlink <- RMINC::vertexTable(with(files, c(nativeRMStlink_right, nativeRMStlink_left)))
tlap <- RMINC::vertexTable(with(files, c(right_tlap_rs, left_tlap_rs)))

## for_modelling <-
##     rbind(RMINC::vertexTable(files$right_tlap_rs), RMINC::vertexTable(files$left_tlap_rs))
              
## saveRDS(for_modelling, "tlaplace-subjects.rds")

vertex_cors <- sapply(seq_len(ncol(tlap)), function(i) cor(tlap[,i], tlink[,i]))

