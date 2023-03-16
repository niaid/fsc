# baseline nettwork correlation 
set.seed(1990)
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(scglmmr))

# Day 7 array data and modules in day 7 data  # load day 7 gene sets 
sig7 = readRDS("signature_curation/core_d7.rds")
module.list = list("CHI_Day7_Response" = sig7$`CHI d7 Response`)

datapath = here("mid_res/baseline_response/dataV3/")

# read CHI array data 
subjects = data.table::fread(here('data/full_metadata/full_sample_metadata.txt')) %>% 
  filter(CITEdata == 1 & vaccine_cohort == 'H1N1') %$% 
  subjectid

array7 = 
  data.table::fread(here("data/CHI_H1N1_data/microarray/CHI_GE_matrix_gene.txt"), data.table = F) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames("gene") %>% 
  select(which(substr(names(.),1,3) %in% subjects)) %>% 
  select(., matches("day0|day1|day7")) %>% 
  select(-matches(c("day70|pre|day1") )) %>% 
  select(-c('207_day0', '209_day0', '279_day0'))

# calculate microarray fold changge (data is already in log space)
t0 =  array7[ ,seq(from=1, to = ncol(array7), by = 2)]
t1 = array7[ ,seq(from=2, to = ncol(array7), by = 2)]
stopifnot(str_sub(colnames(t1), 1,3) == str_sub(colnames(t0), 1,3))
stopifnot(dim(t0) == dim(t1))
fc7 = t1 - t0
fc7 = as.data.frame(fc7)

# calculate the average z score of the day 7 fold change values across samples
d7res = calc_avg_module_zscore(module.list = module.list, average.data.frame = fc7)
names(d7res) = str_sub(names(d7res), 1,3)

# save results 
saveRDS(d7res, file = paste0(datapath, 'd7res.rds'))

sessionInfo()



