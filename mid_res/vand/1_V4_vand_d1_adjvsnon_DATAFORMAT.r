# R version 4.0.5 
# H5 vs H1 day 1 DE contrast pseudobulk model 
# reanalysis of Howard et. al AS03 vs PBS cotrol with sorted lineages rnaseq data
# data from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0167488#pone.0167488.s009
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

###### save paths  
datapath = here("mid_res/vand/generated_data/"); dir.create(datapath)

 # load log CPM data from supp table 2 list of celltypes 
vand_datapath = here("data/vand/data/")
ctd = list.files(path = vand_datapath, full.names = T)
ctd = ctd[-1]
e = lapply(ctd, function(x){read_delim(x, delim = ",")})

# get celltypes to name list 
celltypes = lapply(e,function(x){ str_sub(colnames(x), 3,5) %>% unique}) %>% unlist
celltypes = setdiff(celltypes, "SEM")
names(e) = celltypes

# setup data for limma / dream lme4
e = lapply(e, function(x){
  x = x %>%
    select(-ENSEMBL63_GENE_ID) %>% 
    select(matches('D000|D001|ENSEMBL63_GENE_NAME')) %>% 
    column_to_rownames("ENSEMBL63_GENE_NAME") 
})


# create metadata 
md = lapply(e, function(x){
  colnames(x) %>%
    as.data.frame() %>% 
    mutate(group = str_sub(., -4,-1)) %>% 
    mutate(group = if_else(group == "_PBS", "xPBS", group )) %>% 
    mutate(subjectid = str_sub(., 1,1)) %>% 
    mutate(timepoint = str_sub(., 7,10)) %>% 
    mutate(timepoint = ifelse(timepoint == 'D000', 'd0', 'd1')) %>% 
    mutate(time.group = paste(timepoint, group,sep = "_")) %>% 
    column_to_rownames(".")
})

# rm samples with missing data.
lapply(md, function(x) { x$subjectid})
### 
# md 1 = P 
# md 5 = P 
# md 6 = C J 

## remove missing data from sample data 
md[[1]] = md[[1]] %>% 
  rownames_to_column("sample") %>% 
  filter(!subjectid %in% "P") %>%
  column_to_rownames("sample")
md[[5]] = md[[5]] %>%
  rownames_to_column("sample") %>% 
  filter(!subjectid %in% "P") %>%
  column_to_rownames("sample")
md[[6]] = md[[6]] %>% 
  rownames_to_column("sample") %>% 
  filter(!subjectid %in% c("C", "J")) %>%
  column_to_rownames("sample")

#  remove missing data from RNAseq data 
e[[1]] = e[[1]][ ,rownames(md[[1]])]
e[[5]] = e[[5]][ ,rownames(md[[5]])]
e[[6]] = e[[6]][ ,rownames(md[[6]])]

#############################################
## create model matrix and check model rank 
d1m = lapply(md, function(x){
  x =  x %>% mutate_if(is.character, as.factor) %$% time.group ; 
  x = model.matrix(~0 + x) ; 
  colnames(x) = str_sub(colnames(x), start = -9, end = -1)
  return(x)
})

# re QC model 
for (i in 1:length(d1m)) {
  model = d1m[[i]] ; print(i)
  stopifnot(Matrix::rankMatrix(model) == ncol(model)) ; stopifnot(any(colSums(model) == 0) == FALSE)
}
stopifnot(all.equal(
  lapply(md, rownames), lapply(e,colnames)
))
  
#  Confirm data has been normalized: 
# lapply(e, boxplot)
saveRDS(md, file = paste0(datapath, 'md.rds'))
saveRDS(e, file = paste0(datapath, 'e.rds'))

sessionInfo()
# R version 4.0.5 Patched (2021-03-31 r80136)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Mojave 10.14.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] magrittr_2.0.1  forcats_0.5.1   stringr_1.4.0   dplyr_1.0.4     purrr_0.3.4     readr_1.4.0     tidyr_1.1.2     tibble_3.0.6   
# [9] ggplot2_3.3.3   tidyverse_1.3.0 here_1.0.1     
# 
# loaded via a namespace (and not attached):
# [1] Rcpp_1.0.6        cellranger_1.1.0  pillar_1.4.7      compiler_4.0.5    dbplyr_2.1.0      tools_4.0.5       lattice_0.20-41  
# [8] jsonlite_1.7.2    lubridate_1.7.9.2 lifecycle_1.0.0   gtable_0.3.0      pkgconfig_2.0.3   rlang_0.4.10      Matrix_1.3-2     
# [15] reprex_1.0.0      cli_2.5.0         rstudioapi_0.13   DBI_1.1.1         haven_2.3.1       withr_2.4.1       xml2_1.3.2       
# [22] httr_1.4.2        fs_1.5.0          generics_0.1.0    vctrs_0.3.6       hms_1.0.0         rprojroot_2.0.2   grid_4.0.5       
# [29] tidyselect_1.1.0  glue_1.4.2        R6_2.5.0          readxl_1.3.1      modelr_0.1.8      backports_1.2.1   scales_1.1.1     
# [36] ellipsis_0.3.1    rvest_0.3.6       assertthat_0.2.1  colorspace_2.0-0  stringi_1.5.3     munsell_0.5.0     broom_0.7.5      
# [43] crayon_1.4.1     






