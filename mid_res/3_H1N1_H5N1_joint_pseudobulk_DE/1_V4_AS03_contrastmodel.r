# R version 4.0.5 
# H5 vs H1 day 1 DE contrast pseudobulk model 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(variancePartition))
#suppressMessages(library(scglmmr))
source("functions/analysis_functions.R")
source('functions/scglmmr.functions.R')

# make output directories 
datapath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/"); 
dir.create(datapath, recursive = TRUE)
figpath = here("mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/figuresV4/"); 
dir.create(figpath, recursive = TRUE)

# parallel options for dream lme4 fits 
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# read metadata to exract sample names in analysis 
samplemd = readRDS(file = here('data/samplemd.rds')) 
d1sx = samplemd %>% filter(time_cohort =='d1') %$% sample

# read processed pseudobulk data 
pb = readRDS(file = here('mid_res/variance_partition/generated_data/pb_vp.rds'))

# remove cell type string from sample names and subset to day 1 cohort
cnames = gsub("~.*","",colnames(pb[[1]]))
pb = lapply(pb, function(x){ 
  x %>% as.data.frame() %>% 
    setNames(nm = cnames) %>% 
    select(all_of(d1sx)) %>% 
    as.matrix() 
  })

# sample metadata for contrast model 
samplemd =readRDS(file = here('data/samplemd.rds')) 
samplemd = 
  samplemd %>% 
  filter(! time_cohort == 'd7') %>% 
  mutate(scaledage = (age - mean(age)) / sd(age)) %>% 
  rename('subjectid' = sampleid) %>% 
  rename('group' = adjmfc.group) %>%  
  mutate(group = ifelse(group %in% c('high', ' low'), "NOAS03", "AS03")) %>% 
  mutate(time.group = paste(timepoint, group, sep = "_"))
# get rid of space
samplemd$time.group = str_replace_all(
  string = samplemd$time.group,
  pattern = ' ',
  replacement = ''
)
# re-level time.group into ordered combined factor 
samplemd$time.group = factor(
  samplemd$time.group,
  levels = c('d0_AS03', 'd1_AS03', 'd0_NOAS03', 'd1_NOAS03')
  )
# format metadata 
samplemd = samplemd %>% 
  remove_rownames() %>% 
  column_to_rownames('sample')

# designmat 
met = samplemd[ ,c('gender', 'scaledage', 'time.group')]
mat = model.matrix( ~ 0 + time.group +  gender + scaledage, data = met)

################################
# specify random intercept model formula with combined factor 
f1 <- ~ 0 + time.group + gender + scaledage + (1|subjectid) 

# specify contrast matrix to test the fold change difference 
# based on levels of time.group this should be cmat = c(-1, 1, 1, -1, 0, 0)
L2 = makeContrastsDream(formula = f1, data = samplemd,
  contrasts = c(delta = "(time.groupd1_AS03 - time.groupd0_AS03) - (time.groupd1_NOAS03 - time.groupd0_NOAS03)")
  )
plotContrasts(L2) + ggsave(filename = paste0(figpath,'contrastmodel.pdf'), width = 7, height = 4)

# fit model on each subset 
# init store 
fit1 = v1 = list()
for (i in 1:length(pb)) {
  
  # init data 
  meta = samplemd
  form = f1 
  contrast_matrix = L2
  counts = pb[[i]]
  
  # dge list 
  d = edgeR::DGEList(counts = counts, samples = meta)
  
  # filter cell type specific lowly expressed genes and calc norm factors 
  gtable = edgeR::filterByExpr(y = d$counts, min.count = 3, design = mat)
  print(names(pb)[i]);print(table(gtable))
  d = d[gtable, keep.lib.sizes=FALSE]
  d = edgeR::calcNormFactors(object = d)
  
  
  # get voom observation level weights 
  v = voomWithDreamWeights(counts = d,
                           formula = form,
                           data = meta,
                           BPPARAM = pparam,
                           plot = TRUE, save.plot = TRUE)
  
  # fit contrast mixed model 
  fitmm = dream(exprObj = v, 
                formula = form,
                data = meta,
                L = contrast_matrix,
                BPPARAM = pparam, 
                useWeights = TRUE, REML = TRUE)
  
  # save results 
  v1[[i]] = v
  fit1[[i]] = fitmm
}
names(v1) = names(pb)


# save day 1 contrast fit 
saveRDS(object = fit1, file = paste0(datapath, 'fit12.rds'))
saveRDS(object = v1, file = paste0(datapath, 'v12.rds'))


# apply eBayes to model fit 
fit12 = readRDS(file = here('mid_res/3_H1N1_H5N1_joint_pseudobulk_DE/dataV4/fit12.rds'))

# over prev fit1e - compare 
fit12e  = lapply(fit12, variancePartition::eBayes) 
names(fit12e) = names(pb)
saveRDS(object = fit12e, file = paste0(datapath, 'fit12e.rds'))


# save model fitting data 
saveRDS(object = samplemd, file = paste0(datapath, 'samplemd12.rds'))
saveRDS(object = pb, file = paste0(datapath, 'pb12.rds'))
saveRDS(object = L2, file = paste0(datapath, 'L212.rds'))

# R 4.2.0 
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] viridis_0.6.3            viridisLite_0.4.2        simr_1.0.7               lme4_1.1-33             
# [5] Matrix_1.5-4             fgsea_1.24.0             ggpubr_0.6.0             magrittr_2.0.3          
# [9] variancePartition_1.28.9 BiocParallel_1.32.6      limma_3.54.2             here_1.0.1              
# [13] lubridate_1.9.2          forcats_1.0.0            stringr_1.5.0            dplyr_1.1.2             
# [17] purrr_1.0.1              readr_2.1.4              tidyr_1.3.0              tibble_3.2.1            
# [21] ggplot2_3.4.2            tidyverse_2.0.0         
# 
# loaded via a namespace (and not attached):
#   [1] utf8_1.2.3              spatstat.explore_3.2-1  reticulate_1.28         RUnit_0.4.32           
# [5] tidyselect_1.2.0        htmlwidgets_1.6.2       grid_4.2.0              Rtsne_0.16             
# [9] devtools_2.4.5          munsell_0.5.0           codetools_0.2-19        ragg_1.2.5             
# [13] ica_1.0-3               future_1.32.0           miniUI_0.1.1.1          withr_2.5.0            
# [17] spatstat.random_3.1-5   colorspace_2.1-0        progressr_0.13.0        Biobase_2.58.0         
# [21] Superpower_0.2.0        knitr_1.42              rstudioapi_0.14         Seurat_4.3.0           
# [25] ROCR_1.0-11             ggsignif_0.6.4          tensor_1.5              listenv_0.9.0          
# [29] emmeans_1.8.5           Rdpack_2.4              labeling_0.4.2          RLRsim_3.1-8           
# [33] polyclip_1.10-4         pheatmap_1.0.12         farver_2.1.1            rprojroot_2.0.3        
# [37] parallelly_1.35.0       vctrs_0.6.2             generics_0.1.3          afex_1.3-0             
# [41] clusterGeneration_1.3.7 xfun_0.39               timechange_0.2.0        R6_2.5.1               
# [45] doParallel_1.0.17       locfit_1.5-9.8          bitops_1.0-7            spatstat.utils_3.0-3   
# [49] cachem_1.0.8            promises_1.2.0.1        scales_1.2.1            nnet_7.3-19            
# [53] gtable_0.3.3            egg_0.4.5               globals_0.16.2          processx_3.8.1         
# [57] goftest_1.2-3           rlang_1.1.1             slanter_0.2-0           systemfonts_1.0.4      
# [61] splines_4.2.0           rstatix_0.7.2           lazyeval_0.2.2          checkmate_2.2.0        
# [65] spatstat.geom_3.2-1     broom_1.0.4             BiocManager_1.30.20     yaml_2.3.7             
# [69] reshape2_1.4.4          abind_1.4-5             backports_1.4.1         httpuv_1.6.10          
# [73] Hmisc_5.1-0             tools_4.2.0             usethis_2.2.0           gplots_3.1.3           
# [77] ellipsis_0.3.2          RColorBrewer_1.1-3      BiocGenerics_0.44.0     sessioninfo_1.2.2      
# [81] ggridges_0.5.4          Rcpp_1.0.10             plyr_1.8.8              base64enc_0.1-3        
# [85] progress_1.2.2          ps_1.7.5                prettyunits_1.1.1       rpart_4.1.19           
# [89] remaCor_0.0.11          deldir_1.0-6            pbapply_1.7-0           cowplot_1.1.1          
# [93] urlchecker_1.0.1        zoo_1.8-12              SeuratObject_4.1.3      ggrepel_0.9.3          
# [97] cluster_2.1.4           fs_1.6.2                data.table_1.14.8       scattermore_1.2        
# [101] lmerTest_3.1-3          lmtest_0.9-40           RANN_2.6.1              mvtnorm_1.1-3          
# [105] fitdistrplus_1.1-11     matrixStats_0.63.0      pkgload_1.3.2           hms_1.1.3              
# [109] patchwork_1.1.2         mime_0.12               evaluate_0.21           xtable_1.8-4           
# [113] pbkrtest_0.5.2          RhpcBLASctl_0.23-42     gridExtra_2.3           compiler_4.2.0         
# [117] KernSmooth_2.23-21      crayon_1.5.2            minqa_1.2.5             htmltools_0.5.5        
# [121] mgcv_1.8-42             later_1.3.1             tzdb_0.3.0              Formula_1.2-5          
# [125] snow_0.4-4              DBI_1.1.3               MASS_7.3-60             boot_1.3-28.1          
# [129] car_3.1-2               cli_3.6.1               rbibutils_2.2.13        parallel_4.2.0         
# [133] igraph_1.4.2            pkgconfig_2.0.3         numDeriv_2016.8-1.1     foreign_0.8-84         
# [137] sp_1.6-0                binom_1.1-1.1           plotly_4.10.2           spatstat.sparse_3.0-1  
# [141] foreach_1.5.2           estimability_1.4.1      callr_3.7.3             digest_0.6.31          
# [145] sctransform_0.3.5       RcppAnnoy_0.0.20        spatstat.data_3.0-1     rmarkdown_2.21         
# [149] leiden_0.4.3            fastmatch_1.1-3         htmlTable_2.4.1         edgeR_3.40.2           
# [153] uwot_0.1.14             curl_5.0.0              gtools_3.9.4            shiny_1.7.4            
# [157] nloptr_2.0.3            lifecycle_1.0.3         nlme_3.1-162            jsonlite_1.8.4         
# [161] aod_1.3.2               carData_3.0-5           desc_1.4.2              fansi_1.0.4            
# [165] pillar_1.9.0            ggsci_3.0.0             lattice_0.21-8          plotrix_3.8-2          
# [169] fastmap_1.1.1           httr_1.4.6              pkgbuild_1.4.1          survival_3.5-5         
# [173] glue_1.6.2              remotes_2.4.2           png_0.1-8               iterators_1.0.14       
# [177] stringi_1.7.12          profvis_0.3.8           textshaping_0.3.6       caTools_1.18.2         
# [181] memoise_2.0.1           irlba_2.3.5.1           future.apply_1.10.0   