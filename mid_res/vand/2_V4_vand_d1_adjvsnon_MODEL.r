# R version 4.0.5 
# vanderilt validation cohort of contrast model results from CITE-seq 
# analysis within analagous main immune cell lineage. 
suppressMessages(library(tidyverse))
suppressMessages(library(here))
suppressMessages(library(magrittr))
suppressMessages(library(variancePartition))
suppressMessages(library(scglmmr))

###### save paths  
datapath = here("mid_res/vand/generated_data/")

# parallel options 
register(SnowParam(4))
pparam = SnowParam(workers = 4, type = "SOCK", progressbar = TRUE)

# load data 
e = readRDS(file = here('mid_res/vand/generated_data/e.rds'))
md = readRDS(file = here('mid_res/vand/generated_data/md.rds'))

# process metadata for contrast model 
# order time.group factor 
md_lmer = lapply(md, function(x){
  x = x %>% mutate(time.group = factor(time.group, 
      levels = c("d0_AS03", "d1_AS03", "d0_xPBS", "d1_xPBS")))
  })

# specify formula for lme4 / dream and set up contrasts 
f1 <- ~ 0 + time.group + (1|subjectid)  

# specify contrast matrix to test the fold change difference 
# based on levels of time.group this should be cmat = c(-1, 1, 1, -1)
L2 = makeContrastsDream(
  formula = f1,
  data = md_lmer[[1]],
  contrasts = c(delta = "(time.groupd1_AS03 - time.groupd0_AS03) - (time.groupd1_xPBS - time.groupd0_xPBS)")
)
plotContrasts(L2) + ggsave(filename = paste0(figpath,'contrastmodel.pdf'), width = 7, height = 4)


# fit model on each subset 
fit1 = fitne = list()
for (i in 1:length(e)) {
  
  # init data 
  norm_dat = e[[i]]
  meta = md_lmer[[i]]
  form = f1 
  contrast_matrix = L2

    # fit contrast mixed model on prenormalized values 
  fitmm = dream(exprObj = norm_dat, 
                formula = form,
                data = meta,
                L = contrast_matrix,
                BPPARAM = pparam, 
                useWeights = FALSE, 
                REML = TRUE)
  # save results 
  fitne[[i]] = fitmm
  fit1[[i]] = variancePartition::eBayes(fitmm) 
}
names(fit1) = names(fitne) = names(e)

# Save 
# day 1 contrast fit 
saveRDS(object = fit1, file = paste0(datapath, 'fit1.rds'))
saveRDS(object = fitne, file = paste0(datapath, 'fit1.rds'))
# model fitting data 
saveRDS(object = md_lmer, file = paste0(datapath, 'md_lmer.rds'))
saveRDS(object = L2, file = paste0(datapath, 'L2.rds'))

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
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scglmmr_0.1.0            variancePartition_1.25.6 BiocParallel_1.24.1      limma_3.46.0             magrittr_2.0.1          
# [6] here_1.0.1               forcats_0.5.1            stringr_1.4.0            dplyr_1.0.4              purrr_0.3.4             
# [11] readr_1.4.0              tidyr_1.1.2              tibble_3.0.6             ggplot2_3.3.3            tidyverse_1.3.0         
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.1.0            lme4_1.1-26                 RSQLite_2.2.7               AnnotationDbi_1.52.0       
# [5] grid_4.0.5                  scatterpie_0.1.7            munsell_0.5.0               codetools_0.2-18           
# [9] statmod_1.4.35              withr_2.4.1                 colorspace_2.0-0            GOSemSim_2.16.1            
# [13] Biobase_2.50.0              rstudioapi_0.13             stats4_4.0.5                ggsignif_0.6.0             
# [17] DOSE_3.16.0                 labeling_0.4.2              MatrixGenerics_1.2.1        Rdpack_2.1.1               
# [21] emmeans_1.5.4               GenomeInfoDbData_1.2.4      polyclip_1.10-0             pheatmap_1.0.12            
# [25] bit64_4.0.5                 farver_2.0.3                rprojroot_2.0.2             downloader_0.4             
# [29] coda_0.19-4                 vctrs_0.3.6                 generics_0.1.0              TH.data_1.0-10             
# [33] R6_2.5.0                    doParallel_1.0.16           GenomeInfoDb_1.26.7         graphlayouts_0.7.2         
# [37] locfit_1.5-9.4              bitops_1.0-6                cachem_1.0.4                fgsea_1.16.0               
# [41] DelayedArray_0.16.3         assertthat_0.2.1            scales_1.1.1                multcomp_1.4-16            
# [45] ggraph_2.0.5                enrichplot_1.10.2           gtable_0.3.0                egg_0.4.5                  
# [49] tidygraph_1.2.0             sandwich_3.0-0              rlang_0.4.10                splines_4.0.5              
# [53] rstatix_0.7.0               broom_0.7.5                 BiocManager_1.30.10         reshape2_1.4.4             
# [57] abind_1.4-5                 modelr_0.1.8                backports_1.2.1             qvalue_2.22.0              
# [61] clusterProfiler_3.18.1      tools_4.0.5                 ellipsis_0.3.1              gplots_3.1.1               
# [65] RColorBrewer_1.1-2          BiocGenerics_0.36.1         Rcpp_1.0.6                  plyr_1.8.6                 
# [69] progress_1.2.2              zlibbioc_1.36.0             RCurl_1.98-1.3              prettyunits_1.1.1          
# [73] ggpubr_0.4.0                viridis_0.5.1               cowplot_1.1.1               S4Vectors_0.28.1           
# [77] zoo_1.8-8                   SummarizedExperiment_1.20.0 haven_2.3.1                 ggrepel_0.9.1              
# [81] fs_1.5.0                    data.table_1.14.0           lmerTest_3.1-3              DO.db_2.9                  
# [85] openxlsx_4.2.3              reprex_1.0.0                mvtnorm_1.1-1               matrixStats_0.58.0         
# [89] hms_1.0.0                   GSVA_1.38.2                 xtable_1.8-4                pbkrtest_0.5-0.1           
# [93] RhpcBLASctl_0.21-247.1      XML_3.99-0.6                rio_0.5.16                  readxl_1.3.1               
# [97] IRanges_2.24.1              gridExtra_2.3               compiler_4.0.5              KernSmooth_2.23-18         
# [101] crayon_1.4.1                shadowtext_0.0.9            minqa_1.2.4                 ggfun_0.0.4                
# [105] snow_0.4-3                  lubridate_1.7.9.2           DBI_1.1.1                   tweenr_1.0.2               
# [109] dbplyr_2.1.0                MASS_7.3-53.1               boot_1.3-27                 Matrix_1.3-2               
# [113] car_3.0-10                  cli_2.5.0                   rbibutils_2.0               parallel_4.0.5             
# [117] igraph_1.2.6                GenomicRanges_1.42.0        pkgconfig_2.0.3             rvcheck_0.1.8              
# [121] numDeriv_2016.8-1.1         foreign_0.8-81              xml2_1.3.2                  foreach_1.5.1              
# [125] annotate_1.68.0             XVector_0.30.0              estimability_1.3            rvest_0.3.6                
# [129] digest_0.6.27               graph_1.68.0                cellranger_1.1.0            fastmatch_1.1-0            
# [133] edgeR_3.32.1                GSEABase_1.52.1             curl_4.3                    gtools_3.8.2               
# [137] nloptr_1.2.2.2              lifecycle_1.0.0             nlme_3.1-152                jsonlite_1.7.2             
# [141] aod_1.3.1                   carData_3.0-4               viridisLite_0.3.0           pillar_1.4.7               
# [145] lattice_0.20-41             fastmap_1.1.0               httr_1.4.2                  survival_3.2-10            
# [149] GO.db_3.12.1                glue_1.4.2                  zip_2.1.1                   iterators_1.0.13           
# [153] bit_4.0.4                   ggforce_0.3.3               stringi_1.5.3               blob_1.2.1                 
# [157] org.Hs.eg.db_3.12.0         caTools_1.18.1              memoise_2.0.0     


